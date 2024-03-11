#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <deque>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <limits.h>
#include <thread>
#include <sys/stat.h>
#include <condition_variable>
#include <mutex>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/gzreader.h>
#include "common.h"
#include "species_kmers.h"
#include "kmsuftree.h"

KSEQ_INIT(gzFile, gzread);

using std::cout;
using std::endl;
using namespace std;

// ===== species_kmers.cpp
// This file contains functions used to scan reads for species-specific k-mers
// Used by demux_species.

species_kmer_counter::species_kmer_counter(int nt, 
    int k, 
    int ns,
    bc_whitelist* wl,
    robin_hood::unordered_map<unsigned long, map<short, int> >* bsc,
    int umi_start,
    int umi_len){
    
    this->terminate_threads = false;
    this->num_threads = nt;
    this->k = k;
    this->wl = wl;
    this->num_species = ns;
    this->initialized = false;
    this->on = false;
    this->bc_species_counts = bsc;
    this->umi_start = umi_start;
    this->umi_len = umi_len;
}

species_kmer_counter::~species_kmer_counter(){
    if (this->on){
        this->close_pool();
        this->on = false;
    }
    if (this->initialized){
        kmsuftree_destruct(kt, 0);
    }
    // Free all UMI counters
    for (robin_hood::unordered_map<unsigned long, umi_set*>::iterator x = bc_species_umis.begin();
        x != bc_species_umis.end(); ){
        delete x->second;
        x = bc_species_umis.erase(x);
    }
}

void species_kmer_counter::init(short species_idx, string& kmerfile){
    this->terminate_threads = false;

    this_species = species_idx;

    if (this->initialized){
        kmsuftree_destruct(kt, 0);
    }
    kt = init_kmer_tree(k);
    parse_kmer_counts_serial(kmerfile, species_idx);
    
    this->initialized = true;
}

/**
 * Only using canonical k-mers = first in sort order relative to reverse complement.
 * Check if a given k-mer or its reverse complement comes first in alphabetical sort order.
 */
bool species_kmer_counter::is_rc_first(char* kmer, int k){
    for (int i = 0; i < k; ++i){
        char char_f = kmer[k];
        switch(kmer[k]){
            case 'a':
            case 'A':
                char_f = 'A';
                break;
            case 'c':
            case 'C':
                char_f = 'C';
                break;
            case 'g':
            case 'G':
                char_f = 'G';
                break;
            case 't':
            case 'T':
                char_f = 'T';
                break;
            default:
                char_f = 'N';
                break;
        }
        char char_rc = 'N';
        switch(kmer[k-1-i]){
            case 'A':
            case 'a':
                char_rc = 'T';
                break;
            case 'c':
            case 'C':
                char_rc = 'G';
                break;
            case 'g':
            case 'G':
                char_rc = 'C';
                break;
            case 't':
            case 'T':
                char_rc = 'A';
                break;
            default:
                char_rc = 'N';
                break;
        }
        if (char_f == 'N' || char_rc == 'N'){
            // No matter; neither will be in the lookup table.
            return false;
        }
        if (char_rc < char_f){
            return true;
        }
        else if (char_f < char_rc){
            return false;
        }   
    }
    return false;
}

/**
 * Adds sequences from reads to data that will be retrieved by worker threads.
 */
void species_kmer_counter::process_gex_files(string& r1filename, 
    string& r2filename){

    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    f_fp = gzopen(r1filename.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1filename.c_str());
        exit(1);
    }    
    r_fp = gzopen(r2filename.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2filename.c_str());
        exit(1);
    }
   
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R3 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        rp_info params;
        params.seq_f = seq_f->seq.s;
        params.seq_r = seq_r->seq.s;
        
        add_rp_job(params);        
    
    }
}

// Stop all running threads and then destroy them
void species_kmer_counter::close_pool(){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->terminate_threads = true;   
    }
    this->has_jobs.notify_all();
    for (int i = 0; i < this->num_threads; ++i){
        this->threads[i].join();
    }
    this->threads.clear();
    this->on = false;

    // Free all UMI counters
    for (robin_hood::unordered_map<unsigned long, umi_set*>::iterator x = bc_species_umis.begin();
        x != bc_species_umis.end(); ){
        delete x->second;
        x = bc_species_umis.erase(x);
    }
}

/**
 * Add to the queue of jobs to read paired reads
 */
void species_kmer_counter::add_rp_job(rp_info& info){
    {
        // Add job to processing queue
        unique_lock<mutex> lock(this->queue_mutex);
        this->rp_jobs.push_back(info);
    }
    this->has_jobs.notify_one();
}

/**
 * Add to the queue of jobs to read read triplets
 */
void species_kmer_counter::add_rt_job(rt_info& info){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->rt_jobs.push_back(info);
    }
    this->has_jobs.notify_one();
}

// Count k-mers for one species in a specific read.
int species_kmer_counter::scan_seq_kmers(string& seq){
    
    int species_count = 0;
    
    // Only one species per job.
    short this_species;

    for (int i = 0; i < seq.length() - kt->k; ++i){
        int rc = 0;
        if (species_kmer_counter::is_rc_first(&seq[i], kt->k)){
            rc = 1;
        }
        void* dat = kmer_tree_lookup(&seq[i], kt, rc);
        if (dat != NULL){
            //this_species = *(short*)dat;
            ++species_count;
        }
    }
    
    return species_count;
}

void species_kmer_counter::scan_gex_data(rp_info& params){
    
    unsigned long bc_key;
    if (wl->lookup(params.seq_f.c_str(), bc_key)){
        
        bool dup_read = false;
        if (umi_start >= 0 && umi_len > 0){
            unique_lock<mutex> lock(this->umi_mutex);
            umi u(params.seq_f.c_str() + umi_start, umi_len);
            if (bc_species_umis.count(bc_key) == 0){
                umi_set* s = new umi_set(umi_len);
                bc_species_umis.emplace(bc_key, s);
            }
            if (bc_species_umis[bc_key]->add(u)){
                dup_read = true;
            }
        }
        
        if (dup_read){
            return;
        }

        // In 10x scRNA-seq data, only the reverse read contains information
        // forward read is barcode
        int nk = scan_seq_kmers(params.seq_r);

        if (nk > 0){    
            unique_lock<mutex> lock(this->counts_mutex);
            if (bc_species_counts->count(bc_key) == 0){
                map<short, int> m;
                bc_species_counts->emplace(bc_key, m);
            } 
            if ((*bc_species_counts)[bc_key].count(this_species) == 0){
                (*bc_species_counts)[bc_key].insert(make_pair(this_species, 0));
            }
            (*bc_species_counts)[bc_key][this_species]++;
        }
    }
}

void species_kmer_counter::parse_kmer_counts_serial(string& countsfilename,
    short species_idx){
    
    gzreader reader(countsfilename);
    while (reader.next()){
        short* si = (short*) malloc(sizeof(short));
        *si = species_idx;
        kmer_tree_add(reader.line, kt, (void*)si, 0);
    } 
    /*
    ifstream infile(countsfilename);
    string line;
    while (infile >> line){
        short* si = (short*) malloc(sizeof(short));
        *si = species_idx;
        kmer_tree_add(line.c_str(), kt, (void*)si, 0); 
    }
    */
}

/**
 * Worker function for a GEX read scanning job
 */
void species_kmer_counter::gex_thread(){
    
     while(true){
        if (this->rp_jobs.size() == 0 && this->terminate_threads){
            return;
        }
        rp_info params;
        {
            unique_lock<mutex> lock(this->queue_mutex);
            this->has_jobs.wait(lock, [this]{ return rp_jobs.size() > 0 ||
                terminate_threads;});
            if (this->rp_jobs.size() == 0 && this->terminate_threads){
                return;
            }
            params = this->rp_jobs[0];
            this->rp_jobs.pop_front();
        }
        scan_gex_data(params);
     }
 }

void species_kmer_counter::launch_gex_threads(){
    
    if (!this->initialized){
        fprintf(stderr, "ERROR: not initialized\n");
        exit(1);
    }

    this->terminate_threads = false;

    for (int i = 0; i < this->num_threads; ++i){
        
        this->threads.push_back(thread(&species_kmer_counter::gex_thread,
            this));
    }
    this->on = true;
}


