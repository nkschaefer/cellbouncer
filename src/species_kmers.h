#ifndef _CELLBOUNCER_SPECIES_KMERS_H
#define _CELLBOUNCER_SPECIES_KMERS_H
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
#include <condition_variable>
#include <mutex>
#include <sys/stat.h>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/khashtable.h>
#include "common.h"

// ===== species_kmers.h
// Contains functions used to scan reads for species-specific kmers
// used by demux_species.

// Information to represent read pairs.
struct rp_info{
    //std::string seq_f;
    //std::string seq_r;
    //char[FQ_READ_LEN] seq_f;
    //char[FQ_READ_LEN] = seq_r;
    char* seq_f;
    char* seq_r;
    int len_f;
    int len_r;
    rp_info(const char* seq_f, int seq_f_len, const char* seq_r, int seq_r_len); 
    rp_info(const rp_info& other);
    ~rp_info();    
};

// Information to represent read triplets.
struct rt_info{
    std::string seq_1;
    std::string seq_2;
    std::string seq_3;
};
/*
struct kmer_node_ptr{
    kmer_node_ptr* f_A;
    kmer_node_ptr* f_C;
    kmer_node_ptr* f_G;
    kmer_node_ptr* f_T;
    kmer_node_ptr* r_A;
    kmer_node_ptr* r_C;
    kmer_node_ptr* r_G;
    kmer_node_ptr* r_T;
    
    bool f_A_flip;
    bool f_C_flip;
    bool f_G_flip;
    bool f_T_flip;
    bool r_A_flip;
    bool r_C_flip;
    bool r_G_flip;
    bool r_T_flip;
    
    short species;
};
*/

// Define a class for parallel processing stuff
class species_kmer_counter{
    protected:
        static bool is_rc_first(const char* kmer, int k);
    private:
        
        // Common to all jobs
        std::mutex queue_mutex;
        std::condition_variable has_jobs;
        bool terminate_threads;
        int num_threads;
        std::vector<std::thread> threads;
        std::vector<std::vector<int> > species_counts;
        std::deque<khashkey> khashkeys;

        bool initialized;
        bool on;
        std::string out_base;
        int num_species;
        khashtable<short> tab;
        int k;
        
        char* kmer_buf;
        bool kmer_buf_init;
         
        // -1 (all k-mers) or a number to sample from each species
        int n_samp;
                
        int umi_start;
        int umi_len;

        std::mutex counts_mutex;
        std::mutex umi_mutex;
        robin_hood::unordered_map<unsigned long, std::map<short, int> >* bc_species_counts;
        
        robin_hood::unordered_map<unsigned long, umi_set_exact* > bc_species_umis;

        bool use_umis;
        bc_whitelist* wl;
        
        void close_pool();
         
        void launch_gex_threads();
        
        // Task-specific
        // Parameter queue for kmer file parsing jobs
        std::deque<std::pair<std::string, short> > kmer_parse_jobs;

        // Parameter queue for paired-read (GEX) jobs
        std::deque<rp_info> rp_jobs;
        
        // Parameter queue for three-read (ATAC) jobs
        std::deque<rt_info> rt_jobs;
        
        //void add_node_to_tree_simple(char* kmer, char* buf, kmer_tree_p kt, short species_idx);
        //void add_node_to_tree_canonical(char* kmer, char* buf, kmer_tree_p kt, short species_idx);

        // Function to process RNA-seq reads
        void gex_thread(int thread_idx);
        
        void scan_seq_kmers(const char* seq, int len, int* species_counts, khashkey& key);
        
        void scan_gex_data(const char* seq_f, int seq_f_len, const char* seq_r, int seq_r_len, int thread_idx=0);
        
        void add_rp_job(const char* seq_f, int seq_f_len, const char* seq_r, int seq_r_len);

        void add_rt_job(rt_info& info);

        void parse_kmer_counts_serial(std::string& filename,
            short species_idx);

   
    public:
        
        species_kmer_counter(int nt, int k, int ns,
            bc_whitelist* wl, 
            robin_hood::unordered_map<unsigned long, map<short, int> >* bsc,
            int umi_start = 16, int umi_len = 12);

        ~species_kmer_counter();

        void init(short species_idx, std::string& kmerfile);
        void add(short species_idx, std::string& kmerfile);

        void disable_umis(); 
        void enable_umis();
        
        void set_n_samp(int ns);

        void process_gex_files(std::string& r1filename, std::string& r2filename);
        
};

#endif
