#ifndef _CELLID_SPECIES_KMERS_H
#define _CELLID_SPECIES_KMERS_H
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
#include "common.h"
#include "kmsuftree.h"

// ===== kmers_in_reads.h
// Contains functions used to scan reads for species-specific kmers
// used by demux_species.

// Information to represent read pairs.
struct rp_info{
    std::string seq_f;
    std::string seq_r;
};

// Information to represent read triplets.
struct rt_info{
    std::string seq_1;
    std::string seq_2;
    std::string seq_3;
};

// Define a class for parallel processing stuff
class species_kmer_counter{
    protected:
        static bool is_rc_first(char* kmer, int k);
    private:
        
        // Common to all jobs
        std::mutex queue_mutex;
        std::condition_variable has_jobs;
        bool terminate_threads;
        int num_threads;
        std::vector<std::thread> threads;
        bool initialized;
        bool on;
        std::string out_base;
        int num_species;
        short this_species;
        kmer_tree_p kt;
        int k;
        
        int umi_start;
        int umi_len;

        std::mutex counts_mutex;
        std::mutex umi_mutex;
        robin_hood::unordered_map<unsigned long, std::map<short, int> >* bc_species_counts;
        robin_hood::unordered_map<unsigned long, umi_set* > bc_species_umis;
        bc_whitelist* wl;

        // Task-specific
        // Parameter queue for kmer file parsing jobs
        std::deque<std::pair<std::string, short> > kmer_parse_jobs;

        // Parameter queue for paired-read (GEX) jobs
        std::deque<rp_info> rp_jobs;
        
        // Parameter queue for three-read (ATAC) jobs
        std::deque<rt_info> rt_jobs;
        
        // Function to process RNA-seq reads
        void gex_thread();
        
        int scan_seq_kmers(std::string& seq);
        
        void scan_gex_data(rp_info& info);
        
        void add_rp_job(rp_info& info);

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

        void close_pool();
         
        void launch_gex_threads();
        
        void process_gex_files(std::string& r1filename, std::string& r2filename);
        
};

#endif
