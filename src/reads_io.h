#ifndef READSIO_H
#define READSIO_H
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
#include "bc.h"
#include "common.h"
#include "kmsuftree.h"
#include "serialize.h"
#include "robin_hood.h"

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
class readsThreadPool{
    
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
        //std::vector<gzFile> out_files;
        outstream_info out_files[100];
        std::vector<std::string> out_filenames;

        // Task-specific
        // Parameter queue for kmer file parsing jobs
        std::deque<std::pair<std::string, short> > kmer_parse_jobs;

        // Parameter queue for paired-read (GEX) jobs
        std::deque<rp_info> rp_jobs;
        // Parameter queue for three-read (ATAC) jobs
        std::deque<rt_info> rt_jobs;
        
    public:
        readsThreadPool(int nt, int ns, std::string& ob);

        ~readsThreadPool();

        void init();

        void close_pool();
         
        void load_counts_from_tmp(
            robin_hood::unordered_map<unsigned long, std::map<short, double> >& bc_species_counts,
            robin_hood::unordered_map<unsigned long, double>& bc_tots);
     
        void add_kmer_parse_job(std::string& fn, short species);

        void add_rp_job(rp_info& info);

        void add_rt_job(rt_info& info);
        
        // Function to parse a kmer file
        void kmer_parse_thread(kmer_tree_p kt, 
            std::recursive_mutex& ktmutex);

        // Function to process RNA-seq reads
        void gex_thread(kmer_tree_p kt, 
            outstream_info& os,
            //std::map<unsigned long, std::map<short, float> >& bc_species_counts,
            //std::map<unsigned long, int>& bc_species_tots,
            //std::recursive_mutex& counts_mutex,
            //std::recursive_mutex& tots_mutex,
            bcset& rna_bc2idx,
            std::vector<unsigned long>& whitelist_rna);

        // Function to process ATAC-seq reads
        void atac_thread(kmer_tree_p kt,
            outstream_info& os,
            //std::map<unsigned long, std::map<short, float> >& bc_species_counts,
            //std::map<unsigned long, int>& bc_species_tots,
            //std::recursive_mutex& counts_mutex,
            //std::recursive_mutex& tots_mutex,
            bcset& atac_bc2idx,
            std::vector<unsigned long>& whitelist_rna);

        void launch_kmer_parse_threads(kmer_tree_p kt,
            std::recursive_mutex& ktmutex);

        void launch_gex_threads(kmer_tree_p kt,
            //std::map<unsigned long, std::map<short, float> >& bc_species_counts,
            //std::map<unsigned long, int>& bc_species_tots,
            //std::recursive_mutex& counts_mutex,
            //std::recursive_mutex& tots_mutex,
            bcset& rna_bc2idx,
            std::vector<unsigned long>& whitelist_rna);
        
        void launch_atac_threads(kmer_tree_p kt,
            //std::map<unsigned long, std::map<short, float> >& bc_species_counts,
            //std::map<unsigned long, int>& bc_species_tots,
            //std::recursive_mutex& counts_mutex,
            //std::recursive_mutex& tots_mutex,
            bcset& atac_bc2idx,
            std::vector<unsigned long>& whitelist_rna);

};

// Function to trim the path off of a filename
std::string filename_nopath(std::string& filename);

void parse_whitelists(std::string& atac_filename,
    std::string& gex_filename,
    bcset& atac_bc2idx,
    std::vector<unsigned long>& whitelist_rna,
    bcset& rna_bc2idx);

void process_gex_files(std::string& r1filename,
    std::string& r2filename,
    readsThreadPool& pool);

void demux_atac_reads(std::string& r1filename, 
    std::string& r2filename, 
    std::string& r3filename, 
    bcset& bc2species,
    std::map<short, std::string>& idx2species,
    std::string& outdir,
    bcset& atac_bc2idx,
    std::vector<unsigned long>& whitelist_rna,
    long int& atac_reads_tot,
    std::map<short, long int>& atac_reads_pass);
 
void demux_rna_reads(std::string& r1filename, 
    std::string& r2filename,
    std::string& file_prefix,
    bcset& bc2species,
    std::map<short, std::string>& idx2species,
    std::string& outdir,
    long int& rna_reads_tot,
    std::map<short, long int>& rna_reads_pass);
/* 
// Function to scan GEX reads (in a multithreaded way) for species-specific k-mers
void scan_gex_reads(std::string& r1filename, 
    std::string& r2filename,
    bcset& rna_bc2idx,
    std::vector<unsigned long>& whitelist_rna,
    int bclen,
    int k,
    kmer_tree_p kt,
    std::map<unsigned long, std::map<short, float> >& bc_species_counts,
    std::map<unsigned long, int>& bc_species_tots);
*/
void parse_kmer_counts_serial(std::string& filename,
    short species_idx,
    kmer_tree_p kt);

#endif
