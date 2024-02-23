#include <getopt.h>
#include <argp.h>
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
#include <sys/stat.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <limits.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <htswrapper/bam.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include "demux_species_io.h"
#include "reads_demux.h"
#include "species_kmers.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_species [OPTIONS]\n");
    fprintf(stderr, "Given reads from a multi-species pooled experiment and lists\n");
    fprintf(stderr, "   of species-specific k-mers, demultiplexes the reads by species\n");
    fprintf(stderr, "   and creates library files for running cellranger-arc\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n   ===== GENERAL OPTIONS =====\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    fprintf(stderr, "   --doublet_rate -D What is the prior expected doublet rate?\n");
    fprintf(stderr, "       (OPTIONAL; default = 0.1). Must be a decimal between 0 and 1,\n");
    fprintf(stderr, "       exclusive.\n");
    fprintf(stderr, "   --output_directory -o The directory in which to place output files.\n");
    fprintf(stderr, "       Read file names will be extracted from read group tags.\n");
    fprintf(stderr, "   --num_threads -t The number of threads to use for parallel processing\n");
    fprintf(stderr, "       (default 1)\n");
    fprintf(stderr, "   --dump -d Only dump per-barcode data (barcode, then count of reads\n");
    fprintf(stderr, "       per species (tab separated)) and barcode-to-species assignments\n");
    fprintf(stderr, "       instead of demultiplexing reads. These files are created in\n");
    fprintf(stderr, "       <output_directory>/species_counts.txt and\n");
    fprintf(stderr, "       <output_directory>/species_assignments.txt regardless.\n");
    fprintf(stderr, "       Standard behavior is to create this file and then demultiplex\n");
    fprintf(stderr, "       reads; this option causes the program to quit after generating\n");
    fprintf(stderr, "       the file.\n");
    fprintf(stderr, "\n   ===== READ FILE INPUT OPTIONS =====\n");
    fprintf(stderr, "   --atac_r1 -1 ATAC R1 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "   --atac_r2 -2 ATAC R2 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "   --atac_r3 -3 ATAC R3 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "   --rna_r1 -r Forward RNA-seq reads to demultiplex (can specify multiple\n");
    fprintf(stderr, "       times)\n");
    fprintf(stderr, "   --rna_r2 -R Reverse RNA-seq reads to demultiplex (can specify multiple\n");
    fprintf(stderr, "       times)\n");
    fprintf(stderr, "   --custom_r1 -x Forward other (i.e. sgRNA or antibody capture) reads\n");
    fprintf(stderr, "       to demultiplex (can specify multiple times). Assumes barcodes are\n");
    fprintf(stderr, "       at the beginning of R1.\n");
    fprintf(stderr, "   --custom_r2 -X Reverse other (i.e. sgRNA or antibody capture) reads\n");
    fprintf(stderr, "       to demultiplex (can specify multiple times). Assumes barcodes are\n");
    fprintf(stderr, "       at the beginning of R1.\n");
    fprintf(stderr, "   --names_custom -n Name of data type in custom reads file, in same\n");
    fprintf(stderr, "       number and order as read files. Presets: CRISPR = CRISPR sgRNA\n");
    fprintf(stderr, "       capture, Ab = antibody capture. Names will be appended to the\n");
    fprintf(stderr, "       beginning of demultiplexed FASTQ files and inserted into 10X\n");
    fprintf(stderr, "       library files. For example, if providing sgRNA capture files\n");
    fprintf(stderr, "       sgRNA_R1.fq.gz and sgRNA_R2.fq.gz along with antibody capture\n");
    fprintf(stderr, "       files anti_R1.fq.gz and anti_R2.fq.gz, you could specify:\n");
    fprintf(stderr, "       -x sgRNA_R1.fq.gz -X sgRNA_R2.fq.gz\n");
    fprintf(stderr, "       -x anti_R1.fq.gz -X anti_R2.fq.gz\n");
    fprintf(stderr, "       -n CRISPR -n Antibody.\n");
    fprintf(stderr, "\n   ===== BARCODE WHITELIST OPTIONS =====\n");
    fprintf(stderr, "   --whitelist_atac -w If multiome data and demultiplexing ATAC-seq\n");
    fprintf(stderr, "       reads, provide both the ATAC-seq barcode whitelist (here) and the\n");
    fprintf(stderr, "       RNA-seq barcode whitelist (-W) (REQUIRED). If not multiome or\n");
    fprintf(stderr, "       RNA-seq only, this whitelist is not required.\n");
    fprintf(stderr, "   --whitelist_rna -W If multiome data and demultiplexing ATAC-seq reads,\n");
    fprintf(stderr, "       provide both the ATAC-seq barcode whitelist (-w) and the RNA-seq\n");
    fprintf(stderr, "       barcode whitelist (here) (REQUIRED). If not multiome or RNA-seq only,\n");
    fprintf(stderr, "       these are not required.\n");
    fprintf(stderr, "\n   ===== OTHER INPUT OPTIONS =====\n");
    fprintf(stderr, "   --k -k File containing species-specific k-mers. Specify once per species.\n");
    fprintf(stderr, "   --species -s Name of species corresponding to entries of -k. Specify\n");
    fprintf(stderr, "       multiple times, the same number and order as -k.\n");
    fprintf(stderr, "   --countsfile -C If a previous run already completed and you want to run\n");
    fprintf(stderr, "       again using the same counts, specify the counts file\n");
    fprintf(stderr, "       (<output_directory>/species_counts.txt) here. This will skip the\n");
    fprintf(stderr, "       expensive step of iterating through the BAM file to count reads.\n");
    fprintf(stderr, "       You can then either demultiplex reads (default) or dump data and quit\n");
    fprintf(stderr, "       (-d). If using this parameter, then --speciesfile/-S is also required.\n");
    fprintf(stderr, "   --speciesfile -S If loading data from a previous run with -C, then -S is\n");
    fprintf(stderr, "       also required. This parameter is a tab separated file that lists index,\n");
    fprintf(stderr, "       followed by human-readable species name. It will be generated by\n");
    fprintf(stderr, "       previous runs as <output_directory>/species_names.txt.\n");
    fprintf(stderr, "\n ===== NOTES =====\n");
    fprintf(stderr, "   This program works by counting k-mers in RNA-seq data exclusively. The other\n");
    fprintf(stderr, "   types of reads are provided to be demultiplexed only, by sharing of barcodes\n");
    fprintf(stderr, "   with the RNA-seq data. If you provide other types of data (i.e. ATAC, sgRNA\n");
    fprintf(stderr, "   capture), this program will attempt to create a 10X Genomics-format library\n");
    fprintf(stderr, "   file to help run data. Any feature barcoding data will need an accompanying\n");
    fprintf(stderr, "   feature reference file, though, which must be created manually (see\n");
    fprintf(stderr, "   https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).\n");
    fprintf(stderr, "   Once k-mers are counted once, it creates a counts file and species names\n");
    fprintf(stderr, "   file in the output directory. These can be loaded to demultiplex reads\n");
    fprintf(stderr, "   instead of repeating k-mer counting, which is the most expensive step.\n");
    fprintf(stderr, "   When counting k-mers, species-specific k-mer files (-k), species names\n");
    fprintf(stderr, "   (-s), RNA-seq reads (-r and -R), and an RNA-seq barcode whitelist (-W) are\n");
    fprintf(stderr, "   all required.\n");
    fprintf(stderr, "   When demultiplexing based on previously-computed k-mer counts, a counts file\n");
    fprintf(stderr, "   (-C), species names file (-S), and all reads to demultiplex\n");
    fprintf(stderr, "   (-r/-R, -1/-2/-3, -x/-X/-n) are required.\n");
    exit(code);
}

/**
 * Use a mixture of multinomial distributions to model k-mer counts from each species
 * in each cell. Each species of origin will be a source distribution, as will
 * each possible doublet combination.
 *
 * Posterior log likelihood ratios of each species identity come from fit distribution
 * parameters and the provided prior estimate of doublet rate.
 *
 * Cell assignments are placed into the given maps - one for singlets, and another for
 * doublets. Log likelihood ratios of best to second best assignment are also stored.
 */ 
void fit_model(robin_hood::unordered_map<unsigned long, map<short, int> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, short>& bc2species, 
    robin_hood::unordered_map<unsigned long, pair<unsigned int, unsigned int> >& bc2doublet,
    robin_hood::unordered_map<unsigned long, double>& bc2llr,
    map<short, string>& idx2species,
    double doublet_rate){
    
    int n_species = idx2species.size();
    
    vector<mixtureDist> dists;
    
    // Which component distributions represent doublets? 
    map<int, bool> dist_doublet;
    // What are the idenities of component distributions representing doublets?
    map<int, pair<int, int> > dist2doublet_comb;
    // What are the identities of component distributions representing singlets?
    map<int, int> dist2singlet;

    int doublet_dist_count = 0;
    int singlet_dist_count = 0;

    // How will we initialize multinomial parameters for components that match the expected species? 
    double target_weight = 0.9;

    // Each species will be represented by a multinomial distribution, with one component for
    // each species. 
    for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); ++i2s){
        
        vector<double> params;
        for (int i = 0; i < n_species; ++i){
            if (i == i2s->first){
                // Target species
                params.push_back(target_weight);
            }
            else{
                // Non-target species
                params.push_back((1.0-target_weight)/((double)(n_species-1)));
            }
        }
        
        mixtureDist dist("multinomial", params);
        dist.name = i2s->second;
        dist.set_num_inputs(n_species);
        dist_doublet.insert(make_pair(dists.size(), false));
        dist2singlet.insert(make_pair(dists.size(), i2s->first));
        dists.push_back(dist);
        singlet_dist_count++;
        
        
        // Create doublet distribution
        for (int j = i2s->first+1; j < n_species; ++j){
            vector<double> doublet_params;
            for (int i = 0; i < n_species; ++i){
                if (i == i2s->first || i == j){
                    // Target species
                    doublet_params.push_back(target_weight/2.0);
                }
                else{
                    // Non-target species
                    doublet_params.push_back((1.0-target_weight)/((double)(n_species-2)));
                }
            }
            mixtureDist dist("multinomial", doublet_params);
            string dname;
            if (i2s->second < idx2species[j]){
                dname = i2s->second + "." + idx2species[j];
            }
            else{
                dname = idx2species[j] + "." + i2s->second;
            }
            dist.name = dname;
            dist.set_num_inputs(n_species);
            dist_doublet.insert(make_pair(dists.size(), true));
            dist2doublet_comb.insert(make_pair(dists.size(), make_pair(i2s->first, j)));
            dists.push_back(dist);
            doublet_dist_count++;
        }
        
    }
    
    // Prepare input data
    vector<vector<double> > obs;
    vector<unsigned long> bcs;
    for (robin_hood::unordered_map<unsigned long, map<short, int> >::iterator x = 
        bc_species_counts.begin(); x != bc_species_counts.end(); ++x){
        bcs.push_back(x->first);
        vector<double> row;
        for (map<short, int>::iterator sc = x->second.begin(); sc != x->second.end(); ++sc){
            row.push_back((double)sc->second);
        }
        obs.push_back(row);
    }

    vector<double> dist_weights;
    for (map<int, bool>::iterator dd = dist_doublet.begin(); dd != dist_doublet.end();
        ++dd){
        double weight;
        if (dd->second){
            weight = doublet_rate / (double)doublet_dist_count;
        }
        else{
            weight = (1.0-doublet_rate)/(double)singlet_dist_count;
        }
        dist_weights.push_back(weight);
    }
    
    // Create and fit mixture model 
    mixtureModel mod(dists, dist_weights); 
    mod.fit(obs);
    mod.print();

    // Use fit distributions and doublet rate prior to assign identities
    // and log likelihood ratios to cell barcodes

    for (int i = 0; i < obs.size(); ++i){
        
        vector<pair<double, int> > lls;
        for (int j = 0; j < mod.n_components; ++j){
            double ll = mod.dists[j].loglik(obs[i]);
            // Incorporate prior prob of doublet rate
            if (dist_doublet[j]){
                ll += log2(doublet_rate);
            }
            else{
                ll += log2(1.0 - doublet_rate);
            }
            lls.push_back(make_pair(-ll, j));
        }
        sort(lls.begin(), lls.end());

        double llr = -lls[0].first - -lls[1].first;
        int maxmod = lls[0].second;
        bool is_doublet = dist_doublet[maxmod];
        if (is_doublet){
            bc2doublet.emplace(bcs[i], dist2doublet_comb[maxmod]);            
        }
        else{
            bc2species.emplace(bcs[i], dist2singlet[maxmod]);
        }
        bc2llr.emplace(bcs[i], llr);
    }
}

int main(int argc, char *argv[]) {    
    
    // Define long-form program options 
    static struct option long_options[] = {
       {"output_directory", required_argument, 0, 'o'},
       {"atac_r1", required_argument, 0, '1'},
       {"atac_r2", required_argument, 0, '2'},
       {"atac_r3", required_argument, 0, '3'},
       {"rna_r1", required_argument, 0, 'r'},
       {"rna_r2", required_argument, 0, 'R'},
       {"custom_r1", required_argument, 0, 'x'},
       {"custom_r2", required_argument, 0, 'X'},
       {"names_custom", required_argument, 0, 'n'},
       {"whitelist_atac", required_argument, 0, 'w'},
       {"whitelist_rna", required_argument, 0, 'W'},
       {"dump", no_argument, 0, 'd'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"k", required_argument, 0, 'k'},
       {"species", required_argument, 0, 's'},
       {"countsfile", required_argument, 0, 'C'},
       {"speciesfile", required_argument, 0, 'S'},
       {"num_threads", required_argument, 0, 't'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string outdir = "";
    vector<string> atac_r1files;
    vector<string> atac_r2files;
    vector<string> atac_r3files;
    vector<string> rna_r1files;
    vector<string> rna_r2files;
    vector<string> custom_r1files;
    vector<string> custom_r2files;
    vector<string> custom_names;
    string whitelist_atac_filename;
    string whitelist_rna_filename;
    bool countsfile_given = false;
    string countsfilename;
    bool speciesfile_given = false;
    string speciesfilename;
    int num_threads = 1;
    double doublet_rate = 0.1;
    vector<string> kmerfiles;
    vector<string> speciesnames;
    bool dump = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "t:o:1:2:3:r:R:x:X:n:k:s:C:S:w:W:D:dh", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'o':
                outdir = optarg;
                break;           
            case 'k':
                kmerfiles.push_back(optarg);
                break;
            case 's':
                speciesnames.push_back(optarg);
                break;
            case '1':
                atac_r1files.push_back(optarg);
                break;
            case '2':
                atac_r2files.push_back(optarg);
                break;
            case '3':
                atac_r3files.push_back(optarg);
                break;
            case 'r':
                rna_r1files.push_back(optarg);
                break;
            case 'R':
                rna_r2files.push_back(optarg);
                break;
            case 'x':
                custom_r1files.push_back(optarg);
                break;
            case 'X':
                custom_r2files.push_back(optarg);
                break;
            case 'n':
                custom_names.push_back(optarg);
                break;
            case 'w':
                whitelist_atac_filename = optarg;
                break;
            case 'W':
                whitelist_rna_filename = optarg;
                break;
            case 'C':
                countsfile_given = true;
                countsfilename = optarg;
                break;
            case 'S':
                speciesfile_given = true;
                speciesfilename = optarg;
                break;
            case 'd':
                dump = true;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (outdir == "" && !dump){
        fprintf(stderr, "ERROR: output_directory / -o is required\n");
        exit(1);
    }
    else if (outdir == "."){
        outdir = "";
    } 
    else if (outdir != "" && outdir[outdir.length() - 1] != '/'){
        outdir += "/";
    }
    if (atac_r1files.size() == 0 && rna_r1files.size() == 0 && custom_r1files.size() == 0){
        // No reads given. Only okay if dumping data and species/counts files are given.
        if (!(dump && speciesfilename.length() > 0 && countsfilename.length() > 0)){
            fprintf(stderr, "ERROR: one or more of ATAC, RNA-seq, or custom \
(feature barcoding) types of reads are required.\n");
            fprintf(stderr, "This requirement can be avoided if you have chosen to dump data and \
are loading species-specific k-mer counts per barcode from files (with the species file and \
counts file options)\n");
            exit(1);
        }
    }
    if (kmerfiles.size() == 0 && (countsfilename == "" || speciesfilename == "")){
        fprintf(stderr, "ERROR: you must either load counts from a prior run using \
-C and -S, or specify k-mer count files using -k.\n");
        exit(1);
    }
    if (kmerfiles.size() > 0 && speciesnames.size() != kmerfiles.size()){
        fprintf(stderr, "ERROR: you must provide the same number and order of kmer files \
-k and species names -s.\n");
        exit(1);
    }
    if (atac_r1files.size() != atac_r2files.size() || atac_r1files.size() != atac_r3files.size()){
        fprintf(stderr, "ERROR: non-matching numbers of R1, R2, and/or R3 ATAC input files.\n");
        exit(1);
    }
    if (rna_r1files.size() != rna_r2files.size()){
        fprintf(stderr, "ERROR: non-matching numbers of R1 and R2 RNA-seq input files.\n");
        exit(1);
    }
    if (custom_r1files.size() != custom_r2files.size()){
        fprintf(stderr, "ERROR: non-matching numbers of R1 and R2 custom input files.\n");
        exit(1);
    }
    if (custom_r1files.size() != custom_names.size()){
        fprintf(stderr, "ERROR: you must provide a name/data type for each custom read \
file to demultiplex\n");
        exit(1);
    }
    if (atac_r1files.size() > 0 && (whitelist_atac_filename == "" || whitelist_rna_filename == "")){
        fprintf(stderr, "ERROR: if demultiplexing ATAC-seq reads from multiome data, \
both -w and -W are required.\n");
        exit(1);
    }
    if (rna_r1files.size() > 0 && whitelist_rna_filename == "" && (!countsfile_given || !speciesfile_given)){
        fprintf(stderr, "ERROR: if demultiplexing RNA-seq data and you are not loading pre-computed \
            k-mer counts, you must provide an RNA-seq bc whitelist (-W)\n");
        exit(1);
    }
    if (atac_r1files.size() > 0 && rna_r1files.size() > 0 && 
        (whitelist_atac_filename == "" || whitelist_rna_filename == "")){
        fprintf(stderr, "ERROR: demultiplexing combined ATAC and RNA-seq data (multiome), \
but whitelist files (-w and -W) not provided. Exiting.\n");
        exit(1);
    }
    if ((countsfile_given && !speciesfile_given) || (!countsfile_given && speciesfile_given)){
        fprintf(stderr, "ERROR: if loading data from a previous run, both --countsfile/-C \
and --speciesfile/-S are required.\n");
        exit(1);
    }
    if (num_threads > 100){
        fprintf(stderr, "ERROR: maximum number of threads (100) exceeded.\n");
        exit(1);
    }
    if (doublet_rate <= 0 || doublet_rate >= 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1, exclusive.\n");
        exit(1);
    }

    // Read whitelists
    // Not necessary if reading pre-generated species counts from a file AND
    // not demultiplexing reads.
    bc_whitelist wl;
    bool multiome = false;

    if (!(dump && countsfile_given)){
        if (whitelist_rna_filename != "" && whitelist_atac_filename != ""){
            wl.init(whitelist_rna_filename, whitelist_atac_filename);
            multiome = true;
        }
        else if (whitelist_rna_filename != ""){
            wl.init(whitelist_rna_filename);
            multiome = false;
        }
        else{
            wl.init(whitelist_atac_filename);
            multiome = false;
        }
    }
    
    // Map each species name to a numeric index
    map<short, string> idx2species;
    map<string, short> species2idx;

    // Peek at one of the k-mer count files to obtain k.
    int k;
    if (kmerfiles.size() > 0){
        ifstream peek(kmerfiles[0]);
        string line;
        peek >> line;
        k = line.length();
        fprintf(stderr, "Using k = %d\n", k);
    }
    
    robin_hood::unordered_map<unsigned long, map<short, int> > bc_species_counts;
    
    if (!countsfile_given){
        
        // Init species k-mer counter 
        species_kmer_counter counter(num_threads, k, kmerfiles.size(), &wl, &bc_species_counts);

        // i = species index
        for (int i = 0; i < kmerfiles.size(); ++i){
            
            // Parse k-mer file 
            fprintf(stderr, "Loading %s-specific k-mers\n", speciesnames[i].c_str());
            counter.init(i, kmerfiles[i]);
            fprintf(stderr, "done\n");
            
            idx2species.insert(make_pair(i, speciesnames[i]));
            species2idx.insert(make_pair(speciesnames[i], i));
            
            // Process RNA-seq data.
            counter.launch_gex_threads();    
            
            for (int i = 0; i < rna_r1files.size(); ++i){
                // Read through this set of files and launch jobs.
                counter.process_gex_files(rna_r1files[i], rna_r2files[i]); 
            }
            
            // Shut off running threads.
            counter.close_pool();
            
        }

        // Create file listing species names (in case we need to re-run without
        // loading the BAMs)
        string species_names_out = outdir;
        if (species_names_out[species_names_out.length()-1] != '/'){
            species_names_out += "/";
        }
        species_names_out += "species_names.txt";
        FILE* sn_out = fopen(species_names_out.c_str(), "w");
        for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end();
            ++i2s){
            fprintf(sn_out, "%d\t%s\n", i2s->first, i2s->second.c_str());
        }
        fclose(sn_out);
    
        // Create a counts file so we don't have to do the expensive process of counting 
        // k-mers next time, if we need to do something over.

        string countsfile_out = outdir;
        if (countsfile_out[countsfile_out.length()-1] != '/'){
            countsfile_out += "/";
        }
        countsfile_out += "species_counts.txt";
        FILE* countsfile = fopen(countsfile_out.c_str(), "w");         
        // Dump species counts
        print_bc_species_counts(bc_species_counts, idx2species, countsfile); 
        fclose(countsfile);
    
    }
    else{
        // Obtain species-specific k-mer counts from a data file (from a previous run), 
        // instead of counting k-mers in reads, which will take a long time.
        load_from_files(countsfilename, speciesfilename, idx2species, species2idx, 
            bc_species_counts);
    }
    
    // Now assign bcs to species.
    robin_hood::unordered_map<unsigned long, short> bc2species;
    robin_hood::unordered_map<unsigned long, pair<unsigned int, unsigned int> > bc2doublet;
    robin_hood::unordered_map<unsigned long, double> bc2llr;

    // Fit a mixture model to the data
    fit_model(bc_species_counts, bc2species, bc2doublet, bc2llr, idx2species,
        doublet_rate);

    FILE* bc_out;
    if (dump){
        bc_out = stdout;
    }
    else{
        string bc_out_name = outdir;
        if (bc_out_name[bc_out_name.length()-1] != '/'){
            bc_out_name += "/";
        }   
        bc_out_name += "bc_assignments.txt";
        bc_out = fopen(bc_out_name.c_str(), "w");
    }
    
    print_assignments(bc_out, bc2species, bc2doublet, bc2llr, idx2species);

    if (dump){
        // Our job is done here
        exit(0);
    }
    else{
        fclose(bc_out);
    }
    
    // See if we can/should create "library files" for multiome data, to save headaches later
    create_library_file(rna_r1files, atac_r1files, custom_r1files, 
        custom_names, idx2species, outdir);
    
    // Now go through reads and demultiplex by species.
    reads_demuxer demuxer(wl, bc2species, idx2species, outdir);

    for (int i = 0; i < atac_r1files.size(); ++i){
        fprintf(stderr, "Processing ATAC files %s, %s, and %s\n", 
            atac_r1files[i].c_str(), atac_r2files[i].c_str(), atac_r3files[i].c_str());
        demuxer.init_atac(atac_r1files[i], atac_r2files[i], atac_r3files[i]);
        demuxer.scan_atac();
    } 
    for (int i = 0; i < rna_r1files.size(); ++i){
        fprintf(stderr, "Processing RNA-seq files %s and %s\n", 
            rna_r1files[i].c_str(), rna_r2files[i].c_str());
        demuxer.init_rna(rna_r1files[i], rna_r2files[i]);
        demuxer.scan_rna();
    }
    for (int i = 0; i < custom_r1files.size(); ++i){
        fprintf(stderr, "Processing custom read files %s and %s\n", 
            custom_r1files[i].c_str(), custom_r2files[i].c_str());
        demuxer.init_custom(custom_names[i], custom_r1files[i], custom_r2files[i]);
        demuxer.scan_custom();
    } 

    return 0;
}