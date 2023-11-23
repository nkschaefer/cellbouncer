#include <getopt.h>
#include <argp.h>
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
#include <htswrapper/bc_hash.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include "reads_io.h"
#include "robin_hood.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "kmers_demux_species [OPTIONS]\n");
    fprintf(stderr, "Given reads from a multi-species pooled experiment and lists of species-specific k-mers, \
demultiplexes the reads by species and creates library files for running cellranger-arc\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n   ===== GENERAL OPTIONS =====\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    fprintf(stderr, "    --doublet_rate -D What is the prior expected doublet rate? (OPTIONAL; \
default = 0.1). Must be a decimal between 0 and 1, non-inclusive.\n");
    fprintf(stderr, "    --output_directory -o The directory in which to place output files. Read file names will be extracted from read group tags.\n");
    fprintf(stderr, "    --num_threads -t The number of threads to use for parallel processing (default 1)\n");
     fprintf(stderr, "    --dump -d Only dump per-barcode data (barcode, then count of reads per species \
(tab separated)) and barcode-to-species assignments instead of demultiplexing reads. \
These files are created in <output_directory>/species_counts.txt and \
<output_directory>/species_assignments.txt regardless. \
Standard behavior is to create this file and then demultiplex reads; this option causes \
the program to quit after generating the file.\n");
    fprintf(stderr, "\n   ===== READ FILE INPUT OPTIONS =====\n");
    fprintf(stderr, "    --atac_r1 -1 ATAC R1 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "    --atac_r2 -2 ATAC R2 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "    --atac_r3 -3 ATAC R3 reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "    --rna_r1 -r Forward RNA-seq reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "    --rna_r2 -R Reverse RNA-seq reads to demultiplex (can specify multiple times)\n");
    fprintf(stderr, "    --custom_r1 -x Forward other (i.e. sgRNA or antibody capture) reads to demultiplex \
(can specify multiple times). Assumes barcodes are at the beginning of R1.\n");
    fprintf(stderr, "    --custom_r2 -X Reverse other (i.e. sgRNA or antibody capture) reads to demultiplex \
(can specify multiple times). Assumes barcodes are at the beginning of R1.\n");
    fprintf(stderr, "    --names_custom -n Name of data type in custom reads file, in same number and order \
as read files. Presets: CRISPR = CRISPR sgRNA capture, Ab = antibody capture. Names will be appended \
to the beginning of demultiplexed FASTQ files and inserted into 10X library files. For example, if providing \
sgRNA capture files sgRNA_R1.fq.gz and sgRNA_R2.fq.gz along with antibody capture files anti_R1.fq.gz and anti_R2.fq.gz, \
you could specify -x sgRNA_R1.fq.gz -X sgRNA_R2.fq.gz -x anti_R1.fq.gz -X anti_R2.fq.gz -n CRISPR -n Antibody.\n");
    fprintf(stderr, "\n   ===== BARCODE WHITELIST OPTIONS =====\n");
    fprintf(stderr, "    --whitelist_atac -w If multiome data and demultiplexing ATAC-seq reads, provide both the ATAC-seq barcode whitelist (here) \
and the RNA-seq barcode whitelist (-W) (REQUIRED). If not multiome or RNA-seq only, this whitelist is not required.\n");
    fprintf(stderr, "    --whitelist_rna -W If multiome data and demultiplexing ATAC-seq reads, provide both the ATAC-seq barcode whitelist (-w) \
and the RNA-seq barcode whitelist (here) (REQUIRED). If not multiome or RNA-seq only, these are not required.\n");
    fprintf(stderr, "\n   ===== OTHER INPUT OPTIONS =====\n");
    fprintf(stderr, "    --k -k File containing species-specific k-mers. Specify once per species.\n");
    fprintf(stderr, "    --species -s Name of species corresponding to entries of -k. Specify \
multiple times, the same number and order as -k.\n");
    fprintf(stderr, "    --countsfile -C If a previous run already completed and you want to \
run again using the same counts, specify the counts file (<output_directory>/species_counts.txt) \
here. This will skip the expensive step of iterating through the BAM file to count reads. \
You can then either demultiplex reads (default) or dump data and quit (-d). If using this parameter, \
then --speciesfile/-S is also required.\n");
    fprintf(stderr, "    --speciesfile -S If loading data from a previous run with -C, then \
-S is also required. This parameter is a tab separated file that lists index, followed by \
human-readable species name. It will be generated by previous runs as \
<output_directory>/species_names.txt.\n");
    fprintf(stderr, "\n ===== NOTES =====\n");
    fprintf(stderr, "   This program works by counting k-mers in RNA-seq data exclusively. The other \
types of reads are provided to be demultiplexed only, by sharing of barcodes with the RNA-seq data. \
If you provide other types of data (i.e. ATAC, sgRNA capture), this program will attempt to create a \
10X Genomics-format library file to help run data. Any feature barcoding data will need an accompanying \
feature reference file, though, which must be created manually \
(see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).\n");
    fprintf(stderr, "   Once k-mers are counted once, it creates a counts file and species names \
file in the output directory. These can be loaded to demultiplex reads instead of repeating \
k-mer counting, which is the most expensive step.\n");
    fprintf(stderr, "   When counting k-mers, species-specific k-mer files (-k), species names \
(-s), RNA-seq reads (-r and -R), and an RNA-seq barcode whitelist (-W) are all required.\n");
    fprintf(stderr, "   When demultiplexing based on previously-computed k-mer counts, a counts file (-C), \
species names file (-S), and all reads to demultiplex (-r/-R, -1/-2/-3, -x/-X/-n) are required.\n");
    exit(code);
}

/**
 * Is the given character a digit between 0 and 9?
 */
bool isdigit(char ch){
    static char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    for (int i = 0; i < 10; ++i){
        if (ch == digits[i]){
            return true;
        }
    }
    return false;
}

/**
 * Truncate a file name to come up with a "base" file name that 
 * can be used to access/write other files in the same group.
 */
bool filename_base(string& filename, string& fnbase){
    stringstream splitter(filename);
    string chunk;
    vector<string> chunks;
    bool sfound = false;
    bool lfound = false;
    bool rfound = false;
    bool numfound = false;
    
    fnbase = "";    
    while (getline(splitter, chunk, '_')){
        chunks.push_back(chunk);
    }
    
    for (vector<string>::reverse_iterator chunk = chunks.rbegin(); chunk != chunks.rend(); ++chunk){
        if (!numfound){
            string suffix = chunk->substr(chunk->length()-8-1, 9);
            string prefix = chunk->substr(0, chunk->length()-8-1);
            bool pass = true;
            if (suffix == ".fastq.gz"){
                for (int i = 0; i < prefix.length(); ++i){
                    if (!isdigit(prefix[i])){
                        pass = false;
                        break;
                    }
                }
            }
            else{
                pass = false;
            }
            if (pass){
                numfound = true;
            }
            else{
                return false;
            }
        }
        else{
            if (!rfound){
                if ((*chunk)[0] == 'R' && isdigit((*chunk)[1])){
                    rfound = true;
                }
                else{
                    return false;
                }
            }
            else{
                if (!lfound){ 
                    if ((*chunk)[0] == 'L'){
                        bool pass = true;
                        for (int i = 1; i < chunk->length(); ++i){
                            if (!isdigit((*chunk)[i])){
                                pass = false;
                                break;
                            }
                        }
                        if (pass){
                            lfound = true;
                        }
                        else{
                            return false;
                        }
                    }
                }
                else{
                    if (!sfound){
                        if ((*chunk)[0] == 'S'){
                            bool pass = true;
                            for (int i = 1; i < chunk->length(); ++i){
                                if (!isdigit((*chunk)[i])){
                                    pass = false;
                                    break;
                                }
                            }
                            if (pass){
                                sfound = true;
                            }
                            else{
                                return false;
                            }
                        } 
                    }
                    else{
                        // Made it.
                        if (fnbase != ""){
                            fnbase = *chunk + "_" + fnbase;
                        }
                        else{
                            fnbase = *chunk;
                        }
                    }
                }
            }
        } 
    }
    return true;
}

/**
 * After demultiplexing, write a "stats" file containing some
 * summary data that can be reviewed.
 */
void create_stats_file(string& outdir,
    map<short, string>& idx2species,
    bcset& bc2species,
    robin_hood::unordered_map<unsigned long, map<short, double> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, double>& bc_species_tots,
    long int atac_reads_tot,
    map<short, long int>& atac_reads_pass,
    long int rna_reads_tot,
    map<short, long int>& rna_reads_pass,
    long int custom_reads_tot,
    map<short, long int>& custom_reads_pass){

    // Print data about how many barcodes were mapped to each species
    string statsfilename = outdir;
    if (outdir[outdir.length()-1] == '/'){
        statsfilename += "stats.txt";
    }
    else{
        statsfilename += "/stats.txt";
    }
    FILE* statsf = fopen(statsfilename.c_str(), "w");
    fprintf(statsf, "species\tnBarcodes\tnBarcodes.perc.assigned\tnBarcodes.perc.all\tatac_pairs\tatac_pairs.perc\trna_pairs\trna_pairs.perc\tcustom_pairs\tcustom_pairs.perc\n");
    map<short, int> species_cell_counts;
    for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); ++i2s){
        species_cell_counts.insert(make_pair(i2s->first, 0));
    }
    for (bcset::iterator b2s = bc2species.begin(); b2s != bc2species.end(); ++b2s){
        species_cell_counts[b2s->second]++;
        int bctot = 0;
        if (bc_species_tots.count(b2s->first) > 0){
            bctot = bc_species_tots[b2s->first];
        }
    }
    int species_all_tot = 0;
    for (map<short, int>::iterator scc = species_cell_counts.begin(); scc != species_cell_counts.end(); ++scc){
        species_all_tot += scc->second;
    }
    if (atac_reads_tot == 0){
        atac_reads_tot = 1;
    }
    if (rna_reads_tot == 0){
        rna_reads_tot = 1;
    }
    if (custom_reads_tot == 0){
        custom_reads_tot = 1;
    }
    for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); ++i2s){
        fprintf(statsf, "%s\t%d\t%.4f\t%.4f\t%ld\t%.4f\t%ld\t%.4f\t%ld\t%.4f\n", 
            i2s->second.c_str(), 
            species_cell_counts[i2s->first], 
            (float)species_cell_counts[i2s->first] / (float)species_all_tot,
            (float)species_cell_counts[i2s->first] / (float)bc_species_counts.size(),
            atac_reads_pass[i2s->first],
            (float)atac_reads_pass[i2s->first] / (float)atac_reads_tot,
            rna_reads_pass[i2s->first],
            (float)rna_reads_pass[i2s->first] / (float)rna_reads_tot,
            custom_reads_pass[i2s->first],
            (float)custom_reads_pass[i2s->first] / (float)custom_reads_tot);
    }
    fclose(statsf);
}

/**
 * Write out all species-specific k-mer counts per barcode in table format.
 */
void print_bc_species_counts(
    robin_hood::unordered_map<unsigned long, map<short, double> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, double>& bc_species_tots,
    map<short, string>& idx2species,
    FILE* countsfile){
    
    for (robin_hood::unordered_map<unsigned long, map<short, double> >::iterator bci = bc_species_counts.begin(); 
        bci != bc_species_counts.end(); ++bci){
        bc this_bc(bci->first);
        string bcstr = bc2str(this_bc, 16);
        fprintf(countsfile, "%s", bcstr.c_str());
        for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
            double count = 0;
            if (bci->second.count(spec->first) > 0){
                count = bci->second[spec->first];
            }   
            fprintf(countsfile, "\t%d", (int)round(count));
        }
        fprintf(countsfile, "\t%d\n", (int)round(bc_species_tots[bci->first]));
    }
}

/**
 * Instead of iterating through the BAM file again, can just load data from a previous run.
 * Requires counts file and species file both generated from previous runs.
 */
void load_from_files(string& countsfilename,
    string& speciesfilename,
    map<short, string>& idx2species,
    map<string, short>& species2idx,
    robin_hood::unordered_map<unsigned long, map<short, double> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, double>& bc_species_tots){
    
    fprintf(stderr, "Loading barcode-species counts\n"); 
    // Load species names
    ifstream speciesfile(speciesfilename);
    string species_idx_str;
    string species_name;
    while (speciesfile >> species_idx_str >> species_name){
        short species_idx = atoi(species_idx_str.c_str());
        idx2species.insert(make_pair(species_idx, species_name));
        species2idx.insert(make_pair(species_name, species_idx));
    }
    
    // Load count data
    string line;
    ifstream countsfile(countsfilename);
    while (getline(countsfile, line)){
        istringstream splitter(line);
        string field;
        short field_idx = 0;
        
        bc bc_bin;
        while(getline(splitter, field, '\t')){
            if (field_idx == 0){
                // Barcode.
                if (field[field.length()-2] == '-' && field[field.length()-1] == '1'){
                    field = field.substr(0, field.length()-2);
                }
                str2bc(field.c_str(), bc_bin, 16);   
                map<short, double> m;
                bc_species_counts.emplace(bc_bin.to_ulong(), m);
            }
            else if (field_idx <= species2idx.size()){
                // Field index == species ID + 1
                float count = atof(field.c_str());
                bc_species_counts[bc_bin.to_ulong()].insert(make_pair(field_idx-1, count));
            }
            else{
                // Optional last column = total # reads
                int tot = atoi(field.c_str());
                bc_species_tots.emplace(bc_bin.to_ulong(), tot);
            }
            field_idx++;
        }
    }
    fprintf(stderr, "done\n");
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
void fit_model(robin_hood::unordered_map<unsigned long, map<short, double> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, double>& bc_species_tots,
    bcset& bc2species, 
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
    vector<double> obs_weights; // weight observations by total read count
    vector<vector<double> > obs;
    vector<unsigned long> bcs;
    for (robin_hood::unordered_map<unsigned long, map<short, double> >::iterator x = 
        bc_species_counts.begin(); x != bc_species_counts.end(); ++x){
        obs_weights.push_back(bc_species_tots[x->first]);
        bcs.push_back(x->first);
        vector<double> row;
        for (map<short, double>::iterator sc = x->second.begin(); sc != x->second.end(); ++sc){
            row.push_back(sc->second);
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
    mod.fit(obs, obs_weights);
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

/**
 * After assigning barcodes to species, print assignments to an output file that
 * can be reviewed.
 */
void print_assignments(FILE* outf, 
    bcset& bc2species,
    robin_hood::unordered_map<unsigned long, pair<unsigned int, unsigned int> >& bc2doublet,
    robin_hood::unordered_map<unsigned long, double>& bc2llr,
    map<short, string>& idx2species){
    
    for (robin_hood::unordered_map<unsigned long, double>::iterator b2llr = bc2llr.begin();
        b2llr != bc2llr.end(); ++b2llr){
        
        char type = 'S';
        string name;
        if (bc2doublet.count(b2llr->first) > 0){
            type = 'D';
            int idx1 = bc2doublet[b2llr->first].first;
            int idx2 = bc2doublet[b2llr->first].second;
            if (idx2species[idx1] < idx2species[idx2]){
                name = idx2species[idx1] + "+" + idx2species[idx2];
            }
            else{
                name = idx2species[idx2] + "+" + idx2species[idx1];
            }
        }
        else if (bc2species.count(b2llr->first) > 0){
            type = 'S';
            name = idx2species[bc2species[b2llr->first]];
        }
        bc as_bitset(b2llr->first);
        string bc_str = bc2str(as_bitset, 16);
        fprintf(outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type,
            b2llr->second);
    }
}

/**
 * Create a 10x Genomics-style "Library file" to make it easier to run
 * cellranger-arc, if desired.
 */
void create_library_file(vector<string>& rna_r1files,
    vector<string>& atac_r1files,
    vector<string>& custom_r1files,
    vector<string>& custom_names,
    map<short, string>& idx2species,
    const string& outdir){

    char fullpath[200];
    if (rna_r1files.size() > 0 && (atac_r1files.size() > 0 || custom_r1files.size() > 0)){
        // Multiple data types. Create a "library file" for each species
        for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); ++i2s){
            realpath(outdir.c_str(), &fullpath[0]);
            string fnprefix = fullpath;
            if (fnprefix[fnprefix.length()-1] != '/'){
                fnprefix += "/";
            }
            string libfilename = fnprefix + i2s->second + ".library";
            fnprefix += i2s->second;
            FILE* libfile = fopen(libfilename.c_str(), "w");
            fprintf(libfile, "fastqs,sample,library_type\n");
            
            set<string> atac_base_pre;
            for (int i = 0; i < atac_r1files.size(); ++i){
                // Get base name
                string fn_trunc = filename_nopath(atac_r1files[i]); 
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                if (fn_base.length() < 5 || fn_base.substr(0, 5) != "ATAC_"){
                    fn_base = "ATAC_" + fn_base;
                }
                if (atac_base_pre.find(fn_base) == atac_base_pre.end()){
                    fprintf(libfile, "%s,%s,Chromatin Accessibility\n", fnprefix.c_str(),
                        fn_base.c_str());  
                }
                atac_base_pre.insert(fnprefix);
            }
            set<string> rna_base_pre;
            for (int i = 0; i < rna_r1files.size(); ++i){
                string fn_trunc = filename_nopath(rna_r1files[i]);
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                if (fn_base.length() < 4 || fn_base.substr(0, 4) != "GEX_"){
                    fn_base = "GEX_" + fn_base;
                }
                // Ensure only unique entries are written
                if (rna_base_pre.find(fn_base) == rna_base_pre.end()){
                    fprintf(libfile, "%s,%s,Gene Expression\n", fnprefix.c_str(), 
                        fn_base.c_str());
                }
                rna_base_pre.insert(fn_base);
            }
            set<string> custom_base_pre;
            set<string> custom_type_pre;
            for (int i = 0; i < custom_r1files.size(); ++i){
                string fn_trunc = filename_nopath(custom_r1files[i]);
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                // Determine what to append
                string prefix = custom_names[i] + "_";
                if (fn_base.length() < prefix.length() || 
                    fn_base.substr(0, prefix.length()) != prefix){
                    fn_base = prefix + fn_base;
                }
                string type_lib = custom_names[i];
                if (custom_names[i] == "CRISPR"){
                    type_lib = "CRISPR Guide Capture";
                }
                else if (custom_names[i] == "Ab"){
                    type_lib = "Antibody Capture";
                }
                if (custom_base_pre.find(fn_base) == custom_base_pre.end() ||
                    custom_type_pre.find(type_lib) == custom_type_pre.end()){
                    fprintf(libfile, "%s,%s,%s\n", fnprefix.c_str(), 
                        fn_base.c_str(), type_lib.c_str());
                }
                custom_base_pre.insert(fn_base);
                custom_type_pre.insert(type_lib);
            }
            fclose(libfile);
        }
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
    bcset atac_bc2idx;
    bcset rna_bc2idx;
    
    vector<unsigned long> whitelist_rna;
    if (!(dump && countsfile_given)){
        parse_whitelists(whitelist_atac_filename, whitelist_rna_filename,
           atac_bc2idx, whitelist_rna, rna_bc2idx);
    }
    
    // Map each species name to a numeric index (more space efficient)
    map<short, string> idx2species;
    map<string, short> species2idx;

    // Peek at one of the k-mer count files to obtain k.
    int k;
    kmer_tree_p kt;
    if (kmerfiles.size() > 0){
        ifstream peek(kmerfiles[0]);
        string line;
        peek >> line;
        k = line.length();
        
        kt = init_kmer_tree(k);
        fprintf(stderr, "Using k = %d\n", k);
    }
    
    // Load k-mer counts.
    if (kmerfiles.size() > 0){
        for (int i = 0; i < kmerfiles.size(); ++i){
            fprintf(stderr, "Parsing k-mer counts file %d\n", i);
            parse_kmer_counts_serial(kmerfiles[i], i, kt);
            fprintf(stderr, "done\n");
            idx2species.insert(make_pair(i, speciesnames[i]));
            species2idx.insert(make_pair(speciesnames[i], i));
        }
    }

    // Map a numeric representation of each barcode sequence to a another map that
    // ties species index to counts of species-specific k-mers from that species.
    robin_hood::unordered_map<unsigned long, map<short, double> > bc_species_counts;

    // Map a numeric representation of each barcode sequence to a total number
    // of reads matching that barcode.
    robin_hood::unordered_map<unsigned long, double> bc_tots;
    
    if (!countsfile_given){
        // Obtain species-specific k-mer counts from reads
        
        // Determine tmp output file name
        string tmp_out = outdir;
        if (tmp_out[tmp_out.length()-1] != '/'){
            tmp_out += "/";
        }
        tmp_out += "counts";

        // Init thread pool
        readsThreadPool pool(num_threads, idx2species.size(), tmp_out);
        pool.init();
    
        /*
         * NOTE: species-specific k-mer sets will be different for ATAC vs RNA-seq data,
         * because of splicing. Additionally, k-mer sets for ATAC will be massive and take
         * forever to load. For that reason (for now), we will only demultiplex based on
         * RNA-seq, but split ATAC-files based on what is learned about barcodes based on
         * the RNA-seq data.
         *
        if (atac_r1files.size() > 0){
            // Process ATAC data.
            pool.launch_rt_threads();
            for (int i = 0; i < atac_r1files.size(); ++i){

            }
        }
        */

        if (rna_r1files.size() > 0){
            
            // Process RNA-seq data.
            pool.launch_gex_threads(kt, 
                rna_bc2idx, whitelist_rna);
            
            for (int i = 0; i < rna_r1files.size(); ++i){
                // Read through this set of files and launch jobs.
                process_gex_files(rna_r1files[i], rna_r2files[i], pool); 
            }
            
            // Shut off running threads.
            pool.close_pool();
        }
        
        // Destroy k-mer to species mapping data to save memory
        kmsuftree_destruct(kt, 0);

        // Load counts we just computed from tmp files on disk
        pool.load_counts_from_tmp(bc_species_counts, bc_tots);

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
    }
    else{
        // Obtain species-specific k-mer counts from a data file (from a previous run), 
        // instead of counting k-mers in reads, which will take a long time.
        load_from_files(countsfilename, speciesfilename, idx2species, species2idx, 
            bc_species_counts, bc_tots);
    }
    
    if (!countsfile_given){
        
        // Create a counts file so we don't have to do the expensive process of counting 
        // k-mers next time, if we need to do something over.

        string countsfile_out = outdir;
        if (countsfile_out[countsfile_out.length()-1] != '/'){
            countsfile_out += "/";
        }
        countsfile_out += "species_counts.txt";
        FILE* countsfile = fopen(countsfile_out.c_str(), "w");         
        // Dump species counts
        print_bc_species_counts(bc_species_counts, bc_tots, idx2species, countsfile); 
        fclose(countsfile);
    }
    
    // Now assign bcs to species.
    bcset bc2species;
    robin_hood::unordered_map<unsigned long, pair<unsigned int, unsigned int> > bc2doublet;
    robin_hood::unordered_map<unsigned long, double> bc2llr;

    // Fit a mixture model to the data
    fit_model(bc_species_counts, bc_tots, bc2species, bc2doublet, bc2llr, idx2species, doublet_rate);

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
    
    long int atac_reads_tot = 0;
    map<short, long int> atac_reads_pass;
    long int rna_reads_tot = 0;
    map<short, long int> rna_reads_pass;
    long int custom_reads_tot = 0;
    map<short, long int> custom_reads_pass;
    for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); 
        ++i2s){
        atac_reads_pass.insert(make_pair(i2s->first, 0));
        rna_reads_pass.insert(make_pair(i2s->first, 0));
        custom_reads_pass.insert(make_pair(i2s->first, 0));
    }

    for (int i = 0; i < atac_r1files.size(); ++i){
        fprintf(stderr, "Processing ATAC files %s, %s, and %s\n", 
            atac_r1files[i].c_str(), atac_r2files[i].c_str(), atac_r3files[i].c_str());
        demux_atac_reads(atac_r1files[i], atac_r2files[i], atac_r3files[i], bc2species, 
            idx2species, outdir, atac_bc2idx, whitelist_rna, atac_reads_tot, 
            atac_reads_pass);    
    } 
    for (int i = 0; i < rna_r1files.size(); ++i){
        fprintf(stderr, "Processing RNA-seq files %s and %s\n", 
            rna_r1files[i].c_str(), rna_r2files[i].c_str());
        string rna_prefix = "GEX";
        demux_rna_reads(rna_r1files[i], rna_r2files[i], rna_prefix, 
            bc2species, idx2species, outdir, rna_reads_tot, rna_reads_pass);
    }
    for (int i = 0; i < custom_r1files.size(); ++i){
        fprintf(stderr, "Processing custom read files %s and %s\n", 
            custom_r1files[i].c_str(), custom_r2files[i].c_str());
        demux_rna_reads(custom_r1files[i], custom_r2files[i], custom_names[i],
            bc2species, idx2species, outdir, custom_reads_tot, custom_reads_pass);
    } 

    // Finally, create file with demultiplexing-related stats.
    create_stats_file(outdir, idx2species, bc2species, bc_species_counts,
        bc_tots, atac_reads_tot, atac_reads_pass, rna_reads_tot,
        rna_reads_pass, custom_reads_tot, custom_reads_pass);

    return 0;
}
