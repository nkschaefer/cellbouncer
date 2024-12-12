#define STRINGIZE_(x) #x
#define STRINGIZE(x) STRINGIZE_(x)
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
#include <map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/seq_fuzzy_match.h>
#include <htswrapper/bc_scanner.h>
#include <mixtureDist/functions.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <optimML/functions.h>
#include <optimML/multivar_sys.h>
#include <optimML/brent.h>
#include "common.h"
#include <htswrapper/gzreader.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_tags [OPTIONS]\n");
    fprintf(stderr, "Given reads containing cell hashing/MULTIseq OR CRISPR sgRNA capture\n");
    fprintf(stderr, "  data, counts the occurrences of each label or sgRNA per cell barcode\n");
    fprintf(stderr, "  and assigns each cell barcode an identity. Cells can receive anywhere\n");
    fprintf(stderr, "  from zero identities to a combination of every possible identity.\n");
    fprintf(stderr, "  In the case of cell hashing/MULTIseq, multiple identities correspond to\n");
    fprintf(stderr, "  droplets containing multiple cells. In the case of sgRNA data, they\n");
    fprintf(stderr, "  can reflect either multiplets OR cells that received multiple sgRNAs.\n");
    fprintf(stderr, "A single run will produce a .counts file, an .assignments file, an .ids\n");
    fprintf(stderr, "  file, and a .table file.\n");
    fprintf(stderr, "\n   ===== OUTPUTS =====\n");
    fprintf(stderr, "  The .counts file is a table of counts of each identity for each cell\n");
    fprintf(stderr, "    barcode. If you run again with the same output prefix, this file will\n");
    fprintf(stderr, "    be loaded rather than re-counting labels in reads.\n");
    fprintf(stderr, "  The .assignments file is a table of cell barcode identities. Without the\n");
    fprintf(stderr, "    --sgRNA option set, columns are cell barcode, identity, S/D/M (singlet,\n");
    fprintf(stderr, "    doublet, or multiplet), and log likelihood ratio of best to second best\n");
    fprintf(stderr, "    assignment. With the --sgRNA option set, columns are cell barcode, identity,\n");
    fprintf(stderr, "    number of assignments (1 or more), and log likelihood ratio of best to second\n");
    fprintf(stderr, "    best assignment.\n");
    fprintf(stderr, "  The .table file is created only when the --sgRNA option is set. It is a\n");
    fprintf(stderr, "    differently-formatted table of cell barcode identities, more useful for high\n");
    fprintf(stderr, "    MOI experiments, where multiple IDs are expected. Columns are cell barcode,\n");
    fprintf(stderr, "    followed by every possible sgRNA assignment. Each cell receives a 0 for\n");
    fprintf(stderr, "    each guide it was not assigned, and a 1 for each guide it was assigned.\n");
    fprintf(stderr, "  The .bg file is a table of cell barcode and inferred percent ambient (background)\n");
    fprintf(stderr, "    tags. This will include cells that were removed from final analysis.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n   ===== REQUIRED =====\n");
    fprintf(stderr, "   --output_prefix -o The prefix for output file names. Will create\n");
    fprintf(stderr, "       a file ending in .counts and another ending in .assignments.\n");
    fprintf(stderr, "\n   ===== TO COUNT TAGS IN FASTQ FILES =====\n");
    fprintf(stderr, "     ----- REQUIRED -----\n");
    fprintf(stderr, "   --read1 -1 Forward read file for cell hashing/MULTIseq/sgRNA capture data. Can\n");
    fprintf(stderr, "       specify multiple times.\n");
    fprintf(stderr, "   --read2 -2 Reverse read file for cell hashing/MULTIseq/sgRNA capture data. Can\n");
    fprintf(stderr, "       specify multiple times.\n");
    fprintf(stderr, "   --whitelist -w Cell barcode whitelist file (i.e. included with 10X Genomics\n");
    fprintf(stderr, "       Cellranger software in cellranger-x.y.z/lib/python/cellranger/barcodes/\n");
    fprintf(stderr, "   --exact_only -e Only accept cell barcode matches to allowed barcode list when\n");
    fprintf(stderr, "       they are exact. This will speed up tag counting in reads and significantly\n");
    fprintf(stderr, "       decrease memory usage, at the cost of picking up slightly fewer tag counts\n");
    fprintf(stderr, "       per cell barcode.\n");
    fprintf(stderr, "   --seqs -s File mapping sequences (i.e. MULTI-seq barcodes or sgRNA capture\n");
    fprintf(stderr, "       sequences) to IDs. These can either be final assignments (i.e. sgRNA names)\n");
    fprintf(stderr, "       or intermediate IDs, which can then be mapped to final IDs via the --mapping\n");
    fprintf(stderr, "       file. For example, you can use a file that maps MULTIseq barcodes to well IDs\n");
    fprintf(stderr, "       for all experiments as the --seqs file, and then provide a second --mapping\n");
    fprintf(stderr, "       file for each experiment that assigns a name/identity to each well. If you do\n");
    fprintf(stderr, "       not provide a --mapping file, names in the --seqs file will be in final output.\n");
    fprintf(stderr, "       File should be tab-separated, with one sequence and one name per line.\n");
    fprintf(stderr, "   --names -N (OPTIONAL) If you are using a technology where there is a consistent mapping\n");
    fprintf(stderr, "       of sequences to intermediate names (i.e. MULTIseq barcodes to well names), but where\n");
    fprintf(stderr, "       the meaning of each intermediate name changes in every experiment (i.e. the MULTIseq\n");
    fprintf(stderr, "       well represents a different individual or treatment in every experiment), then provide\n");
    fprintf(stderr, "       the sequence to intermediate name mapping via --seqs/-s and the intermediate to final\n");
    fprintf(stderr, "       name mapping here. This should be a 2-column tab separated file with lines containin\n");
    fprintf(stderr, "       the intermediate name followed by the final name.\n");
    fprintf(stderr, "       If you do not provide an argument here, the names in the --seqs/-s file will be used\n");
    fprintf(stderr, "       in output files.\n");
    fprintf(stderr, "       One advantage of using this argument is that the program will check to see if any\n");
    fprintf(stderr, "       sequences in the --seqs/-s file not given names in this file are more highly represented\n");
    fprintf(stderr, "       in your data than expected. This can help identify when you have mixed up sample wells.\n");
    fprintf(stderr, "       ----- OPTIONAL -----\n");
    fprintf(stderr, "   --mismatches -m The number of allowable mismatches to a MULTIseq barcode/HTO/sgRNA\n");
    fprintf(stderr, "       capture sequence to accept it. Default = 2; set to -1 to take best match overall,\n");
    fprintf(stderr, "       if there are no ties.\n");
    fprintf(stderr, "   --filt -f Attempt to find \"negative\" cell barcodes that should not receive an ID.\n");
    fprintf(stderr, "       Default: attempt to assign an ID to all cell barcodes, unless in --sgRNA mode,\n");
    fprintf(stderr, "       in which case filtering will happen even without this option set.\n");
    fprintf(stderr, "   --umi_len -u This program assumes that forward reads begin with a cell barcode followed\n");
    fprintf(stderr, "       by a UMI. Cell barcode length is set during compilation with the -D BC_LENX2\n");
    fprintf(stderr, "        argument (default = 32, meaning 16 bp barcodes). This sets the length of the UMI\n");
    fprintf(stderr, "       (default = 12 bp; cannot exceed 16 bp). Some versions of 10X kits use 10 bp UMIs.\n");
    fprintf(stderr, "   --exact -e Require cell barcodes to exactly match those in the whitelist. Default:\n");
    fprintf(stderr, "       allow one mismatch. This will improve speed, especially for large whitelists.\n");
    fprintf(stderr, "\n     ===== ALTERNATIVE INPUT OPTION: MEX format =====\n");
    fprintf(stderr, "   If you have already counted tag barcodes using another program (i.e. 10X Genomics\n");
    fprintf(stderr, "       CellRanger or kallisto/bustools kite), you can load market exchange (MEX) format\n");
    fprintf(stderr, "       output via these options instead of using this program to count tags in reads.\n");
    fprintf(stderr, "   MEX format consists of three files, to be provided below:\n");
    fprintf(stderr, "   --mtx -M The matrix/.mtx file (optionally gzipped)\n");
    fprintf(stderr, "   --features -F The features/\"genes\" file (optionally gzipped)\n");
    fprintf(stderr, "   --cell_barcodes -B (same argument as below): the list of cell barcodes (optionally\n");
    fprintf(stderr, "       gzipped)\n");
    fprintf(stderr, "   --feature_type -t If the MEX-format files contain multiple data types (i.e. if you\n");
    fprintf(stderr, "       processed RNA-seq and sgRNA capture data in the same CellRanger run, then specify\n");
    fprintf(stderr, "       the feature type for tag data to extract. For 10X cell hashing, for example,\n");
    fprintf(stderr, "       specify \"Multiplexing\\ Capture\", for antibody capture specify \"Antibody\\ Capture\",\n");
    fprintf(stderr, "       and for sgRNA capture specify \"CRISPR\\ Guide\\ Capture\".\n");
    fprintf(stderr, "\n   ===== STRONGLY RECOMMENDED =====\n");
    fprintf(stderr, "   --cell_barcodes -B The filtered list of valid cell barcodes, i.e. from cellranger or\n");
    fprintf(stderr, "       STARsolo. Only cells in this list will be processed.\n");
    fprintf(stderr, "       NOTE: if using market exchange-format input (see above), this argument is required.\n");
    fprintf(stderr, "\n   ===== OPTIONAL =====\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    fprintf(stderr, "   --sgRNA -g Specifies that data is from sgRNA capture. This affects how sequence\n");
    fprintf(stderr, "       matching in reads is done and slightly changes the output files. Also makes\n");
    fprintf(stderr, "       the default separator for multiple IDs a comma instead of +.\n");
    fprintf(stderr, "   --prob -p If you are processing sgRNA capture data (--sgRNA option), this parameter\n");
    fprintf(stderr, "       sets the minimum probability a guide is present in a cell to be assigned. Because\n");
    fprintf(stderr, "       of how the algorithm works, this is actually the conditional probability the guide\n");
    fprintf(stderr, "       is present, conditional on the presence of any likelier guides. Default = 0.5\n");
    fprintf(stderr, "   --comma -c By default, cells assigned multiple identities will have these\n");
    fprintf(stderr, "       identities separated by + in the .assignments and .ids files. If\n");
    fprintf(stderr, "       identities contain + within them, though, this option switches to a\n");
    fprintf(stderr, "       comma separator.\n");
    print_libname_help(); 
    exit(code);
}

/**
 * Returns barcode length.
 */
int load_seq2well(string& filename, 
    map<string, string>& seq2well, 
    vector<string>& seqlist){
    
    ifstream inf(filename);
    string seq_str;
    string well;
    int seq_len = -1;
    while (inf >> seq_str >> well){
        if (seq_len == -1){
            seq_len = seq_str.length();
        }
        else if (seq_len != seq_str.length()){
            fprintf(stderr, "ERROR: mismatching sequence lengths in file %s:\n",
                filename.c_str());
            fprintf(stderr, "%d vs %ld\n", seq_len, seq_str.length());
            exit(1);
        }
        seq2well.insert(make_pair(seq_str, well));
        seqlist.push_back(seq_str);
    }
    return seq_len;
}

void load_well_mapping(string& filename, map<string, string>& well2id){
    ifstream inf(filename);
    string well;
    string id;
    while (inf >> well >> id){
        if (well2id.count(well) > 0){
            fprintf(stderr, "ERROR: well %s is assigned to at least 2 identities:\n", well.c_str());
            fprintf(stderr, "%s\n", well2id[well].c_str());
            fprintf(stderr, "%s\n", id.c_str());
            exit(1);
        }
        well2id.insert(make_pair(well, id));
    }
}

/**
 * Write counts to disk for future inspection / loading on a subsequent run
 *
 * Populates bc_tag_counts (data structure used by other functions to represent
 * final counts) from bc_ms_umis (data structure used by counting function to 
 * deduplicate UMIs)
 */
void dump_counts(string& filename, 
    robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >& bc_ms_umis,
    vector<string>& seq_names,
    robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    string& libname,
    bool cellranger,
    bool seurat,
    bool underscore){

    FILE* f = fopen(filename.c_str(), "w");
    
    // Print header line
    fprintf(f, "barcode");
    for (vector<string>::iterator n = seq_names.begin(); n != seq_names.end(); ++n){
        if (n->length() > 0){
            fprintf(f, "\t%s", n->c_str());
        }
    }
    fprintf(f, "\n");

    for (robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >::iterator x =
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){

        map<int, long int> m;
        bc_tag_counts.emplace(x->first, m);

        string bc_str = bc2str(x->first);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore);
        
        fprintf(f, "%s", bc_str.c_str());
        // j == index in bc_tag_counts
        int j = 0;
        for (int i = 0; i < x->second.size(); ++i){
            if (seq_names[i].length() > 0){
                // This MULTIseq barcode was intended to be present in the data set
                if (x->second[i] == NULL){
                    fprintf(f, "\t0");
                }
                else{
                    bc_tag_counts[x->first].insert(make_pair(j, x->second[i]->count()));
                    fprintf(f, "\t%d", x->second[i]->count());
                }
                ++j;
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

/**
 * Alternate method for writing counts to disk, in case they were loaded from
 * MEX-format files.
 */
void dump_counts2(string& countsfn, vector<string>& seqnames,
    robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    string& libname,
    bool cellranger,
    bool seurat,
    bool underscore){
    
    FILE* f = fopen(countsfn.c_str(), "w");
    fprintf(f, "barcode");
    for (int i = 0; i < seqnames.size(); ++i){
        fprintf(f, "\t%s", seqnames[i].c_str());
    }
    fprintf(f, "\n");

    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator btc = 
        bc_tag_counts.begin(); btc != bc_tag_counts.end(); ++btc){
        string bcstr = bc2str(btc->first);
        mod_bc_libname(bcstr, libname, cellranger, seurat, underscore);
        fprintf(f, "%s", bcstr.c_str());
        for (int i = 0; i < seqnames.size(); ++i){
            long int count = 0;
            if (btc->second.count(i) > 0){
                count = btc->second[i];
            }
            fprintf(f, "\t%ld", count);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

/**
 * If counting was already done on a previous run and the user just wants to identify cells,
 * retrieve counts from the output file using this function.
 */
int load_counts(string& filename,
    robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    vector<string>& labels){
    
    string line;
    ifstream inf(filename);
    bool first = true;

    vector<string> header;

    while (getline(inf, line)){
        istringstream splitter(line);
        string field;
        int idx = 0;
        unsigned long cur_bc;
        while(getline(splitter, field, '\t')){
            if (first){
                // Ensure names are unique
                for (int x = 0; x < header.size(); ++x){
                    if (header[x] == field){
                        fprintf(stderr, "ERROR: label %s appears multiple times in header.\n", field.c_str());
                        exit(1);
                    }
                }
                header.push_back(field);
            }
            else{
                if (idx == 0){
                    // Trim off any extra text from the end of the barcode sequence
                    cur_bc = bc_ul(field);
                    map<int, long int> m;
                    bc_tag_counts.emplace(cur_bc, m);      
                }
                else{
                    long int val = atol(field.c_str());
                    if (val > 0){
                        bc_tag_counts[cur_bc].insert(make_pair(idx-1, val));
                    }
                }
            }
            ++idx;
        }

        first = false;
    }

    for (int i = 0; i < header.size()-1; ++i){
        labels.push_back(header[i+1]);
    }
    return header.size()-1;
}

/**
 * For use with Brent 1-dim solver: infer proportion background
 * in a cell - adjusted to handle arbitrary number of foreground identities
 */
double ll_mix_multi(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int ndim = data_i.at("ndim");

    vector<double> x;
    vector<double> p;
    
    char buf[50];
    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        if (data_d.count(bufstr) == 0){
            continue;
        }
        double x_i = data_d.at(bufstr);
        
        x.push_back(x_i);

        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgfrac = data_d.at(bufstr);
        
        // Check whether foreground - and get mixing proportion if so
        sprintf(&buf[0], "m_%d", i);
        bufstr = buf;
        double m = data_d.at(bufstr);

        if (m > 0){
            p.push_back(bgfrac * param + (1.0 - param) * m);
        }
        else{
            p.push_back(bgfrac * param);
        }
    }
    return dmultinom(x, p)/log2(exp(1));
}

/**
 * For use with Brent's 1-D solver: derivative function to 
 * find the MLE percent background in a cell
 */
double dll_mix_multi(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int ndim = data_i.at("ndim");

    char buf[50];

    double dL_df = 0.0;

    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        if (data_d.count(bufstr) == 0){
            continue;
        }
        double x_i = data_d.at(bufstr);
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgfrac = data_d.at(bufstr);
        
        sprintf(&buf[0], "m_%d", i);
        bufstr = buf;
        double m = data_d.at(bufstr);
        if (m > 0){
            dL_df += (x_i*(bgfrac - m))/((1-param)*m + param*bgfrac);
        }
        else{
            // Background
            dL_df += x_i / param;
        }
    }
    
    return dL_df;
}

/**
 * Simple linear regression routine given a vector of independent variable (x)
 * and a vector of dependent variable (y).
 *
 * Returns a pair where .first is slope and .second is intercept.
 */
pair<double, double> linreg(const vector<double>& x, const vector<double>& y){
    double ysum = 0.0;
    double ysum_sq = 0.0;
    double xsum = 0.0;
    double xsum_sq = 0.0;
    double xysum = 0.0;

    for (int i = 0; i < x.size(); ++i){
        xysum += x[i] * y[i];
        xsum += x[i];
        ysum += y[i];
        xsum_sq += x[i] * x[i];
        ysum_sq += y[i] * y[i];
    }

    double num = (double)x.size();
   
    double slope = (num * xysum - xsum*ysum)/(num * xsum_sq - xsum*xsum);
    double intercept = (ysum * xsum_sq - xsum * xysum) / 
        (num * xsum_sq - xsum * xsum);
    return make_pair(slope, intercept);

}

/**
 * Weighted version of above.
 */
pair<double, double> linreg_weights(const vector<double>& x, const vector<double>& y, 
    vector<double>& weights){
     
    double weighted_sum_x = 0.0;
    double weighted_sum_y = 0.0;
    double weighted_sum_xy = 0.0;
    double weighted_sum_x2 = 0.0;
    double total_weight = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double w = weights[i];
        weighted_sum_x += w * x[i];
        weighted_sum_y += w * y[i];
        weighted_sum_xy += w * x[i] * y[i];
        weighted_sum_x2 += w * x[i] * x[i];
        total_weight += w;
    }

    double mean_x = weighted_sum_x / total_weight;
    double mean_y = weighted_sum_y / total_weight;

    double numerator = weighted_sum_xy - total_weight * mean_x * mean_y;
    double denominator = weighted_sum_x2 - total_weight * mean_x * mean_x;

    double slope = numerator / denominator;
    double intercept = mean_y - slope * mean_x;
    return make_pair(slope, intercept);

}

/**
 * Given a cell and its likeliest identity, determine the likeliest percent
 * of its reads composed of background/ambient counts via maximum likelihood.
 */
double get_cell_bg(vector<double>& obsrow,
    set<int>& assns,
    vector<double>& model_means,
    vector<double>& bgfracs,
    int ndim){
    
    // Infer background prop
    optimML::brent_solver mix(ll_mix_multi, dll_mix_multi);
    mix.set_maxiter(-1);
    double thistot = 0.0;
    for (set<int>::iterator a = assns.begin(); a != assns.end(); ++a){
        thistot += model_means[*a];
    }
    
    // If we can't solve, then fall back on the sum of all off-target counts
    double this_offtarget = 0.0;
    double this_tot = 0.0;

    vector<vector<double> > n;
    vector<vector<double> > bg;
    vector<vector<double> > m;

    mix.add_data_fixed("ndim", ndim);
    for (int k = 0; k < ndim; ++k){
        if (model_means[k] > 0){
            n.push_back(vector<double>{ obsrow[k] });
            bg.push_back(vector<double>{ bgfracs[k] });
            if (assns.find(k) == assns.end()){
                this_offtarget += obsrow[k];
                m.push_back(vector<double>{ -1.0 });
            }
            else{
                m.push_back(vector<double>{ model_means[k] / thistot });
            }
        }
        else{
            n.push_back(vector<double>{ } );
            bg.push_back(vector<double>{ } );
            m.push_back(vector<double>{ } );
        }
        this_tot += obsrow[k];
    }

    char buf[50];
    for (int k = 0; k < ndim; ++k){
        if (model_means[k] > 0){
            sprintf(&buf[0], "n_%d", k);
            string bufstr = buf;
            mix.add_data(bufstr, n[k]);
            sprintf(&buf[0], "bg_%d", k);
            bufstr = buf;
            mix.add_data(bufstr, bg[k]);
            sprintf(&buf[0], "m_%d", k);
            bufstr = buf;
            mix.add_data(bufstr, m[k]);
        }
    }
    
    mix.constrain_01();
    double p_inferred = mix.solve(0,1);
    //mix.print(0,1,0.001);
    //fprintf(stderr, "P = %f\n", p_inferred);
    if (!mix.root_found){
        return this_offtarget/obsrow[obsrow.size()-1];
    } 
    else{
        return p_inferred;
    }
}

/**
 * Function used by Brent's method solver to find the point of intersection
 * between a negative binomial distribution and an exponential distribution.
 */
double diff_negbin_exp(double x, const map<string, double>& data_d, const map<string, int>& data_i){
    double mu = data_d.at("mu");
    double phi = data_d.at("phi");
    double lambda = data_d.at("lambda");

    // Brent solver is set to maximize, we want to minimize. So make negative.
    double ll1 = dnbinom(x, mu, phi);
    double ll2 = dexp(x, lambda);
    
    return ll1-ll2;
}

/**
 * First step: fit a 2-way mixture model to counts of each label.
 * Background is modeled as a negative binomial over low counts,
 * and foreground is modeled as an exponential over high counts.
 *
 * Stores mu/phi of negative binomial (low) dist in lomeans and lophis
 * Stores mean of exponential in himeans
 * Stores mean count of high distribution (1/exponential parameter) 
 *   in model_means, or -1 if the label is filtered out
 */
void fit_dists_2way(vector<vector<double> >& obscol,
    vector<string>& labels,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    string& output_prefix){
    
    vector<double> loweights;
    vector<double> hiweights;
    
    string outfn = output_prefix + ".dists";
    FILE* outf = fopen(outfn.c_str(), "w");
    
    vector<double> dist_seps;
    double dist_sep_min = 0.0;

    fprintf(stderr, "===== Initial model means: =====\n");
    for (int i = 0; i < labels.size(); ++i){
        if (obscol[i].size() < 3){
            lomeans.push_back(-1);
            lophis.push_back(-1);
            himeans.push_back(-1);
        }
        else{
            sort(obscol[i].begin(), obscol[i].end());
            double lowval = obscol[i][0];
            double lowval2;
            if (obscol[i].size() % 2 == 0){
                double lowval2a = obscol[i][obscol[i].size()/2];
                double lowval2b = obscol[i][obscol[i].size()/2-1];
                lowval2 = (lowval2a + lowval2b)/2.0;
            }
            else{
                lowval2 = obscol[i][(obscol[i].size()-1)/2];
            }
            double highval = obscol[i][obscol[i].size()-1];
            vector<mixtureDist> dists;
            pair<int, double> lomm = nbinom_moments(lowval, pow((lowval2-lowval)/2.0, 2));
            
            mixtureDist low("negative_binomial", vector<double>{ lowval, 1.0 });
            mixtureDist high("exponential", 1.0/highval);
            dists.push_back(low);
            dists.push_back(high);
            mixtureModel mod(dists);
            mod.fit(obscol[i]);
            
            double lomean = mod.dists[0].params[0][0];
            double lophi = mod.dists[0].params[0][1];
            double himean = 1.0/mod.dists[1].params[0][0];
            double w = mod.weights[0];
            
            mixtureDist low2("negative_binomial", vector<double>{ lowval2, 1.0 });
            mixtureDist high2("exponential", 1.0/highval);
            vector<mixtureDist> dists2;
            dists2.push_back(low2);
            dists2.push_back(high2);
            mixtureModel mod2(dists2);
            mod2.fit(obscol[i]);
            
            // If model2 is a better fit, take it.
            if (mod2.loglik > mod.loglik && 
                1.0/mod2.dists[1].params[0][0] > mod2.dists[0].params[0][0] && 
                mod2.weights[0] > 0 && mod2.weights[1] > 0 &&
                1.0/mod2.dists[1].params[0][0] > 1.0/mod.dists[1].params[0][0]){
                
                lomean = mod2.dists[0].params[0][0];
                lophi = mod2.dists[0].params[0][1];
                himean = 1.0/mod2.dists[1].params[0][0];
                w = mod2.weights[0];
            }       
            
            double dist_sep = pnbinom(himean, lomean, lophi);
            if (dist_sep > 1-1e-4){
                dist_sep = 1-1e-4;
            }
            else if (dist_sep < 1e-4){
                dist_sep = 1e-4;
            }
            dist_seps.push_back(dist_sep);
            if (dist_sep_min == 0.0 || dist_sep < dist_sep_min){
                dist_sep_min = dist_sep;
            }
            fprintf(stderr, "%s:\t%f\t%f\n", labels[i].c_str(), lomean-1, himean-1);
            fprintf(outf, "%s\t%f\t%f\t%f\t%f\n", labels[i].c_str(),
                w, lomean, lophi, 1.0/himean);

            lomeans.push_back(lomean);
            lophis.push_back(lophi);
            himeans.push_back(himean);
            hiweights.push_back(1.0-w);
            loweights.push_back(w);
        }
    }

    fclose(outf);
    
    // Now test whether any labels appear to have dropped out. This happens in two steps:
    
    // 1. Compute the "separation" between the two dists fit for each label.
    //  This is the CDF of the lower (negative binomial) dist evaluated at the mean
    //  of the upper dist.
    //  Then fit a two-way Beta mixture model to these values. If it seems bimodal, then
    //  toss out all tags in the lower component (less separated).
    
    // 2. Of tags that survive step 1, compute the mean of all lower components and the
    //  mean of all upper components. Treat each as negative binomial. Then toss out any
    //  tags for which the upper component is more likely to have originated from the mean
    //  lower component dist than the mean upper component dist.

    fprintf(stderr, "Checking for label dropout...\n");
    mixtureDist seplo("beta", 10.0*dist_sep_min, 10.0*(1.0-dist_sep_min));
    mixtureDist sephi("beta", 10.0*0.999, 10.0*(0.001));
    vector<mixtureDist> dists;
    dists.push_back(seplo);
    dists.push_back(sephi);
    mixtureModel sepmod(dists);
    sepmod.fit(dist_seps);
    
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < labels.size(); ++i){
        if (sepmod.assignments[i] == 0){
            count1++;
        }
        else if (sepmod.assignments[i] == 1){
            count2++;
        }
    }
    if (count1 > 0 && count2 > 0){
        double a1 = sepmod.dists[0].params[0][0];
        double b1 = sepmod.dists[0].params[0][1];
        double a2 = sepmod.dists[1].params[0][0];
        double b2 = sepmod.dists[1].params[0][1];
        double mu1 = a1/(a1+b1);
        double mu2 = a2/(a2+b2);
        double p1 = pbeta(mu1, a2, b2);
        double p2 = pbeta(mu2, a1, b1);
        if (p1 <= 0.01 && p2 >= 0.99){
            // Remove lower component.
            for (int i = 0; i < labels.size(); ++i){
                if (sepmod.assignments[i] == 0){
                    fprintf(stderr, "Remove label %s\n", labels[i].c_str());
                    lomeans[i] = -1;
                    lophis[i] = -1;
                    himeans[i] = -1;
                }
            }
        }
    }

    // Finally now compute mean background dist
    double bg_mu_sum = 0.0;
    double bg_var_sum = 0.0;
    int bg_mu_count = 0;
    
    double fg_mu_sum = 0.0;
    double fg_var_sum = 0.0;
    int fg_mu_count = 0;

    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] > 0){
            bg_mu_sum += lomeans[i];
            bg_var_sum += (lomeans[i] + (lomeans[i]*lomeans[i])/lophis[i]);
            bg_mu_count++;
        }
        if (himeans[i] > 0){
            fg_mu_sum += himeans[i];
            fg_var_sum += (himeans[i]*himeans[i]);
            fg_mu_count++;
        }
    }
    bg_mu_sum /= (double)bg_mu_count;
    bg_var_sum /= pow((double)bg_mu_count, 2);
    //double bg_phi_sum = (bg_mu_sum*bg_mu_sum)/(bg_var_sum - bg_mu_sum);
    
    fg_mu_sum /= (double)fg_mu_count;
    fg_var_sum /= pow((double)fg_mu_count, 2);
    pair<double, double> muphi_fg = nbinom_moments(fg_mu_sum, fg_var_sum);
    pair<double, double> muphi_bg = nbinom_moments(bg_mu_sum, bg_var_sum);
    
    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] > 0){
            double ll1 = dnbinom(himeans[i], muphi_bg.first, muphi_bg.second);
            double ll2 = dnbinom(himeans[i], muphi_fg.first, muphi_fg.second);
            if (ll1 > ll2){
                fprintf(stderr, "Remove label %s\n", labels[i].c_str());
                lomeans[i] = -1;
                lophis[i] = -1;
                himeans[i] = -1;
            }
        }
    }
}

// Define a struct that will let us use a set of int as a map key.
struct modelkey{
    set<int> members;
};

bool operator<(const modelkey& k1, const modelkey& k2){
    if (k1.members.size() < k2.members.size()){
        return true;
    }
    else if (k2.members.size() < k1.members.size()){
        return false;
    }
    else{
        // Same length
        set<int>::iterator x1 = k1.members.begin();
        set<int>::iterator x2 = k2.members.begin();
        while (x1 != k1.members.end() && x2 != k2.members.end()){
            if (*x1 < *x2){
                return true;
            }
            else if (*x2 < *x1){
                return false;
            }
            ++x1;
            ++x2;
        }
        // Equal
        return false;
    }
}

/**
 * optimML style log likelihood function for fitting a Dirichlet GLM
 * to total foreground counts (predicts percent of each tag present in a cell)
 */
double ll_dir(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int nmod = data_i.at("nmod");
    double tot = data_d.at("tot");
    double ltot = log(tot);

    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    
    int signp;
    
    char buf[30];
    string bufstr;
    for (int i = 0; i < nmod; ++i){
        double a = exp(params[i*2]*ltot + params[i*2 + 1]);
        
        sprintf(&buf[0], "p_%d", i);
        bufstr = buf;
        double p = data_d.at(bufstr);
        
        term1 += a;
        term2 += lgammaf_r(a, &signp);
        term3 += (a - 1.0)*log(p);
        
    }
    return lgammaf_r(term1, &signp) - term2 + term3;
}

/**
 * optimML style log likelihood function for fitting a multinomial GLM
 * to total foreground counts (predicts percent of each tag present in a cell)
 */
double ll_multn(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int nmod = data_i.at("nmod");
    double tot = data_d.at("tot");
    double ltot = log(tot);
    
    char buf[30];
    string bufstr;
    
    double psum = 0.0;
    vector<double> ps(nmod);
    vector<double> xs(nmod);

    for (int i = 0; i < nmod; ++i){
        double p = expit(params[i*2]*ltot + params[i*2 + 1]);
       
        sprintf(&buf[0], "x_%d", i);
        bufstr = buf;
        double x = data_d.at(bufstr);
        
        xs[i] = x;
        ps[i] = p;
        psum += p; 
    }
    for (int i = 0; i < nmod; ++i){
        ps[i] /= psum;
    }

    return dmultinom(xs, ps)/log2(exp(1));
}

/**
 * optimML style log likelihood derivative function for fitting Dirichlet
 * GLM (predicts fraction of each tag present in a cell based on total
 * foreground counts)
 */
void dll_dir(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
     
    int nmod = data_i.at("nmod");
    double tot = data_d.at("tot");
    double ltot = log(tot);

    vector<double> dll_da;
    vector<double> a_all;

    char buf[30];
    string bufstr;
    double term1 = 0.0;
    for (int i = 0; i < nmod; ++i){
        double a = exp(params[i*2]*ltot + params[i*2 + 1]);
        a_all.push_back(a);

        sprintf(&buf[0], "p_%d", i);
        bufstr = buf;
        double p = data_d.at(bufstr);
        term1 += a;
        
        dll_da.push_back(-digamma(a) + log(p));
    }
    for (int i = 0; i < nmod; ++i){
        dll_da[i] += digamma(term1);
        results[i*2] += dll_da[i] * ltot * a_all[i];
        results[i*2 + 1] += dll_da[i] * a_all[i];
    }
}

/**
 * optimML style log likelihood derivative function for fitting multinomial
 * GLM (predicts fraction of each tag present in a cell based on total
 * foreground counts)
 */
void dll_multn(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    int nmod = data_i.at("nmod");
    double tot = data_d.at("tot");
    double ltot = log(tot);

    char buf[30];
    string bufstr;
    
    vector<double> dll_dp(nmod);
    vector<double> xs(nmod);
    vector<double> e_negx1(nmod);
    vector<double> dp_dt(nmod);

    double psum = 0.0;

    for (int i = 0; i < nmod; ++i){
        double val = params[i*2]*ltot + params[i*2 + 1];
        e_negx1[i] = exp(-val);
        double p = expit(val);
        psum += p;

        sprintf(&buf[0], "x_%d", i);
        bufstr = buf;
        double x = data_d.at(bufstr);
        
        xs[i] = x;
        dll_dp[i] = p;
    }
    
    for (int i = 0; i < nmod; ++i){
        // Calculate dLL/dp
        double p = dll_dp[i] / psum;
        dll_dp[i] = xs[i] / p;
        
        // Calculate dp/dt
        double e_negx1_p1_2 = pow(e_negx1[i] + 1, 2);
        double e_negx1_p1_3 = e_negx1_p1_2*(e_negx1[i] + 1);
        double dp_dt = e_negx1[i] / (e_negx1_p1_2 * psum) - 
            e_negx1[i] / (e_negx1_p1_3 * psum * psum);
        
        // Calculate dt/dparams[i*2]
        double dt_dslope = ltot;

        // Calculate dt/dparams[i*2+1]
        double dt_dintercept = 1.0;

        results[i*2] += (dll_dp[i] * dp_dt * dt_dslope);
        results[i*2+1] += (dll_dp[i] * dp_dt * dt_dintercept);
    }
}

/**
 * Given some parameters and a desired model (tags parameter - a set of 
 * all tags present in the identity), figures out how to predict the
 * distribution of tags present in the cell under the model.
 *
 * This is done by fitting a Dirichlet GLM. Input (independent variable)
 * is total number of foreground tag counts in the cell. Output (dependent variable)
 * is fraction of all tags in the cell coming from each tag.
 *
 * Data are simulated draws from the distributions fit in fit_dists_2way().
 *
 * The goal here is to account for different sequencing depth across cells. If testing
 * a model, we don't want to just look at the predicted total counts, but the 
 * predicted fraction of counts conditional on the number of observed reads.
 */
void get_model(map<int, int>& idx2idx,
    set<int>& tags,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    vector<vector<double> >& bg_samps, 
    vector<double>& bg_tots,
    map<int, vector<double> >& model_samps,
    double& count_mu,
    double& count_phi,
    vector<double>& dircoeffs,
    bool sgrna,
    int nthreads){
    
    int nsamp = bg_samps.size();
    int nmod = idx2idx.size();
    
    bool success = false;
    while (!success){
        vector<double> logtots;
        vector<double> logpercs;
        vector<double> logitpercs;

        vector<double> tots;
        vector<double> bgs;
        vector<double> fgs;
        vector<double> pbgs;

        vector<vector<double> > samps_frac;
        vector<vector<double> > samps_column;
        vector<vector<double> > samps_frac_column_logit;
        vector<double> inset_log;
        vector<double> inset;

        for (int i = 0; i < nmod; ++i){
            vector<double> v;
            samps_column.push_back(v);
            samps_frac_column_logit.push_back(v);
        }
        
        vector<vector<double> > incl_columns;
        vector<vector<double> > incl_columns_logit;

        for (set<int>::iterator t = tags.begin(); t != tags.end(); ++t){
            vector<double> v;
            incl_columns.push_back(v);
            incl_columns_logit.push_back(v);
        }
        
        for (int i = 0; i < nsamp; ++i){
            vector<double> row = bg_samps[i];
            double tot = bg_tots[i];
            double totsamp = 0.0;
             
            map<int, double> percmember;
            
            double il = 0.0;
            int mi = 0;
            for (set<int>::iterator member = tags.begin(); member != tags.end(); ++member){
                double this_samp = model_samps[*member][i];
                totsamp += this_samp;
                int ri = idx2idx[*member];
                tot -= row[ri];
                row[ri] = this_samp;
                //row[ri] += this_samp;
                tot += this_samp;
                il += this_samp;
                incl_columns[mi].push_back(this_samp);
                ++mi;
            }
            inset_log.push_back(log(il));
            inset.push_back(il);
            
            fgs.push_back(totsamp);

            logtots.push_back(log(tot));
            logpercs.push_back(log(1.0-totsamp/tot));
            logitpercs.push_back(log(1.0-totsamp/tot)-log(totsamp/tot));
            bgs.push_back(tot-totsamp);
            pbgs.push_back(1.0-totsamp/tot);
            tots.push_back(tot);

            for (int mi2 = 0; mi2 < mi; ++mi2){
                if (incl_columns[mi2][i] == tot){
                    incl_columns_logit[mi2].push_back(logit((incl_columns[mi2][i]-1)/tot));
                }
                else{
                    incl_columns_logit[mi2].push_back(logit((incl_columns[mi2][i]+1)/(tot+1)));
                }
            }

            for (int j = 0; j < nmod; ++j){
                if (row[j] == tot){
                    samps_frac_column_logit[j].push_back(logit((row[j]-1)/tot));
                }
                else{
                    samps_frac_column_logit[j].push_back(logit((row[j] + 1)/(tot+1)));
                }
                //row[j] /= tot;
                samps_column[j].push_back(row[j]);
            }
            
            samps_frac.push_back(row);
        } 
        
        pair<double, double> mv = welford(fgs);
        pair<double, double> mp = nbinom_moments(mv.first, mv.second);
        count_mu = mp.first;
        count_phi = mp.second;

        try{ 
            if (sgrna){
                // This block should be used if fitting a model that only considers
                // fraction of counts from each foreground tag, and total background counts.
                vector<double> params;
                for (int j = 0; j < tags.size(); ++j){
                    pair<double, double> si = linreg(inset_log, incl_columns_logit[j]);
                    params.push_back(si.first);
                    params.push_back(si.second);
                }
                pair<double, double> sibg = linreg(inset_log, logitpercs);
                params.push_back(sibg.first);
                params.push_back(sibg.second);
                optimML::multivar_ml_solver solv(params, ll_multn, dll_multn);
                solv.add_data_fixed("nmod", (int)(tags.size() + 1));
                solv.add_data("tot", inset);
                solv.set_threads(nthreads);
                char bufx[100];
                string bufxstr;
                for (int j = 0; j < tags.size(); ++j){
                    sprintf(&bufx[0], "x_%d", j);
                    bufxstr = bufx;
                    solv.add_data(bufxstr, incl_columns[j]);
                }
                sprintf(&bufx[0], "x_%d", (int)tags.size());
                bufxstr = bufx;
                solv.add_data(bufxstr, pbgs);
                solv.solve();
                dircoeffs = solv.results;
                success = true;

            }
            else{

                // This block should be used if fitting a model that considers counts of
                // each tag.

                vector<double> params_d;
                for (int j = 0; j < nmod; ++j){
                    pair<double, double> si = linreg(inset_log, samps_frac_column_logit[j]);
                    params_d.push_back(si.first);
                    params_d.push_back(si.second);
                }
                optimML::multivar_ml_solver solv(params_d, ll_multn, dll_multn);
                solv.add_data_fixed("nmod", nmod);
                char buf1[100];
                string buf1str;
                for (int j = 0; j < nmod; ++j){
                    sprintf(&buf1[0], "x_%d", j);
                    buf1str = buf1;
                    solv.add_data(buf1str, samps_column[j]);
                }
                solv.add_data("tot", inset);
                solv.set_threads(nthreads);
                solv.solve();
                dircoeffs = solv.results;
                success = true;
            }
            /* 
            // Block to print predicted vs actual fractions    
            for (int i = 0; i < nsamp; ++i){
                vector<double> preds;
                double predtot = 0.0;
                for (int j = 0; j < nmod; ++j){
                    double pred = expit(solv.results[j*2]*log(inset[i]) + solv.results[j*2+1]);
                    preds.push_back(pred);
                    predtot += pred;
                }
                for (int j = 0; j < nmod; ++j){
                    fprintf(stdout, "%f\t%d\t%f\t%f\n", log(inset[i]), j, samps_column[j][i]/tots[i], preds[j]/predtot);
                }
            }
            exit(0);
            */

        }
        catch (int err){
            // Weird math error with GLM. Re-sample counts and try again.
            fprintf(stderr, "Resampling...\n");
            for (int i = 0; i < nsamp; ++i){
                // Background
                for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
                    bg_samps[i][ix->second] = rnbinom(lomeans[ix->first], lophis[ix->first]);
                }
                // Foreground
                for (set<int>::iterator t = tags.begin(); t != tags.end(); ++t){
                    model_samps[*t][i] = rexp(1.0/himeans[*t]);
                }
            }        
        }
    } 
}

double assign_cells_aux(vector<vector<double> >& obs,
    vector<double>& obstots,
    vector<vector<double> >& samps,
    vector<vector<double> >& samps_frac,
    vector<double>& rowtots,
    map<int, vector<double> >& model_samps,
    int nthreads,
    vector<unsigned long>& bcs,
    vector<string>& labels,
    vector<int>& models_kept,
    map<int, int>& idx2idx,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    modelkey& key_bg,
    double bg_count_mu,
    double bg_count_phi,
    map<modelkey, vector<double> >& model_coeffs,
    map<modelkey, pair<double, double> >& model_count_params,
    map<modelkey, vector<double> >& dirdists,
    map<modelkey, vector<vector<double> > >& dirsamps,
    map<modelkey, vector<double> >& dirmeans,
    map<modelkey, vector<double> >& dirtots,
    vector<double>& bgprops,
    vector<set<int> >& results,
    vector<double>& llrs,
    vector<double>& bgs,
    bool twopass_pass1,
    bool twopass_pass2,
    bool filt,
    bool sgrna){
    
    double prob = 0.5; 
    double llsum = 0.0;

    for (int i = 0; i < obs.size(); ++i){
        vector<pair<double, int> > lls;
        vector<string> names;
        vector<int> nums;

        vector<double> obsrowclean;

        // Sort in decreasing order of enrichment of p over expectation in background
        vector<pair<double, int> > dim_enrich_sort;
        
        double tot = obstots[i];
        if (tot == 0){
            set<int> s; 
            results.push_back(s);
            llrs.push_back(0.0);
            bgs.push_back(0.0);

            continue;
        } 
        bool all_bg = false;
        double bg_llr = 0.0;
        if (filt || sgrna){
            all_bg = true;
        }
        for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
            double p = obs[i][ix->first]/tot - bgprops[ix->second];
            dim_enrich_sort.push_back(make_pair(-p, ix->first));
            obsrowclean.push_back(obs[i][ix->first]);
            if (filt){
                double ll1 = dnbinom(obs[i][ix->first], lomeans[ix->first], lophis[ix->first]);
                double ll2 = dexp(obs[i][ix->first], 1.0/himeans[ix->first]);
                if (ll2 > ll1 & obs[i][ix->first] > lomeans[ix->first]){
                    all_bg = false;
                }
                else{
                    bg_llr += (ll1-ll2);
                }
            }    
        }

        if ((filt || sgrna) && all_bg){
            set<int> s;
            results.push_back(s);
            bgs.push_back(1.0);
            llrs.push_back(bg_llr);
            llsum += dmultinom(obsrowclean, bgprops);
            continue;
        }

        sort(dim_enrich_sort.begin(), dim_enrich_sort.end());
        double tot_counts = 0;
        modelkey included;
        double mean = 0.0;
        double var = 0.0;
        vector<double> llsize;
        
        double llsecond = 0.0;
        double llmax = 0.0;
        double llprev = 0.0;
        
        set<int> result;
        double nonbg = 0.0;
        
        long int maxn = dim_enrich_sort.size();
        if (!sgrna){
            maxn = 2;
        }
        for (int k = 0; k < maxn; ++k){
            int j = dim_enrich_sort[k].second;
            nonbg += obs[i][j];

            included.members.insert(j);
            
            if (twopass_pass2 && dirdists.count(included) == 0){
                continue;
            }    

            vector<double> coeffs_model;
            double count_mu;
            double count_phi;

            if (model_coeffs.count(included) == 0){
                // Generate parameters
                
                get_model(idx2idx,
                    included.members,
                    lomeans,
                    lophis,
                    himeans,
                    samps,
                    rowtots,
                    model_samps,
                    count_mu,
                    count_phi,
                    coeffs_model,
                    sgrna,
                    nthreads);

                // Store so we don't need to re-generate next time
                model_coeffs.insert(make_pair(included, coeffs_model));
                model_count_params.insert(make_pair(included, make_pair(count_mu, count_phi)));
            }   
            else{
                coeffs_model = model_coeffs[included];
                pair<double, double> mp = model_count_params[included];
                count_mu = mp.first;
                count_phi = mp.second;
            }
            

            // Multinomial params            
            vector<double> p(obsrowclean.size());
            double partot = 0.0;

            if (sgrna){
                int mod_idx = 0;
                vector<double> p_fg;
                for (int x = 0; x < included.members.size(); ++x){
                    double this_p = expit(coeffs_model[x*2]*log(nonbg) + coeffs_model[x*2+1]);
                    partot += this_p;
                    p_fg.push_back(this_p);
                    mod_idx++;
                }
                double p_bg = expit(coeffs_model[mod_idx*2]*log(nonbg) + coeffs_model[mod_idx*2+1]);
                partot += p_bg;
                for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
                    p[ix->second] = (p_bg/partot)*bgprops[ix->second];
                }
                mod_idx = 0;
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    p[idx2idx[*m]] = p_fg[mod_idx]/partot;
                    mod_idx++;
                }
            }
            else{
                for (int x = 0; x < obsrowclean.size(); ++x){
                    double this_p = expit(coeffs_model[x*2]*log(nonbg) + coeffs_model[x*2+1]);
                    partot += this_p;
                    p[x] = this_p;
                }
                for (int x = 0; x < obsrowclean.size(); ++x){
                    p[x] /= partot;
                }
            }

            double ll = dmultinom(obsrowclean, p);
            ll += dnbinom(nonbg, count_mu, count_phi);
            
            if (twopass_pass2){
                vector<double> obsrowprop = obsrowclean;
                for (int x = 0; x < obsrowprop.size(); ++x){
                    obsrowprop[x] /= tot;
                }
                ll += ddirichlet(obsrowprop, dirdists[included]);
            }

            // Things are sorted so that once log likelihood starts to dip, we can stop looking
            // (it will only continue to decrease).
            if (llmax != 0.0 && ll < llmax){
                if (llsecond == 0.0 || llsecond < ll){
                    llsecond = ll;
                }
                if (result.size() == 0){
                    fprintf(stderr, "%s\t%f\t%f\t%ld\n", bc2str(bcs[i]).c_str(), ll, llmax, result.size());
                    fprintf(stderr, "%s %f\n", labels[j].c_str(), -dim_enrich_sort[k].first);
                    for (int x = 0; x < obsrowclean.size(); ++x){
                        fprintf(stderr, "%s) %f\t%f\n", labels[models_kept[x]].c_str(), obsrowclean[x]/tot, bgprops[x]);  
                    }
                    fprintf(stderr, "\n");
                }
                break;
            }
            else if (ll == llmax){
                // This should not happen often.
                llsecond = llmax;
                // Don't add in result - keep more parsimonious result (fewer assignments/guides)
            }
            else{
                bool keep = true;
                // Check if we're requiring a minimum probability to keep.
                if (prob > 0 && prob > 0.5 && llprev != 0.0){
                    // Check probability to see if we should bail.
                    double denom = pow(2, llprev - llprev) + pow(2, ll-llprev);
                    double thisprob = pow(2, ll-llprev) / denom;
                    if (thisprob < prob){
                        // Don't accept.
                        keep = false;
                        break;
                    }
                }
                if (keep){
                    llsecond = llmax;
                    llmax = ll;
                    llprev = ll;
                    result.insert(j);
                }
            }
        }
        
        double pbg;
        bool all_bg_mod = false;
        if (sgrna || filt){
            double ll0 = dmultinom(obsrowclean, bgprops);
            ll0 += dnbinom(tot, bg_count_mu, bg_count_phi);
            
            if (twopass_pass2){
                vector<double> obsrowprop = obsrowclean;
                for (int x = 0; x < obsrowprop.size(); ++x){
                    obsrowprop[x] /= tot;
                }
                ll0 += ddirichlet(obsrowprop, dirdists[key_bg]);
            }

            if (ll0 > llmax){
                all_bg_mod = true;
                result.clear();
                llsecond = llmax;
                llmax = ll0;
                pbg = 1.0;
            }
            else if (ll0 > llsecond){
                llsecond = ll0;
            }
        }
        if (!all_bg_mod){
            pbg = get_cell_bg(obs[i], result, himeans, bgprops, himeans.size());
        }
        
        double tot_pseudo = tot + (double)obsrowclean.size();
        if (twopass_pass1){
            modelkey key;
            if (all_bg_mod){
                key = key_bg;
            }
            else{
                key = included;
            }
            if (dirsamps.count(key) == 0){
                vector<vector<double> > v1;
                dirsamps.insert(make_pair(key, v1));
                vector<double> v2;
                dirmeans.insert(make_pair(key, v2));
                dirtots.insert(make_pair(key, v2));
                for (map<int, int>::iterator idx = idx2idx.begin(); idx != idx2idx.end(); ++idx){
                    dirmeans[key].push_back(1.0);
                    dirtots[key].push_back(1.0);
                    vector<double> v3;
                    dirsamps[key].push_back(v3);
                }
            }
            for (map<int, int>::iterator idx = idx2idx.begin(); idx != idx2idx.end(); ++idx){
                dirmeans[key][idx->second] += obsrowclean[idx->second];
                dirtots[key][idx->second]++;
                dirsamps[key][idx->second].push_back((obsrowclean[idx->second]+1)/tot_pseudo);
            }
        }

        llsum += llmax;
        
        results.push_back(result);
        llrs.push_back(llmax-llsecond);
        bgs.push_back(pbg);
    }
    return llsum;
}

/**
 * Main function that assigns all cells to identities.
 */
double assign_cells(vector<vector<double> >& obs,
    vector<unsigned long>& bcs,
    int nsamp,
    vector<string>& labels,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    bool sgrna,
    vector<set<int> >& results,
    vector<double>& llrs,
    vector<double>& bgs,
    bool filt,
    int nthreads,
    vector<double>& bgprops,
    bool twopass){
    
    bgprops.clear();

    map<modelkey, vector<double> > model_coeffs;
    map<modelkey, pair<double, double> > model_count_params;
    
    map<modelkey, vector<double> > dirdists;    
    map<modelkey, vector<vector<double> > > dirsamps;
    map<modelkey, vector<double> > dirmeans;
    map<modelkey, vector<double> > dirtots;
    
    // Get model for all background
    modelkey key_bg;
    vector<int> models_kept;
    map<int, int> idx2idx;
    for (int j = 0; j < lomeans.size(); ++j){
        if (lomeans[j] > 0){
            idx2idx.insert(make_pair(j, models_kept.size()));
            models_kept.push_back(j);
        }
    }
    vector<double> meanconc(models_kept.size());
    for (int j = 0; j < models_kept.size(); ++j){
        meanconc[j] = 0.0;
    }
    vector<vector<double> > samps;
    vector<vector<double> > samps_frac;
    vector<double> rowtots;

    for (int i = 0; i < nsamp; ++i){
        vector<double> samp;
        double tot = 0.0;
        for (int j = 0; j < models_kept.size(); ++j){
            int model = models_kept[j];
            samp.push_back(rnbinom(lomeans[model], lophis[model]));
            if (samp[j] < 0){
                fprintf(stderr, "??? %f\n", samp[j]);
                fprintf(stderr, "%d, %d\n", j, model);
                fprintf(stderr, "%f %f\n", lomeans[model], lophis[model]);
                fprintf(stderr, "%s\n", labels[model].c_str());
                exit(1);
            }
            meanconc[j] += (1/(double)nsamp)*samp[j];
            tot += samp[j];
        }
        rowtots.push_back(tot);

        vector<double> sampfrac = samp;
        for (int j = 0; j < models_kept.size(); ++j){
            sampfrac[j] /= tot;
        }
        samps.push_back(samp);
        samps_frac.push_back(sampfrac);
    }
    
    pair<double, double> bgmv = welford(rowtots);
    pair<double, double> bgmp = nbinom_moments(bgmv.first, bgmv.second);
    double bg_count_mu = bgmp.first;
    double bg_count_phi = bgmp.second;

    double meanconctot = 0.0;
    for (int i = 0; i < meanconc.size(); ++i){
        meanconctot += meanconc[i];
    } 
    for (int i = 0; i < meanconc.size(); ++i){
        meanconc[i] /= meanconctot;
    }

    //vector<double> bgprops;

    fprintf(stderr, "===== Ambient tag profile: =====\n");
    //for (int j = 0; j < dirprops_bg.size(); ++j){
    for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
        bgprops.push_back(meanconc[ix->second]);
        //bgprops.push_back(dirprops_bg[ix->second] / bgtot);
        fprintf(stderr, "%s:\t%f\n", labels[ix->first].c_str(), meanconc[ix->second]); 
        //fprintf(stderr, "%s:\t%f\n", labels[ix->first].c_str(), bgprops[ix->second]);
    }
   
    double prob = 0.5;
    
    /*
    vector<double> hi_nbinom_mu;
    vector<double> hi_nbinom_phi;
    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] <= 0 || himeans[i] <= i){
            hi_nbinom_mu.push_back(-1);
            hi_nbinom_phi.push_back(-1);
        }
        else{
            hi_nbinom_mu.push_back(himeans[i] - lomeans[i]);
            double hivar = himeans[i]*himeans[i];
            double lovar = (lomeans[i]*lomeans[i])/lophis[i];
            double hivar2 = hivar - lovar;
            pair<double, double> mp = nbinom_moments(himeans[i] - lomeans[i], hivar2);
            hi_nbinom_phi.push_back(mp.second);
        }
    }
    */

    // Pre-sample counts for each allowed model.
    map<int, vector<double> > model_samps;
    for (int i = 0; i < models_kept.size(); ++i){
        int model_idx = models_kept[i];
        vector<double> s(nsamp);
        for (int j = 0; j < nsamp; ++j){
            //s[j] = rnbinom(hi_nbinom_mu[model_idx], hi_nbinom_phi[model_idx]); 
            s[j] = rexp(1.0/himeans[model_idx]);
            //double mu = himeans[model_idx];
            //s[j] = rnbinom(mu, (mu*mu)/(mu*mu - mu)); 
        }
        model_samps.insert(make_pair(model_idx, s));
    }

    int ndim = idx2idx.size();
    
    vector<double> obstots;
    for (int i = 0; i < obs.size(); ++i){
        double tot = 0.0;
        for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
            tot += obs[i][ix->first];
        }
        obstots.push_back(tot);
    }
    
    fprintf(stderr, "Assigning cells...\n");
    double llsum = assign_cells_aux(obs,
        obstots,
        samps,
        samps_frac,
        rowtots,
        model_samps,
        nthreads,
        bcs,
        labels,
        models_kept,
        idx2idx,
        lomeans,
        lophis,
        himeans,
        key_bg,
        bg_count_mu,
        bg_count_phi,
        model_coeffs,
        model_count_params,
        dirdists,
        dirsamps,
        dirmeans,
        dirtots,
        bgprops,
        results,
        llrs,
        bgs,
        twopass,
        false,
        filt,
        sgrna);

    if (twopass){
        // Get Dirichlet dists
        for (map<modelkey, vector<double> >::iterator dist = dirmeans.begin(); dist != dirmeans.end(); ++dist){
            for (int i = 0; i < dist->second.size(); ++i){
                dist->second[i] /= dirtots[dist->first][i];
            }
            vector<double> dirparams;
            fit_dirichlet(dist->second, dirsamps[dist->first], dirparams, nthreads);
            dirdists.insert(make_pair(dist->first, dirparams)); 
        }

        fprintf(stderr, "Second round of assignments using Empirical Bayes...\n");
        double llsum = assign_cells_aux(obs,
            obstots,
            samps,
            samps_frac,
            rowtots,
            model_samps,
            nthreads,
            bcs,
            labels,
            models_kept,
            idx2idx,
            lomeans,
            lophis,
            himeans,
            key_bg,
            bg_count_mu,
            bg_count_phi,
            model_coeffs,
            model_count_params,
            dirdists,
            dirsamps,
            dirmeans,
            dirtots,
            bgprops,
            results,
            llrs,
            bgs,
            false,
            true,
            filt,
            sgrna);

    }
    return llsum;
}


/**
 * Returns difference in log likelihood between two Beta distributions.
 * Used for finding intersection point via Brent's root finding method.
 */
double ll_beta_diff(double x, 
    const map<string, double>& data_d, 
    const map<string, int>& data_i){

    double a1 = data_d.at("a1");
    double b1 = data_d.at("b1");
    double a2 = data_d.at("a2");
    double b2 = data_d.at("b2");

    double ll1 = dbeta(x, a1, b1);
    double ll2 = dbeta(x, a2, b2);
    return ll1-ll2;
}

pair<double, double> modes(vector<double>& nums){
    map<double, int> histlo;
    map<double, int> histhi;
    int maxcountlo = -1;
    double maxbinlo = -1;
    int maxcounthi = -1;
    double maxbinhi = -1;
    for (int i = 0; i < nums.size(); ++i){
        double rounded = round(nums[i]*100.0)/100.0;
        if (nums[i] < 0.5){
            if (histlo.count(rounded) == 0){
                histlo.insert(make_pair(rounded, 0));
            }
            histlo[rounded]++;
            if (maxcountlo == -1 || histlo[rounded] > maxcountlo){
                maxcountlo = histlo[rounded];
                maxbinlo = rounded;
            }
        }
        else{
            if (histhi.count(rounded) == 0){
                histhi.insert(make_pair(rounded, 0));
            }
            histhi[rounded]++;
            if (maxcounthi == -1 || histhi[rounded] > maxcounthi){
                maxcounthi = histhi[rounded];
                maxbinhi = rounded;
            }
        }
    }
    return make_pair(maxbinlo, maxbinhi);
}

/**
 * Model percent background counts as 2 dists:
 * low percent background (true positives) and high percent background
 * (false positives).
 *
 * Returns cutoff value between the two distributions.
 */
double sep_bg_dists_percent(vector<double>& bgperc_all_flat){
    
    if (bimod_test(bgperc_all_flat)){
        mixtureModel bgmod;
        mixtureModel bgpercmod; 
        double intptx = -1;
        
        // Get initial values: mean % bg when below 50% and mean % bg when
        // above 50% 
        vector<double> bgperc_flat_low;
        vector<double> bgperc_flat_high;
        
        vector<vector<double> > bgperc_all;

        for (int i = 0; i < bgperc_all_flat.size(); ++i){
            if (bgperc_all_flat[i] > 0 && bgperc_all_flat[i] < 1.0){
                if (bgperc_all_flat[i] < 0.5){
                    bgperc_flat_low.push_back(bgperc_all_flat[i]);
                }
                else{
                    bgperc_flat_high.push_back(bgperc_all_flat[i]);
                }
                bgperc_all.push_back(vector<double>{ bgperc_all_flat[i] });
            }
        }

        pair<double, double> muvar_bgl = welford(bgperc_flat_low);
        pair<double, double> muvar_bgh = welford(bgperc_flat_high);
        pair<double, double> hilomodes = modes(bgperc_all_flat);
        pair<double, double> ab_bgl = beta_moments(muvar_bgl.first, muvar_bgl.second);
        pair<double, double> ab_bgh = beta_moments(muvar_bgh.first, muvar_bgh.second);
        
        mixtureDist lo("beta", ab_bgl.first, ab_bgl.second);
        mixtureDist hi("beta", ab_bgh.first, ab_bgh.second);

        vector<mixtureDist> dists;
        dists.push_back(lo);
        dists.push_back(hi);
        vector<double> w{ 0.5, 0.5};
        bgpercmod.init(dists, w);
        bgpercmod.fit(bgperc_all);
        
        optimML::brent_solver intpt(ll_beta_diff);
        intpt.set_root();
        intpt.add_data_fixed("a1", bgpercmod.dists[0].params[0][0]);
        intpt.add_data_fixed("b1", bgpercmod.dists[0].params[0][1]);
        intpt.add_data_fixed("a2", bgpercmod.dists[1].params[0][0]);
        intpt.add_data_fixed("b2", bgpercmod.dists[1].params[0][1]);
        double mean1 = bgpercmod.dists[0].params[0][0] / 
            (bgpercmod.dists[0].params[0][0] + bgpercmod.dists[0].params[0][1]);
        double mean2 = bgpercmod.dists[1].params[0][0] / 
            (bgpercmod.dists[1].params[0][0] + bgpercmod.dists[0].params[0][1]);
        
        double v1 = bgpercmod.dists[0].params[0][0] + bgpercmod.dists[0].params[0][1];
        double v2 = bgpercmod.dists[1].params[0][0] + bgpercmod.dists[1].params[0][1];

        fprintf(stderr, "Background dists (mu, v): (%f, %f) (%f, %f)\n",
            mean1, v1, mean2, v2);

        // Find intersection point between the two dists; use as cutoff.
        return intpt.solve(mean1, mean2);
    }
    else{
        fprintf(stderr, "Disabled filtering on ambient fraction; distribution not bimodal\n");
        return -1.0;
    }
}

void filt_sgrna(vector<vector<double> >& obs,
    vector<set<int> >& results_all,
    vector<double>& llrs_all,
    vector<double>& bg_all,
    vector<double>& himeans,
    vector<double>& bgprops,
    vector<unsigned long>& bcs,
    vector<string>& labels){
    
    set<int> obsmod;

    vector<int> models_kept;
    for (int i = 0; i < himeans.size(); ++i){
        if (himeans[i] > 0){
            models_kept.push_back(i);
        }
    }
    for (int i = 0; i < models_kept.size(); ++i){
        int mod = models_kept[i];
        vector<double> count_on;
        vector<double> count_off;
        vector<int> ons;
        for (int x = 0; x < obs.size(); ++x){
            if (results_all[x].find(mod) != results_all[x].end()){
                count_on.push_back(obs[x][i]);
                ons.push_back(x);
            }
            else{
                count_off.push_back(obs[x][i]);
            }
        }
        int nrm = 0;
        pair<double, double> mv_on = welford(count_on);
        pair<double, double> mv_off = welford(count_off);
        pair<double, double> mp_on = nbinom_moments(mv_on.first, mv_on.second);
        pair<double, double> mp_off = nbinom_moments(mv_off.first, mv_off.second);
        for (int x = 0; x < ons.size(); ++x){
            int idx = ons[x];
            double ll_off = dnbinom(obs[idx][mod], mp_off.first, mp_off.second);
            double ll_on = dnbinom(obs[idx][mod], mp_on.first, mp_on.second);
            if (ll_off > ll_on){
                // Remove.
                nrm++;
                obsmod.insert(idx);
                results_all[idx].erase(results_all[idx].find(mod));
                llrs_all[idx] += ll_off - ll_on;
            }
        }
        //fprintf(stderr, "%s: (%f %f) (%f %f)\tremoved %d of %ld\n", labels[mod].c_str(), mp_off.first,
        //    mp_off.second, mp_on.first, mp_on.second, nrm, ons.size());
    }
    // Recalculate background percent for those modified
    for (set<int>::iterator o = obsmod.begin(); o != obsmod.end(); ++o){
        bg_all[*o] = get_cell_bg(obs[*o], results_all[*o], himeans, bgprops, himeans.size()); 
    }
}

/**
 * Main function for assigning cell identities from data.
 */
void assign_ids(robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    robin_hood::unordered_map<unsigned long, double>& cell_bg,
    vector<string>& labels,
    string& output_prefix, 
    bool sgrna,
    bool filt,
    double prob,
    int nthreads){

    
    // For initial fitting: let counts for each label be a vector 
    vector<vector<double> > obscol;

    vector<double> meanscol;
    
    for (int i = 0; i < labels.size(); ++i){
        vector<double> v;
        obscol.push_back(v);
        meanscol.push_back(0.0);
    }
    
    // Get data into a format needed by stuff.
    vector<vector<double> > obs;
    vector<unsigned long> bcs;
    vector<double> tots;
    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator x = 
        bc_tag_counts.begin(); x != bc_tag_counts.end(); ++x){    
        vector<double> row;
        for (int i = 0; i < labels.size(); ++i){
            row.push_back(0.0);
            obscol[i].push_back(1.0);
        }
        double tot = 0.0;
        for (map<int, long int>::iterator y = x->second.begin(); y != x->second.end(); 
            ++y){
            row[y->first] += (double)y->second;
            tot += (double)y->second;
            double val = (double)y->second;
            //obscol[y->first].push_back(val);
            obscol[y->first][obscol[y->first].size()-1] += val;
            meanscol[y->first] += val;
        }
        // Skip rows with 0 counts
        if (tot > 0.0){
            tots.push_back(tot);
            row.push_back(tot);
            obs.push_back(row);
            bcs.push_back(x->first);
        }
    }
    
    // Step 1: fit a 2-way mixture model to each label independently.
    // This will allow us to get a rough idea of which cells are which
    // identities, and learn about the background distribution.
    vector<double> lomeans;
    vector<double> lophis;
    vector<double> himeans;
    fit_dists_2way(obscol, labels, lomeans, lophis, himeans, output_prefix);
    
    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] > 1.0){
            lomeans[i] -= 1.0;
            himeans[i] -= 1.0;
        }
    }
    
    fprintf(stderr, "Making preliminary assignments...\n");
    
    int nsamp = 1000;
    vector<set<int> > results_all;
    vector<double> llrs_all;
    vector<double> bg_all;
    vector<double> bgprops;
    
    bool twopass = sgrna;
    twopass = false;

    double loglik = assign_cells(obs,
        bcs,
        nsamp,
        labels,
        lomeans,
        lophis,
        himeans,
        sgrna,
        results_all,
        llrs_all,
        bg_all,
        filt || sgrna,
        nthreads,
        bgprops,
        twopass);
    
    if (false){
    //if (sgrna){
        filt_sgrna(obs, results_all, llrs_all, bg_all, lomeans, bgprops, bcs, labels);
    }

    double bg_cutoff = -1.0;
    if (false){
    //if (sgrna){
        bg_cutoff = sep_bg_dists_percent(bg_all);
        fprintf(stderr, "Remove assignments with %%BG > %.3f\n", bg_cutoff*100);
    } 
    for (int i = 0; i < obs.size(); ++i){
        if (bg_cutoff <= 0 || bg_all[i] <= bg_cutoff){
            assn.emplace(bcs[i], results_all[i]);
            assn_llr.emplace(bcs[i], llrs_all[i]);
            cell_bg.emplace(bcs[i], bg_all[i]);
        }
        else{
            cell_bg.emplace(bcs[i], bg_all[i]);
        }
    }
}

/**
 * Spill cell -> identity assignments to disk (in .assignments file)
 */
void write_assignments(robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    string filename, 
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    vector<string>& labels,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    string sep,
    const string& batch_id,
    bool cellranger,
    bool seurat,
    bool underscore,
    bool sgrna){
    
    FILE* f_out = fopen(filename.c_str(), "w");
    
    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator b = bc_tag_counts.begin();
        b != bc_tag_counts.end(); ++b){
        if (!sgrna && (assn.count(b->first) == 0 || assn[b->first].size() == 0)){
            // Filtered/background. Omit.
            continue;
        }
        else{
            char type = 'S';
            set<int>* a = &assn[b->first];
            if (a->size() == 2){
                type = 'D';
            }
            else if (a->size() > 2){
                type = 'M';
            }
            string name = "";
            set<string> names;
            for (set<int>::iterator i = a->begin(); i != a->end(); ++i){
                names.insert(labels[*i]);
            }
            bool first = true;
            for (set<string>::iterator n = names.begin(); n != names.end(); ++n){
                if (!first){
                    name += sep;
                }
                name += *n;
                first = false;
            }
            if (sgrna && a->size() == 0){
                name = "WT";
            }
            string bc_str = bc2str(b->first);
            mod_bc_libname(bc_str, batch_id, cellranger, seurat, underscore);
            if (sgrna){
                fprintf(f_out, "%s\t%s\t%ld\t%f\n", bc_str.c_str(), name.c_str(), a->size(), 
                    assn_llr[b->first]);
            }
            else{
                fprintf(f_out, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type, assn_llr[b->first]); 
            }
        }
    }
    fclose(f_out);
}

void write_assignments_table(robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    string filename,
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    vector<string>& labels,
    const string& batch_id,
    bool cellranger,
    bool seurat,
    bool underscore){
    
    // Put labels in alphabetical order
    vector<pair<string, int> > labelsort;
    for (int i = 0; i < labels.size(); ++i){
        labelsort.push_back(make_pair(labels[i], i));
    }
    sort(labelsort.begin(), labelsort.end());

    FILE* f_out = fopen(filename.c_str(), "w");

    // Create header
    fprintf(f_out, "barcode");
    for (vector<pair<string, int> >::iterator l = labelsort.begin(); l != labelsort.end(); ++l){
        fprintf(f_out, "\t%s", l->first.c_str());
    }
    fprintf(f_out, "\n");

    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator b = bc_tag_counts.begin();
        b != bc_tag_counts.end(); ++b){
        
        string bcstr = bc2str(b->first);
        mod_bc_libname(bcstr, batch_id, cellranger, seurat, underscore);
        
        set<int>* a = &assn[b->first];
        fprintf(f_out, "%s", bcstr.c_str());
        for (vector<pair<string, int> >::iterator l = labelsort.begin(); l != labelsort.end(); ++l){
            if (a->find(l->second) != a->end()){
                fprintf(f_out, "\t1");
            }
            else{
                fprintf(f_out, "\t0");
            }
        }
        fprintf(f_out, "\n");
    }
    fclose(f_out);
}

void write_bg(string filename,
    robin_hood::unordered_map<unsigned long, double>& cell_bg,
    const string& batch_id){
    
    FILE* outf = fopen(filename.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, double>::iterator b = cell_bg.begin();
        b != cell_bg.end(); ++b){
        string bcstr = bc2str(b->first);
        if (batch_id != ""){
            bcstr += "-" + batch_id;
        }
        fprintf(outf, "%s\t%f\n", bcstr.c_str(), b->second);
    }
    fclose(outf);    
}

/**
 * Checks to see whether the user may have specified the incorrect wells - 
 * are any unspecified wells overrepresented in data?
 */
void check_missing_ms_wells(robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >& bc_ms_umis,
    vector<string>& ms_wells,
    map<string, string>& well2name,
    string& outfilename){
    
    // Count each well across all cells
    vector<pair<int, string> > counts;
    for (int i = 0; i < ms_wells.size(); ++i){
        counts.push_back(make_pair(0, ms_wells[i]));
    }
    for (robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >::iterator x = 
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){
        for (int i = 0; i < x->second.size(); ++i){
            if (x->second[i] != NULL){
                counts[i].first += -x->second[i]->count();
            }
        }
    }
    sort(counts.begin(), counts.end());
    
    int minrank_missing = -1;
    int toprank_present = -1;
    FILE* outf = fopen(outfilename.c_str(), "w");
    for (int i = 0; i < counts.size(); ++i){
        if (well2name.count(counts[i].second) == 0 && well2name[counts[i].second] != ""){
            if (minrank_missing == -1){
                minrank_missing = i;
            }
        }
        else{
            toprank_present = i;
        }
        fprintf(outf, "%d\t%s\t%s\n", -counts[i].first, counts[i].second.c_str(),
            well2name[counts[i].second].c_str());
    }
    fclose(outf);
    if (minrank_missing != -1 && minrank_missing < toprank_present){
        fprintf(stderr, "WARNING: one or more un-named wells has higher counts than one or \
more named/expected wells. You may have mis-specified well-to-ID mappings in provided \
metadata. Please see the output well counts file to evaluate.\n");
    }
}

void count_tags_in_reads(vector<string>& read1fn,
    vector<string>& read2fn,
    vector<string>& seqlist,
    int mismatches,
    string& wlfn,
    int umi_len,
    bool sgrna,
    bool exact_cell_barcodes,
    int seq_len,
    bool has_cell_barcodes,
    set<unsigned long>& cell_barcodes,
    string& cell_barcodesfn,
    map<string, int>& seq2idx,
    robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >& bc_tag_umis){
    // Initiate tag barcode whitelist
    seq_fuzzy_match tagmatch(seqlist, mismatches, true, true);
    map<int, int> matchpos_found;
    int matchcount = 0;
    int matchpos_global = -1;

    // How many matches do we need to find before we stop searching through reads?
    int init_searches = 100;

    if (sgrna){
        tagmatch.set_reverse();
    }
    // Set up object that will scan each pair of read files
    bc_scanner scanner;

    // Arguments: 
    // whitelist file
    // whitelist file 2 (for ATAC)
    // file index of cell barcode
    // cell barcode at end of sequence?
    // cell barcode reverse complemented?
    // is there a second whitelist?
    // file index of UMI
    // start index of UMI
    // UMI length
    // cell barcode length
    int bclen = (int)round((double)BC_LENX2/2.0);
    string wl_to_load = wlfn;
    if (cell_barcodes.size() > 0 && cell_barcodesfn != ""){
        // Load the filtered barcode list instead (much faster)
        wl_to_load = cell_barcodesfn;
    }
    scanner.init(wl_to_load, "", 0, false, false, false, 0, bclen, umi_len, bclen);  
    //scanner.init_multiseq_v3(wlfn);
    if (exact_cell_barcodes){
        scanner.exact_matches_only();
    }

    char seq_buf[seq_len+1];
    for (int i = 0; i < read1fn.size(); ++i){
        // Initiate object to read through FASTQs and find cell barcodes.
        scanner.add_reads(read1fn[i], read2fn[i]);
        while (scanner.next()){
            // At this point, there's a valid cell barcode.
            // scanner.barcode_read holds R1
            // scanner.read_f holds R2
            // scanner.umi holds UMI 
            if (scanner.has_umi){
                int idx = -1;
                if (sgrna){
                    
                    if (matchpos_global == -1){
                        // Search entire read for sequence.
                        // Sequences are assumed to appear in forward orientation
                        idx = tagmatch.match(scanner.read_f);
                        if (idx != -1){
                            if (matchpos_found.count(tagmatch.match_pos) == 0){
                                matchpos_found.insert(make_pair(tagmatch.match_pos, 0));
                            }
                            matchpos_found[tagmatch.match_pos]++;
                            matchcount++;

                            if (matchcount >= init_searches){
                                int matchmax = -1;
                                int posmax = -1;
                                bool tie = false;
                                for (map<int, int>::iterator m = matchpos_found.begin();
                                    m != matchpos_found.end(); ++m){
                                    if (matchmax == -1 || m->second > matchmax){
                                        matchmax = m->second;
                                        posmax = m->first;
                                        tie = false;
                                    }
                                    else if (matchmax != -1 && m->second == matchmax){
                                        tie = true;
                                    }
                                }
                                if (posmax != -1 && !tie){
                                    // Set this position to the place to look for future matches.
                                    matchpos_global = posmax;
                                    tagmatch.set_global();
                                }
                            }
                        }
                    }
                    else{
                        // We already know where to look in the read.
                        strncpy(&seq_buf[0], scanner.read_f + matchpos_global, seq_len);
                        seq_buf[seq_len] = '\0';
                        idx = tagmatch.match(&seq_buf[0]);
                    }
                }
                else{
                    // Cell hashing barcodes are assumed to occur at the beginning of R2, forward orientation
                    strncpy(&seq_buf[0], scanner.read_f, seq_len);
                    seq_buf[seq_len] = '\0';
                    idx = tagmatch.match(&seq_buf[0]);
                }
                
                if (idx != -1){
                    // A matching sequence was found.
                    string seqmatch = seqlist[idx];
                    // Store UMI
                    umi this_umi(scanner.umi, scanner.umi_len);    
                    if (bc_tag_umis.count(scanner.barcode) == 0){
                        vector<umi_set_exact*> v(seqlist.size(), NULL);
                        bc_tag_umis.emplace(scanner.barcode, v);
                    }
                    umi_set_exact* ptr = bc_tag_umis[scanner.barcode][seq2idx[seqmatch]];
                    if (ptr == NULL){
                        // Initialize.
                        //bc_tag_umis[scanner.barcode][seq2idx[seqmatch]] = new umi_set_exact(scanner.umi_len);
                        bc_tag_umis[scanner.barcode][seq2idx[seqmatch]] = new umi_set_exact();
                        ptr = bc_tag_umis[scanner.barcode][seq2idx[seqmatch]];
                    }
                    ptr->add(this_umi); 
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {    
    
    srand(time(NULL));

    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"seqs", required_argument, 0, 's'},
       {"names", required_argument, 0, 'N'},
       {"whitelist", required_argument, 0, 'w'},
       {"cell_barcodes", required_argument, 0, 'B'},
       {"exact_only", no_argument, 0, 'e'},
       {"mismatches", required_argument, 0, 'm'},
       {"filt", no_argument, 0, 'f'},
       {"comma", no_argument, 0, 'c'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"exact", no_argument, 0, 'e'},
       {"umi_len", required_argument, 0, 'u'},
       {"sgRNA", no_argument, 0, 'g'},
       {"mtx", required_argument, 0, 'M'},
       {"features", required_argument, 0, 'F'},
       {"feature_type", required_argument, 0, 't'},
       {"prob", required_argument, 0, 'p'},
       {"num_threads", required_argument, 0, 'T'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string output_prefix = "";
    vector<string> read1fn;
    vector<string> read2fn;
    string mapfn = "";
    string wlfn = "";
    string cell_barcodesfn = "";
    int mismatches = 2;
    string sep = "+";
    string batch_id = "";
    bool sgrna = false;
    bool exact_cell_barcodes = false;
    int umi_len = 12;
    string input_mtx = "";
    string input_features = "";
    string featuretype = "";
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;
    double prob = 0.5;
    bool filt = false;
    int nthreads = 1;

    string pr = STRINGIZE(PROJ_ROOT);
    if (pr[pr.length()-1] != '/'){
        pr += "/";
    }
    string barcodesfn = pr + "data/multiseq_indices.txt";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:n:i:M:F:t:1:2:m:N:w:u:B:s:p:T:fCSUegch", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'M':
                input_mtx = optarg;
                break;
            case 'f':
                filt = true;
                break;
            case 'F':
                input_features = optarg;
                break;
            case 't':
                featuretype = optarg;
                break;
            case 's':
                barcodesfn = optarg;
                break;
            case 'n':
                batch_id = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'p':
                prob = atof(optarg);
                break;
            case 'S':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case '1':
                read1fn.push_back(optarg);
               break;
            case '2':
                read2fn.push_back(optarg);
                break;
            case 'N':
                mapfn = optarg;
                break;
            case 'm':
                mismatches = atoi(optarg);
                break;
            case 'w':
                wlfn = optarg;
                break;
            case 'u':
                umi_len = atoi(optarg);
                break;
            case 'B':
                cell_barcodesfn = optarg;
                break;
            case 'c':
                sep = ",";
                break;
            case 'g':
                sgrna = true;
                break;
            case 'e':
                exact_cell_barcodes = true;
                break;
            case 'T':
                nthreads = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: --output_prefix/-o is required\n");
        exit(1);
    }
    if (is_dir(output_prefix)){
        fprintf(stderr, "ERROR: output_prefix %s is a directory, but should be a file \
name prefix.\n", output_prefix.c_str());
        exit(1);
    }
    if (mismatches < -1){
        fprintf(stderr, "ERROR: mismatches must be >= 0, or -1 to allow all matches\n");
        exit(1);
    }
    if (umi_len < 1 || umi_len > 16){
        fprintf(stderr, "ERROR: UMIs must be positive length, up to 16 bp\n");
    }
    if (read1fn.size() > 0 && read2fn.size() > 0 && input_mtx != ""){
        fprintf(stderr, "ERROR: reads provided along with data in MEX format.\n");
        fprintf(stderr, "If loading MEX-format data there is no need to count tags in reads.\n");
        exit(1);
    }
    if (input_features != "" || input_mtx != ""){
        if (input_features == "" || input_mtx == "" || cell_barcodesfn == ""){
            fprintf(stderr, "ERROR: if loading counts from market exchange (MEX) format files,\n");
            fprintf(stderr, "all three of --mtx/-X, --features/-F, and --cell_barcodes/-B are required.\n");
            exit(1);
        }
        else if (!file_exists(input_mtx)){
            fprintf(stderr, "ERROR: file %s does not exist.\n", input_mtx.c_str());
            exit(1);
        }
        else if (!file_exists(input_features)){
            fprintf(stderr, "ERROR: file %s does not exist.\n", input_features.c_str());
            exit(1);
        }
        else if (!file_exists(cell_barcodesfn)){
            fprintf(stderr, "ERROR: file %s does not exist.\n", cell_barcodesfn.c_str());
            exit(1);
        }
    }
    if (sgrna && sep == "+"){
        sep = ",";
    }
    if (!sgrna && prob != 0.5){
        fprintf(stderr, "ERROR: --prob/-p argument not used in --sgRNA mode.\n");
        exit(1);
    } 
    if (prob < 0.5 || prob >= 1.0){
        fprintf(stderr, "ERROR: --prob/-p must be between 0.5 and 1, exclusive.\n");
        exit(1);
    }

    // Data structure to store counts
    //robin_hood::unordered_map<unsigned long, vector<pair<int, string> > > bc_tag_counts;
    robin_hood::unordered_map<unsigned long, map<int, long int> > bc_tag_counts;

    // How many possible labels are there?
    int n_labels;
    vector<string> labels;

    bool has_mapfile = false;

    string countsfilename = output_prefix + ".counts";
    if (file_exists(countsfilename)){
        if (input_mtx != ""){
            fprintf(stderr, "ERROR: %s already exists, but loading counts from %s.\n", 
                countsfilename.c_str(), input_mtx.c_str());
            fprintf(stderr, "Please run again with a different --output_prefix, or load these counts\n");
            fprintf(stderr, "  instead of the MEX-format data.\n");
            exit(1);
        }

        fprintf(stderr, "Loading previously-computed counts (delete file to avoid this \
next time)...\n");
        n_labels = load_counts(countsfilename, bc_tag_counts, labels);
        
        if (cell_barcodesfn != ""){
            // Filter counts.
            set<unsigned long> cell_barcodes;
            parse_barcode_file(cell_barcodesfn, cell_barcodes);
            vector<unsigned long> rm;
            for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator btc = 
                bc_tag_counts.begin(); btc != bc_tag_counts.end(); ++btc){
                if (cell_barcodes.find(btc->first) == cell_barcodes.end()){
                    rm.push_back(btc->first);
                }
            }
            for (vector<unsigned long>::iterator r = rm.begin(); r != rm.end(); ++r){
                bc_tag_counts.erase(*r);
            }
            fprintf(stderr, "%ld cells loaded\n", bc_tag_counts.size());
        }
    }
    else{

        if (input_mtx != ""){
            
            // Load MEX-format data
            if (!parse_mex(cell_barcodesfn, input_features, input_mtx, 
                bc_tag_counts, labels, featuretype)){
                exit(1);
            } 
            
            // Dump MEX-format data to count table and quit
            string countsfn = output_prefix + ".counts";
            dump_counts2(countsfn, labels, bc_tag_counts, batch_id, cellranger, seurat, underscore);
        }
        else{
            if (read1fn.size() == 0 || read2fn.size() == 0){    
                fprintf(stderr, "ERROR: input FASTQs for reads (-1 and -2) are required.\n");
                exit(1);
            }
            else if (read1fn.size() != read2fn.size()){
                fprintf(stderr, "ERROR: unequal number of forward (-1) and reverse (-2) reads provided\n");
                exit(1);
            }
            else if (wlfn == ""){
                fprintf(stderr, "ERROR: cell barcode whitelist / -w is required\n");
                exit(1);
            }
            
            // Load barcode/sgRNA sequence -> well ID (or name) mapping
            map<string, string> seq2well;
            vector<string> seqlist;
            int seq_len = load_seq2well(barcodesfn, seq2well, seqlist);
            
            // Load cell barcodes, if provided
            set<unsigned long> cell_barcodes;
            bool has_cell_barcodes;
            if (cell_barcodesfn != ""){
                has_cell_barcodes = true;
                parse_barcode_file(cell_barcodesfn, cell_barcodes);
            }    
            // Load intermediate ID -> unique identifier mapping, if provided
            map<string, string> well2id;
            if (mapfn == ""){
                for (map<string, string>::iterator sw = seq2well.begin(); 
                    sw != seq2well.end(); ++sw){
                    // If there is no intermediate ID -> final ID mapping for this intermediate ID,
                    // then the intermediate ID is the final ID.
                    if (well2id.count(sw->second) == 0){
                        well2id.insert(make_pair(sw->second, sw->second));
                    }
                }
            }
            else{
                has_mapfile = true;
                load_well_mapping(mapfn, well2id);
                // Make sure wells specified here have associated valid barcodes.
                set<string> wellnames;
                for (map<string, string>::iterator sw = seq2well.begin(); sw != 
                    seq2well.end(); ++sw){
                    wellnames.insert(sw->second);
                }
                for (map<string, string>::iterator wi = well2id.begin(); wi != well2id.end();
                    ++wi){
                    if (wellnames.find(wi->first) == wellnames.end()){
                        fprintf(stderr, "ERROR: well/intermediate ID -> final ID mapping file contains a well named %s\n",
                            wi->first.c_str());
                        fprintf(stderr, "This well ID is not present in the cell hashing/MULTIseq/sgRNA sequence -> \
        well/intermediate ID mappings\n");
                        exit(1);
                    }
                }
                // If the user has omitted one or more wells that turn up often in the data, we will
                // alert the user about this later.
            }
            n_labels = well2id.size();    
            map<string, int> seq2idx;
            vector<string> seq_names;
            vector<string> seq_wells;
            int ix = 0;
            for (map<string, string>::iterator sw = seq2well.begin(); 
                sw != seq2well.end(); ++sw){
            
                seq2idx.insert(make_pair(sw->first, ix));
                ++ix;
                seq_wells.push_back(sw->second);
                if (well2id.count(sw->second) > 0){
                    seq_names.push_back(well2id[sw->second]);
                    labels.push_back(well2id[sw->second]);
                }
                else{
                    seq_names.push_back("");
                }
            }
            // Data structure to store tag counts per cell barcode
            // counts come from UMIs
            robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> > bc_tag_umis;
            // Process reads and count tags in them
            count_tags_in_reads(read1fn, read2fn, seqlist, mismatches, wlfn,
                umi_len, sgrna, exact_cell_barcodes, seq_len, has_cell_barcodes,
                cell_barcodes, cell_barcodesfn, seq2idx, bc_tag_umis);
            // Write counts to disk
            string countsfn = output_prefix + ".counts";
            dump_counts(countsfn, bc_tag_umis, seq_names, bc_tag_counts, batch_id, cellranger, seurat, underscore);
            if (has_mapfile){
                // Check to see that the desired MULTIseq barcodes are the most common ones.
                // Warn the user if unexpected ones are more common. 
                string wellcountsfn = output_prefix + ".wells";
                check_missing_ms_wells(bc_tag_umis, seq_wells, well2id, wellcountsfn);        
            }
            
            // Free stuff
            for (robin_hood::unordered_map<unsigned long, vector<umi_set_exact*> >::iterator x = 
                bc_tag_umis.begin(); x != bc_tag_umis.end(); ){
                for (int i = 0; i < x->second.size(); ++i){
                    if (x->second[i] != NULL){
                        delete x->second[i];
                    }
                }
                bc_tag_umis.erase(x++);
            }
        }
    }
    // Now we can do the actual demultiplexing
    robin_hood::unordered_map<unsigned long, set<int> > assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    robin_hood::unordered_map<unsigned long, double> cell_bg; 
    assign_ids(bc_tag_counts, assn, assn_llr, cell_bg, labels, output_prefix, sgrna, filt, prob, nthreads);
    
    string assnfilename = output_prefix + ".assignments";
    write_assignments(bc_tag_counts, assnfilename, assn, labels, assn_llr, sep, batch_id, 
        cellranger, seurat, underscore, sgrna);
    if (sgrna){
        string tablefilename = output_prefix + ".table";
        write_assignments_table(bc_tag_counts, tablefilename, assn, labels, batch_id, cellranger, seurat,
            underscore);
    }
    string bgfilename = output_prefix + ".bg";
    write_bg(bgfilename, cell_bg, batch_id);
    return 0;  
}
