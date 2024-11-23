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
 * Assigns the likeliest identity, or set of identities, to a cell
 * given data. Can assign anywhere from zero to all possible
 * identities. Restricts the search space by checking for the likeliest
 * individual to include first, then every possible combination including
 * that individual.
 *
 * Also modifies the bgcount parameter to equal the maximum likelihood 
 *  percent background given the assignment times total number of counts
 *  in this cell.
 */
bool assign_cell(vector<double>& bgprops,
    vector<double>& sizes,
    double bg_mean,
    double bg_var,
    vector<double>& obsrow,
    int ndim,
    double slope_mu,
    double intercept_mu,
    double slope_phi,
    double intercept_phi,
    bool add_bg,
    set<int>& result,
    double& llr,
    double& bgcount,
    bool use_fixed_bg_count, 
    bool sgrna,
    double prob){
    
    // Get expected percent background counts from overall count
    double tot = obsrow[obsrow.size()-1];
    double pbg;
    double alpha;
    double beta;
    double phi;

    if (!use_fixed_bg_count || bg_mean == -1 || tot <= bg_mean){
        pbg = expit(slope_mu * log(tot) + intercept_mu);
        phi = exp(slope_phi * log(tot) + intercept_phi);
        alpha = phi*pbg;
        beta = (1-pbg)*phi;
        //pbg = a_b.first * log(tot) + a_b.second;
        //pbg = exp(pbg)/(exp(pbg) + 1);
    }
    else{
        pbg = bg_mean / tot;
    }
    

    vector<pair<double, int> > lls;
    vector<string> names;
    vector<int> nums;
    
    vector<double> obsrowclean;

    // Sort in decreasing order of enrichment of p over expectation in background
    vector<pair<double, int> > dim_enrich_sort;
   
    for (int j = 0; j < ndim; ++j){
        if (sizes[j] > 0){
            double p = obsrow[j]/tot - bgprops[j];
            dim_enrich_sort.push_back(make_pair(-p, j));
            obsrowclean.push_back(obsrow[j]);
        }
    }

    sort(dim_enrich_sort.begin(), dim_enrich_sort.end());
    double tot_counts = 0;
    set<int> included;
    double mean = 0.0;
    double var = 0.0;
    vector<double> llsize;
    
    double llsecond = 0.0;
    double llmax = 0.0;

    double tot_fg = 0.0;
    
    result.clear();
    
    bool skip = false;

    // Should we check for the case of 100% background? 
    if (add_bg){
        vector<double> params_bg;
        for (int j = 0; j < ndim; ++j){
            if (sizes[j] > 0){
                //params_bg.push_back(pbg * bgprops[j]);
                //params_bg.push_back(bgprops[j]);
                
                params_bg.push_back(phi*bgprops[j]);
            }
        }
        double ll = dmultinom(obsrowclean, params_bg);
        //double ll = ddirichletmultinom(obsrowclean, params_bg);

        if (bg_mean != -1 && bg_var != -1){
            double cbg = obsrow[obsrow.size()-1];

            // Use log-normal distribution
            double m2 = bg_mean*bg_mean;
            double sigma = sqrt(log(bg_var/m2 + 1.0));
            double mu = log(bg_mean) - log(bg_var/m2 + 1.0)/2.0;

            //ll += dnorm(log(cbg), mu, sigma);
        }
        //fprintf(stderr, "-1) %f\n", ll);
        llmax = ll;
        skip = true;
    }
    
    double llprev = llmax;

    long int maxn = dim_enrich_sort.size();
    if (!sgrna){
        maxn = 2;
    }
    
    for (int i = 0; i < maxn; ++i){
        int j = dim_enrich_sort[i].second;
        if (sizes[j] <= 0 || obsrow[j] == 0){
            continue;
        }
        else if (-dim_enrich_sort[i].first < 0){
        //    continue;
        } 
        included.insert(j);
        
        tot_counts += sizes[j];
        
        // Assume exponential variance (and use scaling factor to prevent overflow)
        var += pow(sizes[j], 2)/1e6;
        
        tot_fg += obsrow[j] * (1-bgprops[j]);

        //vector<double> p;
        vector<double> a;
        for (int k = 0; k < ndim; ++k){
            if (sizes[k] > 0){
                double a_this;
                //double p_this;
                if (included.find(k) != included.end()){
                    //p_this = bgprops[k] * pbg + (1.0 - pbg) * (sizes[k]/tot_counts);
                    //a_this = alpha*bgprops[k]*pbg + beta*(1.0-pbg) * (sizes[k]/tot_counts);
                    a_this = (bgprops[k] * pbg + (1.0-pbg)*(sizes[k]/tot_counts))*phi;
                }
                else{
                    //p_this = bgprops[k] * pbg;
                    //a_this = alpha*bgprops[k] * pbg;
                    a_this = bgprops[k] * pbg * phi;
                }
                //p.push_back(p_this);
                a.push_back(a_this);
            }
        }
        //double ll = dmultinom(obsrowclean, p);
        double ll = ddirichletmultinom(obsrowclean, a);

        if (bg_mean != -1 && bg_var != -1){
            double tot_exp = tot_counts + bg_mean;
            var += bg_var/1e6;
            
            // Use log-normal distribution
            double m2 = pow(tot_exp, 2)/1e6;
            double sigma = sqrt(log(var/m2 + 1.0));
            double mu = log(tot_exp) - log(var/m2 + 1.0)/2.0;
            //ll += dnorm(log(obsrow[obsrow.size()-1]), mu, sigma);
        }
        
        //fprintf(stderr, "%d) %f\n", j, ll); 
        // Things are sorted so that once log likelihood starts to dip, we can stop looking
        // (it will only continue to decrease).
        if (llmax != 0.0 && ll < llmax){
            if (llsecond == 0.0 || llsecond < ll){
                llsecond = ll;
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
            skip = false;
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
            if (keep && !skip){
                llsecond = llmax;
                llmax = ll;
                llprev = ll;
                result.insert(j);
            }
        }
    }
    
    vector<double> model_means_this = sizes;

    //fprintf(stderr, "\n");

    llr = llmax - llsecond;
    
    // Result should already be populated
    if (result.size() == 0){
        bgcount = obsrow[obsrow.size()-1];
        return false;
    } 
    else{
        bgcount = get_cell_bg(obsrow, result, model_means_this, bgprops, ndim);
        return true;
    }
}

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
    vector<double>& model_means,
    string& output_prefix){
    
    vector<double> loweights;
    vector<double> hiweights;
    
    string outfn = output_prefix + ".dists";
    FILE* outf = fopen(outfn.c_str(), "w");

    fprintf(stderr, "===== Initial model means: =====\n");
    for (int i = 0; i < labels.size(); ++i){
        if (obscol[i].size() < 3){
            lomeans.push_back(0);
            lophis.push_back(0);
            himeans.push_back(0);
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
            
            vector<mixtureDist> dists2;
            mixtureDist low2("negative_binomial", vector<double>{ lowval2, 1.0 });
            mixtureDist high2("exponential", 1.0/highval);
            dists2.push_back(low2);
            dists2.push_back(high2);
            mixtureModel mod2(dists2);
            mod2.fit(obscol[i]);
            /* 
            fprintf(stderr, "LL %s %f %f\n", labels[i].c_str(), mod.loglik, mod2.loglik);
            double gini1 = 0.0;
            for (int x = 0; x < obscol[i].size(); ++x){
                double gini1b = 1.0 - mod.responsibility_matrix[x][0]*mod.responsibility_matrix[x][0] - 
                    mod.responsibility_matrix[x][1]*mod.responsibility_matrix[x][1];
                gini1 += gini1b;
            }
            gini1 /= (double)obscol[i].size();

            double gini2 = 0.0;
            for (int x = 0; x < obscol[i].size(); ++x){
                double gini2b = 1.0 - mod2.responsibility_matrix[x][0]*mod2.responsibility_matrix[x][0] - 
                    mod2.responsibility_matrix[x][1]*mod2.responsibility_matrix[x][1];
                gini2 += gini2b;
            }
            gini2 /= (double)obscol[i].size();
            fprintf(stderr, "  gini %f %f\n", gini1, gini2);
            *
            fprintf(stderr, "  means %f %f | %f %f\n", mod.dists[0].params[0][0], 1.0/mod.dists[1].params[0][0],
                mod2.dists[0].params[0][0], 1.0/mod2.dists[1].params[0][0]);
            fprintf(stderr, "  weights %f | %f\n", mod.weights[0], mod2.weights[0]);
            */
            double lomean = mod.dists[0].params[0][0];
            double lophi = mod.dists[0].params[0][1];
            double himean = 1.0/mod.dists[1].params[0][0];
            double w = mod.weights[0];
           // fprintf(stderr, "  disps %f %f\n", mod.dists[0].params[0][1], mod2.dists[0].params[0][1]);
            if (1.0/mod2.dists[1].params[0][0] > mod2.dists[0].params[0][0] && 
                mod2.loglik > mod.loglik){
                //fprintf(stderr, "  SWAP\n");
                lomean = mod2.dists[0].params[0][0];
                lophi = mod2.dists[0].params[0][1];
                himean = 1.0/mod2.dists[1].params[0][0];
                w = mod2.weights[0];
            }

            double dist_sep = pnbinom(himean, lomean, lophi);
            fprintf(stderr, "%s:\t%f\t%f\n", labels[i].c_str(), lomean, himean);
            fprintf(outf, "%s\t%f\t%f\t%f\t%f\n", labels[i].c_str(),
                w, lomean, lophi, 1.0/himean);

            if (dist_sep < 0.95){
            //if (mod.weights[1] >= mod.weights[0] || dist_sep < 0.99){
                lomeans.push_back(0);
                lophis.push_back(0);
                himeans.push_back(0);
                hiweights.push_back(0);
                loweights.push_back(1);
            }
            else{
                lomeans.push_back(lomean);
                lophis.push_back(lophi);
                himeans.push_back(himean);
                hiweights.push_back(1.0-w);
                loweights.push_back(w);
            }
        }
    }
    
    fclose(outf);

    // Check to see whether any label appears to have dropped out.
    // We fit normal distributions to the means of low and high count
    // distributions across the set of all labels. If any label's high
    // count mean looks more like a typical low-count than high-count
    // mean, then that label will be removed.
    
    fprintf(stderr, "Checking for label dropout...\n");
    pair<double, double> lomv = welford(lomeans);
    pair<double, double> himv = welford(himeans);
    
    for (int i = 0; i < labels.size(); ++i){
        
        string filtstr = "";
        if (lophis[i] == 0.0 || himeans[i] == 0.0){
            fprintf(stderr, "  Remove label %s, insufficient dist separation\n", 
                labels[i].c_str());
            model_means.push_back(-1);
            lomeans[i] = -1;
            lophis[i] = -1;
            himeans[i] = -1;    
        }
        else if (false){
            double d1 = dnorm(himeans[i], lomv.first, sqrt(lomv.second));
            double d2 = dnorm(himeans[i], himv.first, sqrt(himv.second));
            if (d2 - d1 < 0){
                fprintf(stderr, "  Remove label %s, LLR(ambient-foreground) = %f\n", 
                    labels[i].c_str(),d1-d2);
                model_means.push_back(-1);
                lomeans[i] = -1;
                lophis[i] = -1;
                himeans[i] = -1;
            }
            else{
                model_means.push_back(himeans[i]);
            }
        }
    }
}

double ll_beta_test2(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
   
    double tot = data_d.at("tot");
    double pbg = data_d.at("pbg");
    
    double mu = expit(params[0] * log(tot) + params[1]);
    double phi = exp(params[2] * log(tot) + params[3]);
    //fprintf(stderr, "mu %f phi %f\n", mu, phi);

    double alpha = mu*phi;
    double beta = (1.0-mu)*phi;
    return dbeta(pbg, alpha, beta)/log2(exp(1));
    //return phi*(mu*log(pbg) + log(1-pbg) - mu*log(1-pbg)) - log(pbg) - log(1-pbg) - lgamma(mu*phi) - lgamma(phi-mu*phi) + lgamma(phi);
    //return (alpha-1)*log(pbg) + (beta-1)*log(1-pbg) - lgamma(alpha) - lgamma(beta) + lgamma(alpha+beta);
}

void dll_beta_test2(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double tot = data_d.at("tot");
    double pbg = data_d.at("pbg");
    
    double mu_x = params[0] * log(tot) + params[1];
    double phi_x = params[2] * log(tot) + params[3];

    double mu = expit(mu_x);
    double phi = exp(phi_x);

    double dmu_dx = exp(mu_x)/((exp(mu_x) + 1)*(exp(mu_x) + 1));
    double dphi_dx = exp(phi_x);

    double dx_dp0 = log(tot);
    double dx_dp1 = 1;
    double dx_dp2 = log(tot);
    double dx_dp3 = 1;
    
    double dll_dmu = phi*log(pbg) - phi*log(1-pbg) - phi*digamma(mu*phi) + phi*digamma(phi - mu*phi);
    double dll_dphi = mu*log(pbg) + log(1-pbg) - mu*log(1-pbg) - mu*digamma(mu*phi) - 
        (1-mu)*digamma(phi-mu*phi) + digamma(phi);
    
    //fprintf(stderr, "dll_dmu %f dmu_dx %f dx_dp0 %f | dll_dphi %f dphi_dx %f dx_dp2 %f\n",
    //    dll_dmu, dmu_dx, dx_dp0, dll_dphi, dphi_dx, dx_dp2);

    results[0] += dll_dmu * dmu_dx * dx_dp0;
    results[1] += dll_dmu * dmu_dx * dx_dp1;
    results[2] += dll_dphi * dphi_dx * dx_dp2;
    results[3] += dll_dphi * dphi_dx * dx_dp3;
    
    /*
    double mu = expit(params[0] * log(tot) + params[1]);
    double phi = exp(params[2] * tot + params[3]);
    
    double dmu_dx = exp(
    //double dmu_dt = exp(mu)/((exp(mu) + 1)*(exp(mu) + 1));
    //double dt_dm = log(tot);
    //mu = expit(mu);
    //results[0] += (phi*log(pbg) - phi*log(1-pbg) - phi*digamma(mu*phi) + phi*digamma(phi-phi*mu))*dmu_dt*dt_dm;
    results[1] += (mu*log(pbg) + log(1-pbg) - mu*log(1-pbg) - mu*digamma(mu*phi) - (1-mu)*digamma(phi - mu*phi) + 
           digamma(phi))*tot; 
    results[0] += (log(pbg) - digamma(alpha) + digamma(alpha+beta))*tot;
    //results[1] += (log(1-pbg) - digamma(beta) + digamma(alpha+beta))*tot;
    */
}
double ll_beta_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    
    double tot = data_d.at("tot");
    double pbg = data_d.at("pbg");
    double a = params[0];
    double b = params[1];
    double alpha = a*tot;
    double beta = b*tot;
    return (alpha-1)*log(pbg) + (beta-1)*log(1-pbg) - lgamma(alpha) - lgamma(beta) + lgamma(alpha+beta);
}

void dll_beta_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double tot = data_d.at("tot");
    double pbg = data_d.at("pbg");
    double a = params[0];
    double b = params[1];
    double alpha = a*tot;
    double beta = b*tot;
    results[0] += (log(pbg) - digamma(alpha) + digamma(alpha+beta))*tot;
    results[1] += (log(1-pbg) - digamma(beta) + digamma(alpha+beta))*tot;
}

double ll_binom_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double tot = data_d.at("tot");
    double bg = data_d.at("bg");
    double p = expit(params[0] * log(tot) + params[1]);
    return dbinom(tot, bg, p)/log2(exp(1));
}

void dll_binom_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double tot = data_d.at("tot");
    double bg = data_d.at("bg");
    double p = expit(params[0] * log(tot) + params[1]);
    double dll_dp = (bg - tot*p)/(p - p*p);
    double x = params[0] * log(tot) + params[1];
    double dp_dx = exp(x)/((exp(x) + 1)*(exp(x) + 1));
    double dx_dp0 = log(tot);
    double dx_dp1 = 1.0;
    results[0] += dll_dp * dp_dx * dx_dp0;
    results[1] += dll_dp * dp_dx * dx_dp1;
}

void fit_bg_predictors(vector<double>& logtots,
    vector<double>& logitpercs,
    vector<double>& tots,
    vector<double>& percs,
    double& slope_mu,
    double& intercept_mu,
    double& slope_phi,
    double& intercept_phi){
    
    // Goal: learn how to predict % ambient tags (background) in a cell given its
    // total tag counts.

    // We are going to do a GLM fitting-like procedure to estimate parameters.
    
    // Model: Beta distribution, parameterized as mu = a/(a+b) and phi = a+b.
    // mu is predicted by expit(b_0*x + b_1), where x is the log total tag count.
    // phi is predicted by exp(b_2*x + b_2), where x is the log total tag count.
    
    // To fit the model, we need initial guesses for b_0, b_1, b_2, and b_3.
    // We will estimate these via linear regression.
   
    // First, get b_0 & b_1. These are parameters for estimating logit(mu) from log(total counts)
    pair<double, double> slope_int = linreg(logtots, logitpercs);

    // Next, get b_2 and b_3. These are parameters for estimating log(phi) from log(total counts).
    vector<double> phipred_x;
    vector<double> phipred;
    vector<double> phipred_w;
    
    for (int i = 0; i < logtots.size(); ++i){
        double pred = expit(slope_int.first * logtots[i] + slope_int.second);
        double obs = expit(logitpercs[i]);
        
        // Method of Moments estimator for phi is (mu*(1-mu))/var - 1
        // We already have mu estimate from previous linear regression.
        // Estimate var here as (observation - mu)^2
        double p = (pred*(1.0-pred))/pow(pred-obs, 2) - 1.0;
        // Only use values we can log-transform without issue.
        if (!isnan(log(p))){
            phipred.push_back(log(p));
            phipred_x.push_back(logtots[i]);
            phipred_w.push_back(pow(exp(logtots[i]), 2));
        }
    }
    
    pair<double, double> slope_int_phi = linreg(phipred_x, phipred);
     
    if (isnan(slope_int.first) || isnan(slope_int.second) || isnan(slope_int_phi.first) || isnan(slope_int_phi.second)){
        fprintf(stderr, "HERE\n");
        fprintf(stderr, "%f %f | %f %f\n", slope_int.first, slope_int.second, slope_int_phi.first, slope_int_phi.second);
        exit(1);
    } 
    
    /*
    vector<double> bgs;
    for (int i = 0; i < tots.size(); ++i){
        bgs.push_back(percs[i] * tots[i]);
    } 
    vector<double> params0{ slope_int.first, slope_int.second };
    optimML::multivar_ml_solver solver0(params0, ll_binom_test, dll_binom_test);
    solver0.add_data("tot", tots);
    solver0.add_data("bg", bgs);
    solver0.solve();
    fprintf(stderr, "LL0 %f\n", solver0.log_likelihood);
    fprintf(stderr, "PARAMS %f %f\n", solver0.results[0], solver0.results[1]);
    */

    // We now have estimates of b_0, b_1, b_2, and b_3. Now fit a GLM.
    vector<double> params{ slope_int.first, slope_int.second, slope_int_phi.first, slope_int_phi.second };
    //vector<double> params{ 0.0, 10, 0,0};
    optimML::multivar_ml_solver solver(params, ll_beta_test2, dll_beta_test2);
    solver.add_data("tot", tots);
    solver.add_data("pbg", percs);
   // solver.add_weights(tots);
    solver.solve();
    fprintf(stderr, "LL %f\n", solver.log_likelihood); 
    fprintf(stderr, "PARAMS %f %f %f %f\n", solver.results[0], solver.results[1], solver.results[2], solver.results[3]); 
    // Extract parameters of interest
    slope_mu = solver.results[0];
    intercept_mu = solver.results[1];
    slope_phi = solver.results[2];
    intercept_phi = solver.results[3];
}
/**
 * Do initial assignments to learn the % of background/ambient tags
 * composed of each label, and to find the relationship between
 * log total tag count and logit percent background (returned as
 * linear regression parameters)
 */
void assign_init(vector<vector<double> >& obs,
    vector<string>& labels,
    vector<double>& bgmeans,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    vector<double>& model_means,
    bool sgrna,
    double& slope_mu,
    double& intercept_mu,
    double& slope_phi,
    double& intercept_phi){

    // What does the probability of high-count need to be 
    // for a cell to be assigned that label in the first round?
    double p_thresh = 0.5;
    
    // Figure out mean ambient/background counts per cell
    
    for (int i = 0; i < labels.size(); ++i){
        bgmeans.push_back(0.0);
    }
    
    // This will be used to find the relationship between total
    // counts and percent background in a cell
    vector<double> logtots;
    vector<double> logitpercs;
    vector<double> lrweights;
    vector<unsigned long> logtots_bc;
    double weightfrac = 1.0/(double)obs.size();

    vector<double> tots;
    vector<double> percs;
    vector<double> bgcounts;

    double bgcount = 0.0;

    vector<vector<double> > counts_all;
    vector<vector<double> > weights_all;
    for (int i = 0; i < labels.size(); ++i){
        vector<double> v;
        counts_all.push_back(v);
        weights_all.push_back(v);
    }
    
    // Get initial assignments and use these to learn about background
    // distribution of counts
    double maxtot = 0;

    for (int i = 0; i < obs.size(); ++i){
        set<int> chosen;
        
        double count_tot = obs[i][obs[i].size()-1];
        double countbg = 0.0;
        
        for (int j = 0; j < labels.size(); ++j){
            double p_low = 0.0;
            double p_high = 0.0;
            if (obs[i][j] == 0 || model_means[j] == -1){
                p_low = 1.0;
                p_high = 0.0;
            }
            else{
                double ll_low = dnbinom(obs[i][j], lomeans[j], lophis[j]);
                double ll_high = dexp(obs[i][j], 1.0/himeans[j]);
                double max;
                if (ll_low > ll_high){
                    max = ll_low;
                }
                else{
                    max = ll_high;
                }
                double tot = pow(2, ll_low-max) + pow(2, ll_high-max);
                p_low = pow(2, ll_low-max) / tot;
                p_high = pow(2, ll_high-max) / tot;
            }

            // Also ensure the value is higher than the mean of the lower dist, because sometimes
            // the exponential can scoop up values on the left tail too
            if (p_high > p_thresh && obs[i][j] > lomeans[j]){
                counts_all[j].push_back(obs[i][obs[i].size()-1]);
                chosen.insert(j);
            }
            else{
                countbg += obs[i][j];
                bgmeans[j] += obs[i][j];
            }
        }

        bgcount += countbg;
        
        if (count_tot > 0 && countbg > 0 && countbg < count_tot){
            logtots.push_back(log(count_tot));
            logitpercs.push_back(log(countbg/count_tot) - log(1-countbg/count_tot));
            lrweights.push_back(weightfrac * (count_tot*count_tot));
        
            if (countbg > 0 && countbg < count_tot){
                percs.push_back(countbg/count_tot);     
                tots.push_back((double)count_tot);
                bgcounts.push_back((double) countbg);
                if (count_tot > maxtot){
                    maxtot = count_tot;
                }
                //fprintf(stdout, "%f\t%f\n", countbg, count_tot);
            }   
        }
    }
    
    // Get params
    fit_bg_predictors(logtots, logitpercs, tots, percs, slope_mu,
        intercept_mu, slope_phi, intercept_phi);

    fprintf(stderr, "===== Ambient tag profile: =====\n");
    // Learn the proportions of each label in ambient/background counts
    for (int i = 0; i < labels.size(); ++i){
        bgmeans[i] /= bgcount;
        fprintf(stderr, "%s:\t%f\n", labels[i].c_str(), bgmeans[i]);
    }
    
}

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
    double part1 = lgammaf_r(term1, &signp);
    if (isnan(part1) || isinf(part1) || isnan(term2) || isinf(term2) || isnan(term3) || isinf(term3)){
        fprintf(stderr, "ERR %f %f %f\n", part1, term2, term3);
        return 0.0;
    }
    return lgammaf_r(term1, &signp) - term2 + term3;
}

double ll_dir2(const vector<double>& params,
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
        double a = exp(params[i*3]*ltot*ltot + params[i*3+1]*ltot + params[i*3 + 2]);
        
        sprintf(&buf[0], "p_%d", i);
        bufstr = buf;
        double p = data_d.at(bufstr);
        
        term1 += a;
        term2 += lgammaf_r(a, &signp);
        term3 += (a - 1.0)*log(p);
        
    }
    if (isnan(lgammaf_r(term1, &signp) || isnan(term2) || isnan(term3))){
        fprintf(stderr, "ERR %f %f %f\n", lgammaf_r(term1, &signp), term2, term3);
        exit(1);
    }
    else if (isinf(lgammaf_r(term1, &signp)) || isinf(term2) || isinf(term3)){
        fprintf(stderr, "ERR2 %f %f %f\n", lgammaf_r(term1, &signp), term2, term3);
        for (int i = 0; i < nmod; ++i){
            double a = exp(params[i*3]*ltot*ltot + params[i*3+1]*ltot + params[i*3 + 2]);
            sprintf(&buf[0], "p_%d", i);
            bufstr = buf;
            double p = data_d.at(bufstr);
            fprintf(stderr, "%d) %f %f\n", i, a, p);
        }
        exit(1);
    }
    return lgammaf_r(term1, &signp) - term2 + term3;
}

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
        
        if (isnan(dll_da[dll_da.size()-1]) || isinf(dll_da[dll_da.size()-1]) ||
                isnan(a) || isinf(a) || isnan(term1) || isinf(term1)){
            return;
        }
    }
    for (int i = 0; i < nmod; ++i){
        dll_da[i] += digamma(term1);
        results[i*2] += dll_da[i] * ltot * a_all[i];
        results[i*2 + 1] += dll_da[i] * a_all[i];
    }
}

void dll_dir2(const vector<double>& params,
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
        double a = exp(params[i*3]*ltot*ltot + params[i*3+1]*ltot + params[i*3 + 2]);
        a_all.push_back(a);

        sprintf(&buf[0], "p_%d", i);
        bufstr = buf;
        double p = data_d.at(bufstr);
        term1 += a;
        
        dll_da.push_back(-digamma(a) + log(p));
    }
    for (int i = 0; i < nmod; ++i){
        dll_da[i] += digamma(term1);
        results[i*3] += dll_da[i] * ltot*ltot * a_all[i];
        results[i*3 + 1] += dll_da[i] * ltot * a_all[i];
        results[i*3 + 2] += dll_da[i] * a_all[i];
    }
}

void get_dirichlet(map<int, int>& idx2idx,
    set<int>& tags,
    vector<vector<double> >& bg_samps, 
    vector<double>& bg_tots,
    map<int, vector<double> >& model_samps,
    vector<double>& result,
    double& mu,
    double& phi,
    double& slope,
    double& intercept,
    double& var_slope,
    double& var_intercept,
    double& mean_pbg,
    double& var_pbg,
    vector<double>& dircoeffs,
    vector<double>& betacoeffs){

    int nsamp = bg_samps.size();
    int nmod = idx2idx.size();
    vector<double> means(nmod);
    for (int j = 0; j < nmod; ++j){
        means[j] = 0.0;
    }
    
    vector<double> logtots;
    vector<double> logpercs;
    vector<double> logitpercs;

    vector<double> tots;
    vector<double> bgs;
    vector<double> pbgs;

    vector<vector<double> > samps_frac;
    vector<vector<double> > samps_frac_column;
    vector<vector<double> > samps_frac_column_log;
    vector<double> inset_log;
    vector<double> inset;

    for (int i = 0; i < nmod; ++i){
        vector<double> v;
        samps_frac_column.push_back(v);
        samps_frac_column_log.push_back(v);
    }

    for (int i = 0; i < nsamp; ++i){
        vector<double> row = bg_samps[i];
        double tot = bg_tots[i];
        double totsamp = 0.0;
         
        map<int, double> percmember;
        
        double il = 0.0;
        for (set<int>::iterator member = tags.begin(); member != tags.end(); ++member){
            double this_samp = model_samps[*member][i];
            totsamp += this_samp;
            int ri = idx2idx[*member];
            tot -= row[ri];
            row[ri] = this_samp;
            //row[ri] += this_samp;
            tot += this_samp;
            il += this_samp;
        }
        inset_log.push_back(log(il));
        inset.push_back(il);

        logtots.push_back(log(tot));
        logpercs.push_back(log(1.0-totsamp/tot));
        logitpercs.push_back(log(1.0-totsamp/tot)-log(totsamp/tot));
        bgs.push_back(tot-totsamp);
        pbgs.push_back(1.0-totsamp/tot);
        tots.push_back(tot);
        for (int j = 0; j < nmod; ++j){
            if (row[j] == tot){
                samps_frac_column_log[j].push_back(log((row[j]-1)/tot));
            }
            else{
                samps_frac_column_log[j].push_back(log((row[j] + 1)/(tot+1)));
            }
            row[j] /= tot;
            samps_frac_column[j].push_back(row[j]);
            means[j] += (1.0/(double)nsamp)*row[j];
        }
        
        samps_frac.push_back(row);
    } 
    /*    
    pair<double, double> sitest = linreg(logtots, logitpercs);
    vector<double> params_beta{ sitest.first, sitest.second, 0.0, 0.0 };
    optimML::multivar_ml_solver solv2(params_beta, ll_beta_test2, dll_beta_test2);
    solv2.add_data("tot", tots);
    solv2.add_data("pbg", pbgs);
    solv2.solve();
    betacoeffs = solv2.results;
    */

    vector<double> params_d;
    for (int j = 0; j < nmod; ++j){
        pair<double, double> si = linreg(inset_log, samps_frac_column_log[j]);
        //params_d.push_back(0.0);
        params_d.push_back(si.first);
        params_d.push_back(si.second);
    }
    optimML::multivar_ml_solver solv(params_d, ll_dir, dll_dir);
    solv.add_data_fixed("nmod", nmod);
    char buf1[100];
    string buf1str;
    for (int j = 0; j < nmod; ++j){
        sprintf(&buf1[0], "p_%d", j);
        buf1str = buf1;
        solv.add_data(buf1str, samps_frac_column[j]);
    }
    solv.add_data("tot", inset);
    solv.solve();
    
   /* 
    for (int i = 0; i < tots.size(); ++i){
        vector<double> preds;
        double predtot = 0.0;
        for (int j = 0; j < samps_frac_column.size(); ++j){
            //double pred = exp(logtots[i] * solv.results[j*2] + solv.results[j*2 + 1]);
            //double pred = exp(inset_log[i] * inset_log[i] * solv.results[j*3] + inset_log[i] * solv.results[j*3 + 1] + solv.results[j*3 + 2]);
            double pred = exp(inset_log[i] * solv.results[j*2] + solv.results[j*2+1]);
            predtot += pred;
            preds.push_back(pred);
        }
        for (int j = 0; j < samps_frac_column.size(); ++j){
            fprintf(stdout, "%f\t%d\t%f\t%f\t%f\n", logtots[i], j, samps_frac_column[j][i], preds[j]/predtot, predtot);
            //pair<double, double> mv = linreg(logtots, samps_frac_column_log[j]);
            //fprintf(stderr, "%d) %f %f\n", j, mv.first, mv.second);
            //fprintf(stdout, "%f\t%d\t%f\t%f\n", logtots[i], j, samps_frac_column[j][i], exp(mv.first*logtots[i] + mv.second));
        }
    }
    exit(0);
    */

    dircoeffs = solv.results;

    pair<double, double> muvar_pbg = welford(pbgs);
    mean_pbg = muvar_pbg.first;
    var_pbg = muvar_pbg.second;

    pair<double, double> slope_int = linreg(logtots, logpercs);

    slope = slope_int.first;
    intercept = slope_int.second;
    
    double min = -slope_int.second/slope_int.first;
    vector<double> lt;
    vector<double> var;
    for (int x = 0; x < logtots.size(); ++x){
        if (logtots[x] > min){
            double pred = exp(slope_int.first*logtots[x] + slope_int.second);
            double v = pow(pred - exp(logpercs[x]),2);
            if (!isnan(log(v))){
                lt.push_back(logtots[x]);
                var.push_back(log(v));
            }
        }
    }

    pair<double, double> slope_int_var = linreg(lt, var);
    var_slope = slope_int_var.first;
    var_intercept = slope_int_var.second;

    /* 
    pair<double, double> slint2 = linreg(logtots, logitpercs);
    optimML::multivar_ml_solver solver(vector<double>{ slint2.first, slint2.second, 0.0, 1.0 }, ll_beta_test2, dll_beta_test2);
    solver.add_data("tot", tots);
    solver.add_data("pbg", pbgs);
    solver.solve();
    fprintf(stderr, "BETA %f %f %f %f\n", solver.results[0], solver.results[1], solver.results[2], solver.results[3]);
    */

    fit_dirichlet(means, samps_frac, result, 1);
    pair<double, double> muvar = welford(tots);
    pair<double, double> muphi = nbinom_moments(muvar.first, muvar.second);
    mu = muphi.first;
    phi = muphi.second;
    /* 
    if (tags.size() == 1){
        fprintf(stderr, "[%d]\n", *tags.begin());
        double a0 = 0.0;
        for (int i = 0; i < result.size(); ++i){
            a0 += result[i];
        }
        for (int i = 0; i < result.size(); ++i){
            fprintf(stderr, "%d) %f\n", i, result[i]/a0);
        }
    }
    */
}

void assign_fallback(vector<double>& obsrow,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    int max_assn,
    set<int>& result,
    double& assn_llr){
    
    vector<double> llrs;
    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] > 0 && himeans[i] > 0){
            double ll1 = dnbinom(obsrow[i], lomeans[i], lophis[i]);
            double ll2 = dexp(obsrow[i], 1.0/himeans[i]);
            if (ll2 > ll1){
                result.insert(i);
                llrs.push_back(ll2-ll1);
            }
            else{
                llrs.push_back(ll1-ll2);
            }
            
        }
        else{
            llrs.push_back(0.0);
        }
    }
    if (max_assn > 0){
        if (result.size() > max_assn){
            vector<pair<double, int> > llrsort;
            for (int i = 0; i < llrs.size(); ++i){
                llrsort.push_back(make_pair(llrs[i], i));
            }
            sort(llrsort.begin(), llrsort.end());
            for (int i = 0; i < llrsort.size(); ++i){
                if (result.find(llrsort[i].second) != result.end()){
                    llrs[llrsort[i].second] = -llrs[llrsort[i].second];
                    result.erase(result.find(llrsort[i].second));
                }
                if (result.size() <= max_assn){
                    break;
                }
            }
        }
    }
    assn_llr = 0;
    for (int i = 0; i < llrs.size(); ++i){
        assn_llr += llrs[i];
    }
}

double assign_testnew(vector<vector<double> >& obs,
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
    bool filt){
    
    bool add_bg = false;

    map<modelkey, vector<double> > model_params;
    map<modelkey, pair<double, double> > model_muphi;
    map<modelkey, pair<double, double> > model_slopeint;
    map<modelkey, pair<double, double> > model_varslopeint;
    map<modelkey, pair<double, double> > model_pbgmuvar;
    map<modelkey, vector<double> > model_dircoeffs;
    map<modelkey, vector<double> > model_betacoeffs;

    // Get model for all background
    modelkey key_bg;
    vector<int> models_kept;
    map<int, int> idx2idx;
    for (int j = 0; j < lomeans.size(); ++j){
        fprintf(stderr, "lab %s) %f %f\n", labels[j].c_str(), lomeans[j], himeans[j]);
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
        vector<double> sampfrac = samp;
        for (int j = 0; j < models_kept.size(); ++j){
            sampfrac[j] /= tot;
        }
        samps.push_back(samp);
        samps_frac.push_back(sampfrac);
        rowtots.push_back(tot);
    }
    double meanconctot = 0.0;
    for (int i = 0; i < meanconc.size(); ++i){
        meanconctot += meanconc[i];
    } 
    for (int i = 0; i < meanconc.size(); ++i){
        meanconc[i] /= meanconctot;
    }

    pair<double, double> rowtotmuvar = welford(rowtots);
    pair<double, double> bgmuphi = nbinom_moments(rowtotmuvar.first, rowtotmuvar.second);
    model_muphi.insert(make_pair(key_bg, bgmuphi));
    
    fprintf(stderr, "background count mean %f\n", bgmuphi.first);

    vector<double> dirprops_bg;
    fit_dirichlet(meanconc, samps_frac, dirprops_bg, 1);

    vector<double> bgprops;
    vector<double> bgprops2;

    double bgtot = 0.0;
    for (int j = 0; j < dirprops_bg.size(); ++j){
        bgtot += dirprops_bg[j];
    }
    
    fprintf(stderr, "===== Ambient tag profile: =====\n");
    //for (int j = 0; j < dirprops_bg.size(); ++j){
    for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
        bgprops.push_back(dirprops_bg[ix->second] / bgtot);
        fprintf(stderr, "%s:\t%f\n", labels[ix->first].c_str(), meanconc[ix->second]); 
        //fprintf(stderr, "%s:\t%f\n", labels[ix->first].c_str(), bgprops[ix->second]);
    }
    for (int j = 0; j < lomeans.size(); ++j){
        if (lomeans[j] <= 0){
            bgprops2.push_back(0.0);
        }
        else{
            bgprops2.push_back(dirprops_bg[idx2idx[j]]/bgtot);
        }
    }

    double prob = 0.5;

    model_params.insert(make_pair(key_bg, dirprops_bg));
    
    vector<double> hi_nbinom_mu;
    vector<double> hi_nbinom_phi;
    for (int j = 0; j < lomeans.size(); ++j){
        double mu = -1;
        double phi = -1;
        if (idx2idx.count(j) > 0){
            mu = himeans[j] - lomeans[j];
            double var = mu*mu - (lomeans[j]*lomeans[j])/lophis[j];
            pair<double, double> mm = nbinom_moments(mu, var);
            mu = mm.first;
            phi = mm.second;
        }
        hi_nbinom_mu.push_back(mu);
        hi_nbinom_phi.push_back(phi);
        fprintf(stderr, "HI MU %s %f %f\n", labels[j].c_str(), mu, phi);
    }
    
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
    
    double llsum = 0.0;

    vector<vector<double> > new_low;
    vector<vector<double> > new_high;
    for (int j = 0; j < lomeans.size(); ++j){
        vector<double> v;
        new_low.push_back(v);
        new_high.push_back(v);
    }

    for (int i = 0; i < obs.size(); ++i){
        vector<pair<double, int> > lls;
        vector<string> names;
        vector<int> nums;

        vector<double> obsrowclean;

        // Sort in decreasing order of enrichment of p over expectation in background
        vector<pair<double, int> > dim_enrich_sort;
        
        double tot = obstots[i];
        
        bool all_bg = false;
        double bg_llr = 0.0;
        if (filt){
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
        if (filt && all_bg){
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
        
        double bgmu_max = 0.0;
        double bgphi_max = 0.0;

        double pbg_max = 0.0;
        double pbg_mean = 0.0;
        double pbg_var = 0.0;

        set<int> result;
        double nonbg = 0.0;
        
        bool print = false;
        if (bc2str(bcs[i]) == "ACCTGAAGTGCCTATA"){
        //if (bc2str(bcs[i]) == "AAACCCAAGAAATTCG"){
            print = true;
        }

        long int maxn = dim_enrich_sort.size();
        if (!sgrna){
            maxn = 2;
        }
        for (int k = 0; k < maxn; ++k){
            int j = dim_enrich_sort[k].second;
            nonbg += obs[i][j];

            included.members.insert(j);
            vector<double> params_model;
            pair<double, double> muphi_model;
            pair<double, double> slopeint_model;
            pair<double, double> varslopeint_model;
            vector<double> dircoeffs_model;
            vector<double> betacoeffs_model;

            if (model_params.count(included) == 0){
                // Generate parameters
                double totmu;
                double totphi;
                double slope;
                double intercept;
                double varslope;
                double varintercept;
                get_dirichlet(idx2idx,
                    included.members,
                    samps,
                    rowtots,
                    model_samps,
                    params_model,
                    totmu,
                    totphi,
                    slope,
                    intercept,
                    varslope,
                    varintercept,
                    pbg_mean,
                    pbg_var,
                    dircoeffs_model,
                    betacoeffs_model);

                // Store so we don't need to re-generate next time
                model_params.insert(make_pair(included, params_model));
                muphi_model = make_pair(totmu, totphi);
                model_muphi.insert(make_pair(included, muphi_model));
                slopeint_model = make_pair(slope, intercept);
                model_slopeint.insert(make_pair(included, slopeint_model));
                varslopeint_model = make_pair(varslope, varintercept);
                model_varslopeint.insert(make_pair(included, varslopeint_model));
                model_pbgmuvar.insert(make_pair(included, make_pair(pbg_mean, pbg_var)));
                model_dircoeffs.insert(make_pair(included, dircoeffs_model));
                model_betacoeffs.insert(make_pair(included, betacoeffs_model));
            }   
            else{
                params_model = model_params[included];
                muphi_model = model_muphi[included];
                slopeint_model = model_slopeint[included];
                varslopeint_model = model_varslopeint[included];
                pbg_mean = model_pbgmuvar[included].first;
                pbg_var = model_pbgmuvar[included].second;
                dircoeffs_model = model_dircoeffs[included];
                betacoeffs_model = model_betacoeffs[included];
            }
            
            vector<double> p;
            double totmin = exp(-slopeint_model.second/slopeint_model.first);
            double pbg = exp(slopeint_model.first*log(tot) + slopeint_model.second);
            double var = exp(varslopeint_model.first*log(tot) + varslopeint_model.second);
            double alpha_0 = (pbg*(1.0-pbg))/var - 1.0;
            double alpha_0_fallback = (pbg_mean*(1.0-pbg_mean))/pbg_var - 1.0;
            
            // TEST
            /*
            totmin = 0;
            pbg = expit(betacoeffs_model[0] * log(tot) + betacoeffs_model[1]);
            alpha_0 = exp(betacoeffs_model[2] * log(tot) + betacoeffs_model[3]);
            totmin = 0;
            */

            //set<int> test = result;
            //test.insert(j);
            //double pbg_inferred = get_cell_bg(obs[i], test, himeans, bgprops, himeans.size());
            //pbg_inferred /= tot;
            double alpha = pbg*alpha_0;
            double beta = (1.0-pbg)*alpha_0;
            double pbg_inferred = tot-nonbg;
            /*
            for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
                if (ix->first != j && result.find(ix->first) == result.end()){
                    pbg_inferred += obsrowclean[ix->second];
                } 
            }
            */
            if (pbg_inferred == 0){
                pbg_inferred = 1;
            }
            else if (pbg_inferred == tot){
                pbg_inferred = tot-1.0;
            }
            pbg_inferred /= tot;

            //pbg = pbg_inferred;
            if (false){
            //if (tot > totmin && alpha_0 > 0){
                //fprintf(stderr, "%s\t%f\t%f\t%f", bc2str(bcs[i]).c_str(), tot, pbg, alpha_0);
                
                /* 
                for (int i = 0; i < bgprops.size(); ++i){
                    p.push_back(pbg * bgprops[i]);
                }
                double sizetot = 0.0;
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    sizetot += hi_nbinom_mu[*m];
                }
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    double frac = hi_nbinom_mu[*m]/sizetot;
                    p[idx2idx[*m]] += (1.0-pbg)*frac;
                }
                for (int i = 0; i < p.size(); ++i){
                    params_model[i] = p[i] * alpha_0;
                } 
                */
                pbg = pbg_inferred;

                double t = 0.0;
                double pbgrm = 0.0;
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    //fprintf(stderr, " %s", labels[*m].c_str());
                    pbgrm += bgprops[idx2idx[*m]];
                }
                //fprintf(stderr, "\n");
                double psum = 0.0;
                for (int x = 0; x < params_model.size(); ++x){
                    p.push_back(pbg*(bgprops[x]/(1.0-pbgrm)));
                    t += params_model[x];
                }
                double sizetot = 0.0;
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    sizetot += himeans[*m];
                }
                for (set<int>::iterator m = included.members.begin(); m != included.members.end(); ++m){
                    double frac = himeans[*m]/sizetot;
                    p[idx2idx[*m]] = (1.0-pbg)*frac;
                }
                for (int x = 0; x < params_model.size(); ++x){
                    psum += p[x];
                    params_model[x] = p[x]*alpha_0;
                    //fprintf(stderr, "%s) %f | %f %f | %f\n", labels[models_kept[x]].c_str(), 
                    //    obsrowclean[x]/tot, p[x], params_model[x]/t, bgprops[x]);
                }
                //fprintf(stderr, "psum %f\n", psum);
                //fprintf(stderr, "\n");
                
                

            }
            else if (false){
                alpha = pbg_mean*alpha_0_fallback;
                beta = (1.0-pbg_mean)*alpha_0_fallback;
                alpha_0 = alpha_0_fallback;
                //fprintf(stderr, "B %f %f\n", tot, totmin);
                double t = 0;
                for (int x = 0; x < params_model.size(); ++x){
                    t += params_model[x];
                }
                
                for (int x = 0; x < params_model.size(); ++x){
                    p.push_back(params_model[x]/t);
                }
            }
            
            
            if (print){
                fprintf(stderr, "model");
                for (set<int>::iterator member = included.members.begin(); member != included.members.end(); ++member){
                    fprintf(stderr, " %s", labels[*member].c_str());
                }
                fprintf(stderr, "\n");
            }
           
            double partot = 0.0;   
            for (int x = 0; x < obsrowclean.size(); ++x){
                params_model[x] = exp(dircoeffs_model[x*2]*log(nonbg) + dircoeffs_model[x*2+1]);
                //params_model[x] = exp(dircoeffs_model[x*3]*log(nonbg)*log(nonbg) + dircoeffs_model[x*3+1]*log(nonbg) + dircoeffs_model[x*3+2]);
                partot += params_model[x];
            }
            for (int x = 0; x < params_model.size(); ++x){
                p.push_back(params_model[x]/partot);
                //p[x] = params_model[x] / partot;
                if (print){
                    fprintf(stderr, "%d %s)\t%f\n", models_kept[x], labels[models_kept[x]].c_str(), p[x]);
                }
            }

            
            
            /*
            if (bc2str(bcs[i]) == "GCCAGTGAGGTGATCG"){
                fprintf(stderr, "here\n");
                exit(0);
            }
            */

            //double ll = dmultinom(obsrowclean, p);
            //ll += ddirichlet(p, params_model);
            double ll = ddirichletmultinom(obsrowclean, params_model);
            if (print){
                fprintf(stderr, "partot %f\n", partot);
                fprintf(stderr, "LL %f\n", ll);
                fprintf(stderr, "LL alt %f\n", ddirichletmultinom(obsrowclean, params_model));
                fprintf(stderr, "\n");
            }
            //double ll = ddirichletmultinom(obsrowclean, params_model);
            //ll += dnbinom(tot, muphi_model.first, muphi_model.second);
            
            //ll += dbeta(pbg_inferred, alpha, beta);        
            
            if (isnan(ll)){
                fprintf(stderr, "???\n");
                for (int x = 0; x < obsrowclean.size(); ++x){
                    fprintf(stderr, "%d) %f %f\n", obsrowclean[x], p[x]);
                }
                fprintf(stderr, "%f\n", dmultinom(obsrowclean, p));
                fprintf(stderr, "%f %f %f\n", pbg_inferred, alpha, beta);
                fprintf(stderr, "%f\n", dbeta(pbg_inferred, alpha, beta));
                fprintf(stderr, "tot %f totmin %f\n", tot, totmin);
                fprintf(stderr, "fallback %f %f\n", alpha_0_fallback, pbg_mean);
                fprintf(stderr, "pbg %f\n", pbg);
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
                    
                    bgmu_max = pbg;
                    bgphi_max = alpha_0;
                    
                    pbg_max = pbg_inferred;
                }
            }
        }
        
        // Should we check for the case of 100% background? 
        if (add_bg){
            //double ll = ddirichletmultinom(obsrowclean, dirprops_bg);
            double ll = dmultinom(obsrowclean, bgprops);
            ll += dnbinom(tot, bgmuphi.first, bgmuphi.second);
            if (ll > llmax){
                if (result.size() == 1){
                    fprintf(stderr, "%s %f | %f\n", bc2str(bcs[i]).c_str(), ll, llmax);
                    fprintf(stderr, "%s\n", labels[*result.begin()].c_str());
                    for (int z = 0; z < obsrowclean.size(); ++z){
                        fprintf(stderr, "%s)\t%f\t%f\n", labels[models_kept[z]].c_str(), obsrowclean[z]/tot,
                            bgprops[z]);
                    }
                    fprintf(stderr, "\n");
                }
                llsecond = llmax;
                llmax = ll;
                result.clear();
            }
        }
        llsum += llmax;

        results.push_back(result);
        llrs.push_back(llmax-llsecond);
        
        double pbg = get_cell_bg(obs[i], result, himeans, bgprops2, himeans.size());
        if (pbg > 1){
            fprintf(stderr, "???\n");
            exit(1);
        }
        bgs.push_back(pbg);

        /* 
        if (false){
        //if (bgmu_max > 0 && bgphi_max > 0){
            double alpha = bgmu_max*bgphi_max;
            double beta = (1.0-bgmu_max)*bgphi_max;
            double beta_p = pbeta(pbg, alpha, beta);
            if (pbg > 0.5){
            //    beta_p = 1.0-beta_p;
            }
            if (false){
             //if (true){
            //if (beta_p > 0.9999){
                set<int> result_fallback;
                double llr_fallback;
                int n_assn = -1;
                if (!sgrna){
                    n_assn = 2;
                }
                assign_fallback(obs[i], lomeans, lophis, himeans, n_assn, result_fallback, llr_fallback);
                double pbg_fallback = get_cell_bg(obs[i], result_fallback, himeans, bgprops, himeans.size());
                pbg_fallback /= tot;
                fprintf(stdout, "%s", bc2str(bcs[i]).c_str());
                set<string> nsort1;
                set<string> nsort2;
                for (set<int>::iterator r = result.begin(); r != result.end(); ++r){
                    nsort1.insert(labels[*r]);
                }
                for (set<int>::iterator r = result_fallback.begin(); r != result_fallback.end(); ++r){
                    nsort2.insert(labels[*r]);
                }
                if (nsort1.size() == 2){
                    fprintf(stdout, "\t%s+%s", nsort1.begin()->c_str(), nsort1.rbegin()->c_str());
                }
                else{
                    fprintf(stdout, "\t%s", nsort1.begin()->c_str());
                }
                fprintf(stdout, "\t%f\t%f", pbg*tot, beta_p);
                if (nsort2.size() == 2){
                    fprintf(stdout, "\t%s+%s", nsort2.begin()->c_str(), nsort2.rbegin()->c_str());
                }
                else{
                    fprintf(stdout, "\t%s", nsort2.begin()->c_str());
                }
                fprintf(stdout, "\t%f\t%f\n", pbg_fallback*tot, pbeta(pbg_fallback, alpha, beta));

                //fprintf(stderr, "FALLBACK %f vs %f | %f vs %f\n", pbg_fallback, pbg, llr_fallback, llmax-llsecond);
                if (pbg_fallback < pbg){
                //if (pbg_fallback < pbg && llr_fallback > llmax-llsecond){
                    llrs[llrs.size()-1] = llr_fallback;
                    results[results.size()-1] = result_fallback;
                }  
            }
            //beta_p *= 2;
            //fprintf(stderr, "%f %f %f %f %f) %f\n", pbg, bgmu_max, bgphi_max, alpha, beta, beta_p);
            //fprintf(stdout, "%s\t%f\n", bc2str(bcs[i]).c_str(), beta_p); 
            //if (pbg > bgmu_max && beta_p < 0.001){
            //    fprintf(stdout, "%s\n", bc2str(bcs[i]).c_str());
            //}
        }
        */
        
        //bgs.push_back(pbg_max);
        //double totbg = 0.0;
        for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
            if (result.find(ix->first) == result.end()){
                new_low[ix->first].push_back(obsrowclean[ix->second]);
                //totbg += obsrowclean[ix->second];
            }
            else{
                new_high[ix->first].push_back(obsrowclean[ix->second]);
            }
        }
        //bgs.push_back(totbg/tot);
        //fprintf(stdout, "%s\t%f\t%f\n", bc2str(bcs[i]).c_str(), tot, totbg);
    }

    for (map<int, int>::iterator ix = idx2idx.begin(); ix != idx2idx.end(); ++ix){
        pair<double, double> muvar = welford(new_low[ix->first]);
        pair<double, double> muphi = nbinom_moments(muvar.first, muvar.second);
        pair<double, double> muvar_hi = welford(new_high[ix->first]);
        fprintf(stderr, "New Means %s: %f %f | %f %f\n", labels[ix->first].c_str(),
            muphi.first, muphi.second, muvar_hi.first, sqrt(muvar_hi.second));
        lomeans[muphi.first];
        lophis[muphi.second];
        himeans[muvar_hi.first];
    }
    return llsum;
}

double ll_beta_diff(double x, const map<string, double>& data_d, const map<string, int>& data_i){
    
    double a1 = data_d.at("a1");
    double b1 = data_d.at("b1");
    double a2 = data_d.at("a2");
    double b2 = data_d.at("b2");

    double ll1 = dbeta(x, a1, b1);
    double ll2 = dbeta(x, a2, b2);
    return ll1-ll2;
}

double median(vector<double>& nums){
    sort(nums.begin(), nums.end());
    if (nums.size() % 2 == 0){
        return (nums[nums.size()/2] + nums[nums.size()/2-1])*0.5;
    }
    else{
        return nums[(nums.size()-1) / 2];
    }
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
        
        fprintf(stderr, "Low percent background A %f B %f\n", ab_bgl.first, ab_bgl.second);
        fprintf(stderr, "High percent background A %f B %f\n", ab_bgh.first, ab_bgh.second);

        //mixtureDist lo("exponential", 1.0/muvar_bgl.first);
        mixtureDist lo("beta", ab_bgl.first, ab_bgl.second);
        mixtureDist hi("beta", ab_bgh.first, ab_bgh.second);

        //mixtureDist lo("beta", 0.25, 5.0);
        //mixtureDist hi("beta", 5.0, 0.25);
        //mixtureDist lo("beta", 1.0, 9.0);
        //mixtureDist hi("beta", 9.0, 1.0);
        vector<mixtureDist> dists;
        dists.push_back(lo);
        dists.push_back(hi);
        vector<double> w{ 0.5, 0.5};
        bgpercmod.init(dists, w);
        bgpercmod.fit(bgperc_all);
        bgpercmod.print();
        
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
        // Find intersection point between the two dists; use as cutoff.
        return intpt.solve(mean1, mean2);
    }
    else{
        fprintf(stderr, "Disabled filtering on ambient fraction; distribution not bimodal\n");
        return -1.0;
    }
}

/**
 * Model total background counts as two distributions: low bg count (true positive
 * assignments) and high bg count (false positive assignments).
 *
 * Returns cutoff value for background count.
 */
double sep_bg_dists(vector<double>& bgcount_all){
    if (!bimod_test(bgcount_all)){
        fprintf(stderr, "Disabled filtering step: background count distribution is unimodal\n");
        return -1;
    }
    else{
        
        // Get some initial value guesses
        pair<double, double> muvar_all = welford(bgcount_all);
        sort(bgcount_all.begin(), bgcount_all.end());
        double bglow = percentile(bgcount_all, 0.5);
        double bghigh = bgcount_all[bgcount_all.size()-1];
        
        double mid = (bglow + bghigh)/2.0;
        double sd = (mid-bglow)/2.0;
        double var = sd*sd;
        pair<double, double> nb1 = nbinom_moments(bglow, 2*bglow);
        pair<double, double> nb2 = nbinom_moments(bghigh, 2*bghigh);
        
        vector<mixtureDist> bgdists;
        mixtureDist dist_high("exponential", 1.0/bghigh);
        mixtureDist dist_low("negative_binomial", vector<double>{ bglow, 1.0 });
        bgdists.push_back(dist_low);
        bgdists.push_back(dist_high);
        
        vector<double> w{ 0.5, 0.5 }; 
        mixtureModel bgmod;
        bgmod.init(bgdists, w);
        vector<vector<double> > bgobs;
        for (int i = 0; i < bgcount_all.size(); ++i){
            bgobs.push_back(vector<double>{ bgcount_all[i] });
        }
        bgmod.fit(bgobs);
        
        fprintf(stderr, "===== Unfiltered vs excessive background counts per cell =====\n");
        fprintf(stderr, "Mean counts: %f\t%f\n", bgmod.dists[0].params[0][0], 
            1.0/bgmod.dists[1].params[0][0]);
        fprintf(stderr, "phi %f\n", bgmod.dists[0].params[0][1]);
        fprintf(stderr, "Dist weights: %f\t%f\n", bgmod.weights[0], bgmod.weights[1]);
        
        // Don't filter if we'll lose too much data
        // or if means are too close together.
        
        // Find point of intersection between the two distributions
        optimML::brent_solver intpoint_finder(diff_negbin_exp);
        intpoint_finder.set_root();
        intpoint_finder.add_data_fixed("mu", bgmod.dists[0].params[0][0]);
        intpoint_finder.add_data_fixed("phi", bgmod.dists[0].params[0][1]);
        intpoint_finder.add_data_fixed("lambda", bgmod.dists[1].params[0][0]);
        // Require that point of intersection be between the two means.
        double intpt = intpoint_finder.solve(bgmod.dists[0].params[0][0],
            1.0/bgmod.dists[1].params[0][0]);
        double intpt_p = pnbinom(intpt, bgmod.dists[0].params[0][0], bgmod.dists[0].params[0][1]);
         
        fprintf(stderr, "Intersection point %f, %f of lower dist\n", intpt, 1.0-intpt_p);
        
        // Find point of intersection between the two distributions (the point at which 
        // an observation becomes more likely under the filtered distribution)
        double lambda = bgmod.dists[1].params[0][1];
        double mu = bgmod.dists[0].params[0][0];
        double phi = bgmod.dists[0].params[0][1];

        if (intpt_p < 0.9){
            fprintf(stderr, "  disabled background filtering\n");
            return -1.0;
        }
        else{
            return intpt;
        }
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
    double prob){

    
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
    vector<double> model_means;
    vector<double> lomeans;
    vector<double> lophis;
    vector<double> himeans;
    fit_dists_2way(obscol, labels, lomeans, lophis, himeans, model_means, output_prefix);
    for (int i = 0; i < lomeans.size(); ++i){
        if (lomeans[i] > 1.0){
            lomeans[i] -= 1.0;
        }
    }
    fprintf(stderr, "Making preliminary assignments...\n");
    
    int nsamp = 1000;
    vector<set<int> > results_all;
    vector<double> llrs_all;
    vector<double> bg_all;

    vector<set<int> > results_prev;
    vector<double> llrs_prev;
    vector<double> bg_prev;

    double delta_thresh = 0.01;
    double llprev = 0.0;
    double delta = 999;
    while (delta > delta_thresh){

        results_prev = results_all;
        llrs_prev = llrs_all;
        bg_prev = bg_all;

        results_all.clear();
        llrs_all.clear();
        bg_all.clear();

        double loglik = assign_testnew(obs,
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
            filt || sgrna);
        break;
        if (llprev != 0){
           delta = loglik-llprev;
            if (delta < 0){
                results_all = results_prev;
                llrs_all = llrs_prev;
            }
        }
        fprintf(stderr, "%f -> %f delta %f\n", llprev, loglik, delta);
        llprev = loglik;
    }
    
    for (int i = 0; i < obs.size(); ++i){
        assn.emplace(bcs[i], results_all[i]);
        assn_llr.emplace(bcs[i], llrs_all[i]);
        cell_bg.emplace(bcs[i], bg_all[i]);
    }

    return;

    exit(1);



    // Now make preliminary assignments and learn background dists
    vector<double> bgmeans;
    double slope_mu;
    double intercept_mu;
    double slope_phi;
    double intercept_phi;
    assign_init(obs, labels, bgmeans,   
        lomeans, lophis, himeans, model_means, sgrna, slope_mu, intercept_mu,
        slope_phi, intercept_phi);
    
    // Refine total count -> perc bg relationship 
    vector<double> logtots;
    vector<double> logitpercs;
    vector<double> regtots;
    vector<double> regpercs;

    fprintf(stderr, "Making second round of assignments...\n");
    vector<double> bgcount_all;

    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr;
        double bgcount;
        assign_cell(bgmeans, model_means, -1, -1, obs[i],   
            labels.size(), slope_mu, intercept_mu, slope_phi, intercept_phi, 
            filt, assns, llr, bgcount, false, sgrna, -1);
        
        if (bgcount != -1){
            bgcount_all.push_back(bgcount);    
            if (bgcount > 0 && bgcount < tots[i] && tots[i] > 0){
                double perc = bgcount/tots[i];
                logtots.push_back(log(tots[i]));
                logitpercs.push_back(logit(perc));
                //bgperc_all_flat.push_back(perc);
                regtots.push_back(tots[i]);
                regpercs.push_back(perc);
            }
        }
    }
    // Update % bg from tot prediction
    //slope_int = linreg(logtots, logitpercs);
    fit_bg_predictors(logtots, logitpercs, regtots, regpercs, slope_mu,
        intercept_mu, slope_phi, intercept_phi);
    
    for (int i = 0; i < logtots.size(); ++i){
        double mu = slope_mu*logtots[i] + intercept_mu;
        double phi = slope_phi*logtots[i] + intercept_phi;

        fprintf(stdout, "%f\t%f\t%f\t%f\n", logtots[i], logitpercs[i], mu, exp(phi));
    }
    exit(0);
    
    pair<double, double> muvarbg = welford(bgcount_all);

    fprintf(stderr, "Making final assignments...\n");

    bgcount_all.clear();
    vector<double> bgperc_all;
    
    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr; 
        double bgcount;
        bool assigned = assign_cell(bgmeans, model_means, muvarbg.first, 
            muvarbg.second, obs[i], labels.size(), slope_mu, intercept_mu,
            slope_phi, intercept_phi, filt, assns, llr, bgcount, 
            false, sgrna, prob);
        
        bgcount_all.push_back(bgcount);
        double pbg = bgcount / tots[i];
        bgperc_all.push_back(pbg);
        
        assn.emplace(bcs[i], assns);
        assn_llr.emplace(bcs[i], llr);
        cell_bg.emplace(bcs[i], pbg);
    }

    if (sgrna){
        // Filter assignments
        double cutoff = sep_bg_dists_percent(bgperc_all);
        if (cutoff != -1){
            fprintf(stderr, "Percent background cutoff: %f\n", cutoff);
            for (int i = 0; i < obs.size(); ++i){
                if (bgperc_all[i] > cutoff){
                    // Set to WT
                    assn.erase(bcs[i]);
                    assn_llr.erase(bcs[i]);
                }
            }
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
    while((ch = getopt_long(argc, argv, "o:n:i:M:F:t:1:2:m:N:w:u:B:s:p:fCSUegch", 
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
    assign_ids(bc_tag_counts, assn, assn_llr, cell_bg, labels, output_prefix, sgrna, filt, prob);
    
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
