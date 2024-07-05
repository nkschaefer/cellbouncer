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
    fprintf(stderr, "   --filt -f Perform a filtering step to remove cells that do not fit the model well\n");
    fprintf(stderr, "       (these may correspond to high-order multiplets). Default: no filter\n");
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
    robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_ms_umis,
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

    for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x =
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

        if (m != -1){
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
        if (m != -1){
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

    vector<vector<double> > n;
    vector<vector<double> > bg;
    vector<vector<double> > m;

    mix.add_data_fixed("ndim", ndim);
    for (int k = 0; k < ndim; ++k){
        if (model_means[k] != -1){
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
    }

    char buf[50];
    for (int k = 0; k < ndim; ++k){
        if (model_means[k] != -1){
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
        return this_offtarget;
    } 
    else{
        double n_bg = p_inferred * obsrow[obsrow.size()-1];
        return n_bg;
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
    pair<double, double>& a_b,
    bool add_bg,
    set<int>& result,
    double& llr,
    double& bgcount,
    bool use_fixed_bg_count, 
    bool sgrna){
    
    // Get expected percent background counts from overall count
    double tot = obsrow[obsrow.size()-1];
    double pbg;

    if (!use_fixed_bg_count || bg_mean == -1 || tot <= bg_mean){
        pbg = a_b.first * log(tot) + a_b.second;
        pbg = exp(pbg)/(exp(pbg) + 1);
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
    
    double tot_fg = 0.0;
    
    long int maxn = dim_enrich_sort.size();
    if (!sgrna){
        maxn = 2;
    }
    for (int i = 0; i < maxn; ++i){
        int j = dim_enrich_sort[i].second;
        if (sizes[j] <= 0){
            continue;
        }

        included.insert(j);
        
        tot_counts += sizes[j];
        
        // Assume exponential variance (and use scaling factor to prevent overflow)
        var += pow(sizes[j], 2)/1e6;
        
        tot_fg += obsrow[j] * (1-bgprops[j]);

        vector<double> p;
        for (int k = 0; k < ndim; ++k){
            if (sizes[k] > 0){
                double p_this;
                if (included.find(k) != included.end()){
                    p_this = bgprops[k] * pbg + (1.0 - pbg) * (sizes[k]/tot_counts);
                }
                else{
                    p_this = bgprops[k] * pbg;
                }
                p.push_back(p_this);
            }
        }
        double ll = dmultinom(obsrowclean, p);
        
        if (bg_mean != -1 && bg_var != -1){
            double tot_exp = tot_counts + bg_mean;
            var += bg_var/1e6;
            
            // Use log-normal distribution
            double m2 = pow(tot_exp, 2)/1e6;
            double sigma = sqrt(log(var/m2 + 1.0));
            double mu = log(tot_exp) - log(var/m2 + 1.0)/2.0;
            ll += dnorm(log(obsrow[obsrow.size()-1]), mu, sigma);
        }
        lls.push_back(make_pair(-ll, i));
    }
    
    // Should we check for the case of 100% background? 
    if (add_bg){
        vector<double> params_bg;
        for (int j = 0; j < ndim; ++j){
            if (sizes[j] > 0){
                params_bg.push_back(pbg * bgprops[j]);
            }
        }
        double ll = dmultinom(obsrowclean, params_bg);
        
        if (bg_mean != -1 && bg_var != -1){
            double cbg = obsrow[obsrow.size()-1];

            // Use log-normal distribution
            double m2 = bg_mean*bg_mean;
            double sigma = sqrt(log(bg_var/m2 + 1.0));
            double mu = log(bg_mean) - log(bg_var/m2 + 1.0)/2.0;

            ll += dnorm(log(cbg), mu, sigma);
        }
        lls.push_back(make_pair(-ll, -1));
    }
    sort(lls.begin(), lls.end());
    llr = -lls[0].first - -lls[1].first;
    int chosen = lls[0].second;
    if (chosen == -1){
        bgcount = obsrow[obsrow.size()-1];
        return false;    
    }
    else{
        for (int i = 0; i <= chosen; ++i){
            result.insert(dim_enrich_sort[i].second);
        }
        bgcount = get_cell_bg(obsrow, result, sizes, bgprops, ndim);
        return true;
    }
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
            double highval = obscol[i][obscol[i].size()-1];
            vector<mixtureDist> dists;
            mixtureDist low("negative_binomial", vector<double>{ lowval, 1.0 });
            mixtureDist high("exponential", 1.0/highval);
            dists.push_back(low);
            dists.push_back(high);
            mixtureModel mod(dists);
            mod.fit(obscol[i]);
            double dist_sep = pnbinom(1.0/mod.dists[1].params[0][0],
                mod.dists[0].params[0][0], mod.dists[0].params[0][1]);
            fprintf(stderr, "%s:\t%f\t%f\n", labels[i].c_str(),
                mod.dists[0].params[0][0], 1.0/mod.dists[1].params[0][0]);
            
            fprintf(outf, "%s\t%f\t%f\t%f\t%f\n", labels[i].c_str(),
                mod.weights[0], mod.dists[0].params[0][0], mod.dists[0].params[0][1],
                mod.dists[1].params[0][0]);

            if (mod.weights[1] >= mod.weights[0] || dist_sep < 0.99){
                lomeans.push_back(0);
                lophis.push_back(0);
                himeans.push_back(0);
                hiweights.push_back(0);
                loweights.push_back(1);
            }
            else{
                lomeans.push_back(mod.dists[0].params[0][0]);
                lophis.push_back(mod.dists[0].params[0][1]);
                himeans.push_back(1.0/mod.dists[1].params[0][0]);
                hiweights.push_back(mod.weights[1]);
                loweights.push_back(mod.weights[0]);
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
    pair<double, double> hiweightmv = welford(hiweights);
    pair<double, double> loweightmv = welford(loweights);
    pair<double, double> hiweightbeta = beta_moments(hiweightmv.first, hiweightmv.second);
    pair<double, double> loweightbeta = beta_moments(loweightmv.first, loweightmv.second);
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
        else{
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
                double dbeta1 = dbeta(hiweights[i], loweightbeta.first, loweightbeta.second);
                double dbeta2 = dbeta(hiweights[i], hiweightbeta.first, hiweightbeta.second);
                if (dbeta2 - dbeta1 < 0){
                    fprintf(stderr, "  Remove label %s, LLR(weight) = %f\n", labels[i].c_str(),
                        dbeta1-dbeta2);
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
}

/**
 * Do initial assignments to learn the % of background/ambient tags
 * composed of each label, and to find the relationship between
 * log total tag count and logit percent background (returned as
 * linear regression parameters)
 */
pair<double, double> assign_init(vector<vector<double> >& obs,
    vector<string>& labels,
    vector<double>& bgmeans,
    vector<double>& lomeans,
    vector<double>& lophis,
    vector<double>& himeans,
    vector<double>& model_means,
    bool sgrna){

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
    vector<unsigned long> logtots_bc;
    
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

            if (p_high > p_thresh){
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
        }
    }

    fprintf(stderr, "===== Ambient tag profile: =====\n");
    // Learn the proportions of each label in ambient/background counts
    for (int i = 0; i < labels.size(); ++i){
        bgmeans[i] /= bgcount;
        fprintf(stderr, "%s:\t%f\n", labels[i].c_str(), bgmeans[i]);
    }
    
    // Determine the relationship between total counts and percent background.
    // Model this as a linear fit between log transformed count and logit-transformed
    // percent background.
    pair<double, double> slope_int = linreg(logtots, logitpercs);
    
    return slope_int;
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

double ll_beta_diff(double x, const map<string, double>& data_d, const map<string, int>& data_i){
    
    double a1 = data_d.at("a1");
    double b1 = data_d.at("b1");
    double a2 = data_d.at("a2");
    double b2 = data_d.at("b2");

    double ll1 = dbeta(x, a1, b1);
    double ll2 = dbeta(x, a2, b2);
    return ll1-ll2;
}

/**
 * Main function for assigning cell identities from data.
 */
void assign_ids(robin_hood::unordered_map<unsigned long, map<int, long int> >& bc_tag_counts,
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    robin_hood::unordered_map<unsigned long, double>& cell_bg,
    vector<string>& labels,
    bool filt,
    string& output_prefix, 
    bool sgrna){
    
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
        }
        double tot = 0.0;
        for (map<int, long int>::iterator y = x->second.begin(); y != x->second.end(); 
            ++y){
            row[y->first] += (double)y->second;
            tot += (double)y->second;
            double val = (double)y->second + 1;
            obscol[y->first].push_back(val);
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
    fprintf(stderr, "Making preliminary assignments...\n");

    // Now make preliminary assignments and learn background dists
    vector<double> bgmeans;
    pair<double, double> slope_int = assign_init(obs, labels, bgmeans,   
        lomeans, lophis, himeans, model_means, sgrna);
    
    // Make second-round assignments and store the maximum likelihood-inferred
    // background counts    
    vector<double> bgcount_all; 
    vector<vector<double> > bgperc_all;
    vector<double> bgperc_all_flat;

    fprintf(stderr, "Making second round of assignments to infer ambient counts per cell...\n");
    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr;
        double bgcount; 
        assign_cell(bgmeans, model_means, -1, -1, obs[i],   
            labels.size(), slope_int, true, assns, llr, bgcount, false, sgrna);
        if (bgcount != -1){
            bgcount_all.push_back(bgcount);    
            if (bgcount > 0 && bgcount < tots[i]){
                bgperc_all.push_back(vector<double>{ bgcount / tots[i] }); 
                bgperc_all_flat.push_back(bgcount/tots[i]);
            }
        }
    }
    
    // Use what we learned from round2 of assignments to learn the mean & variance
    // of background counts per cell 
    pair<double, double> muvarbg = welford(bgcount_all);
    
   
    mixtureModel bgmod;
    mixtureModel bgpercmod; 
    double intptx = -1;

    if (sgrna){
        
        // Do not filter the normal way. Filter on percent background counts instead.
        
        if (bimod_test(bgperc_all_flat)){
            filt = false;
            mixtureDist lo("beta", 1.0, 9.0);
            mixtureDist hi("beta", 9.0, 1.0);
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
            double mean1 = bgpercmod.dists[0].params[0][0] / (bgpercmod.dists[0].params[0][0] + bgpercmod.dists[0].params[0][1]);
            double mean2 = bgpercmod.dists[1].params[0][0] / (bgpercmod.dists[1].params[0][0] + bgpercmod.dists[0].params[0][1]);
            intptx = intpt.solve(mean1, mean2);
            fprintf(stderr, "intpt %f\n", intptx);
        }
        else{
            fprintf(stderr, "Disabled filtering on ambient fraction; distribution not bimodal\n");
        }
    }    
    if (filt){ 
        if (!bimod_test(bgcount_all)){
            fprintf(stderr, "Disabled filtering step: background count distribution is unimodal\n");
            filt = false;
        }
    }
    if (filt){ 
        
        fprintf(stderr, "===== Unfiltered vs excessive background counts per cell =====\n");

        pair<double, double> muvar_all = welford(bgcount_all);
        sort(bgcount_all.begin(), bgcount_all.end());
        double bglow = percentile(bgcount_all, 0.1);
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
        bgmod.init(bgdists, w);
        vector<vector<double> > bgobs;
        for (int i = 0; i < bgcount_all.size(); ++i){
            bgobs.push_back(vector<double>{ bgcount_all[i] });
        }
        bgmod.fit(bgobs);
        fprintf(stderr, "Mean counts: %f\t%f\n", bgmod.dists[0].params[0][0], 
            1.0/bgmod.dists[1].params[0][0]);
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
            filt = false;
        }
    }

    int count_filt = 0;
    
    fprintf(stderr, "Making final assignments...\n");
    
    pair<double, double> muphi = nbinom_moments(muvarbg.first, muvarbg.second);
    double p_bg_thresh = 0.01;

    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr; 
        double bgcount;
        bool assigned = assign_cell(bgmeans, model_means, muvarbg.first, 
            muvarbg.second, obs[i], labels.size(), slope_int, true, assns, llr, bgcount, 
            false, sgrna);
        
        double pbg = bgcount / tots[i];
        cell_bg.emplace(bcs[i], pbg);
        
        if (assigned && sgrna){
            if (intptx > 0 && pbg > intptx){
                assns.clear();
                llr = 0.0;
            }
            /*
            if (pbg == 1.0){
                assns.clear();
            }
            else{
                vector<double> row{ pbg };
                double ll_lo = bgpercmod.dists[0].loglik(row);
                double ll_hi = bgpercmod.dists[1].loglik(row);
                if (ll_hi > ll_lo){
                    llr = ll_hi - ll_lo;
                    assns.clear();
                }
            }
            */
        }

        // Because exponential dists can catch low-count stuff too, make sure it's more likely
        // under the high dist AND that its count is high, not very low.
        if (assigned && filt && bgcount > bgmod.dists[0].params[0][0]){
            
            vector<double> row{ bgcount };
            double ll1 = bgmod.dists[0].loglik(row);
            double ll2 = bgmod.dists[1].loglik(row);
            double prob_high;
            if (ll1 < ll2){
                double denom = 1.0 + pow(2, ll2-ll1);
                prob_high = pow(2, ll2-ll1) / denom;
            }
            else{
                double denom = 1.0 + pow(2, ll1-ll2);
                prob_high = 1.0 - pow(2, ll1-ll2) / denom;
            }
            
            if (ll2 > ll1){
                double denom = 1.0 + pow(2, ll2-ll1);
                double prob_high = pow(2, ll2-ll1) / denom;
                if (prob_high > 0.5){
                    assigned = false;
                }
            }
        }
        if (assigned){
            assn.emplace(bcs[i], assns);
            assn_llr.emplace(bcs[i], llr);
        }
        else{
            count_filt++;
        }
    }
    
    if (filt){ 
        fprintf(stderr, "Filtered %d of %ld assignments\n", count_filt, obs.size());
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
void check_missing_ms_wells(robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_ms_umis,
    vector<string>& ms_wells,
    map<string, string>& well2name,
    string& outfilename){
    
    // Count each well across all cells
    vector<pair<int, string> > counts;
    for (int i = 0; i < ms_wells.size(); ++i){
        counts.push_back(make_pair(0, ms_wells[i]));
    }
    for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x = 
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
    robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_tag_umis){
    
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
                    // Cell hashing barcodes are assumed to occur at teh beginning of R2, forward orientation
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
                        vector<umi_set*> v(seqlist.size(), NULL);
                        bc_tag_umis.emplace(scanner.barcode, v);
                    }
                    umi_set* ptr = bc_tag_umis[scanner.barcode][seq2idx[seqmatch]];
                    if (ptr == NULL){
                        // Initialize.
                        bc_tag_umis[scanner.barcode][seq2idx[seqmatch]] = new umi_set(scanner.umi_len);
                        ptr = bc_tag_umis[scanner.barcode][seq2idx[seqmatch]];
                    }
                    ptr->add(this_umi); 
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {    
    
    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"seqs", required_argument, 0, 's'},
       {"names", required_argument, 0, 'N'},
       {"whitelist", required_argument, 0, 'w'},
       {"cell_barcodes", required_argument, 0, 'B'},
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
    bool filt = false;
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
    while((ch = getopt_long(argc, argv, "o:n:i:M:F:t:1:2:m:N:w:u:B:s:CSUegfch", 
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
            case 'f':
                filt = true;
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
            parse_mex(cell_barcodesfn, input_features, input_mtx, 
                bc_tag_counts, labels, featuretype); 
            
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
            robin_hood::unordered_map<unsigned long, vector<umi_set*> > bc_tag_umis;
            
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
            for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x = 
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
    assign_ids(bc_tag_counts, assn, assn_llr, cell_bg, labels, filt, output_prefix, sgrna);
    
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
