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
    fprintf(stderr, "\n   ===== REQUIRED UNLESS LOADING COUNTS FROM A PREVIOUS RUN =====\n");
    fprintf(stderr, "   --read1 -1 Forward read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --read2 -2 Reverse read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --whitelist -w Cell barcode whitelist file (i.e. included with 10X Genomics\n");
    fprintf(stderr, "       Cellranger software in cellranger-x.y.z/lib/python/cellranger/barcodes/\n");
    fprintf(stderr, "\n   ===== STRONGLY RECOMMENDED =====\n");
    fprintf(stderr, "   --mapping -m 2-column tab separated file mapping string interpretation of\n");
    fprintf(stderr, "       MULTIseq label (i.e. unique identifier, or treatment name) to MULTIseq\n");
    fprintf(stderr, "       well identifier.\n");
    fprintf(stderr, "   --cell_barcodes -B The filtered list of valid cell barcodes, i.e. from cellranger or\n");
    fprintf(stderr, "       STARsolo. Only cells in this list will be processed.\n");
    fprintf(stderr, "\n   ===== OPTIONAL =====\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    fprintf(stderr, "   --sgRNA -g Specifies that data is from sgRNA capture. This affects how sequence\n");
    fprintf(stderr, "       matching in reads is done and slightly changes the output files. Also makes\n");
    fprintf(stderr, "       the default separator for multiple IDs a comma instead of +.\n");
    fprintf(stderr, "   --barcodes -b A path to a file listing MULTIseq well IDs and corresponding\n");
    fprintf(stderr, "       barcodes, tab separated. If you do not provide one, the default file (in\n");
    fprintf(stderr, "       the data directory) will be used.\n");
    fprintf(stderr, "   --mismatches -M The number of allowable mismatches to a MULTIseq barcode to accept\n");
    fprintf(stderr, "       it. Default = -1 (take best match overall, if there are no ties)\n");
    fprintf(stderr, "   --filt -f Perform a filtering step to remove cells that do not fit the model well\n");
    fprintf(stderr, "       (these may correspond to high-order multiplets). Default: no filter\n");
    fprintf(stderr, "   --comma -c By default, cells assigned multiple identities will have these\n");
    fprintf(stderr, "       identities separated by + in the .assignments and .ids files. If\n");
    fprintf(stderr, "       identities contain + within them, though, this option switches to a\n");
    fprintf(stderr, "       comma separator.\n");
    fprintf(stderr, "   --batch_id -i Optional string to append to cell barcodes (with - separator).\n");
    fprintf(stderr, "       Default = no batch ID. Use this if you plan to combine multiple data sets\n");
    fprintf(stderr, "       together for analysis. You should process each separately here, but the\n");
    fprintf(stderr, "       batch_id appended here will need to match those you use in your analysis\n");
    fprintf(stderr, "       program (i.e. scanpy).\n");
    exit(code);
}

/**
 * Returns barcode length.
 */
int load_well2bc(string& filename, 
    map<unsigned long, string>& bc2well, 
    vector<string>& bclist){
    
    ifstream inf(filename);
    string bc_str;
    string well;
    bc bc_bitset;
    int bc_len = -1;
    while (inf >> well >> bc_str){
        if (bc_len == -1){
            bc_len = bc_str.length();
        }
        else if (bc_len != bc_str.length()){
            fprintf(stderr, "ERROR: mismatching MULTIseq barcode lengths in file %s:\n",
                filename.c_str());
            fprintf(stderr, "%d vs %ld\n", bc_len, bc_str.length());
            exit(1);
        }
        str2bc(bc_str.c_str(), bc_bitset, bc_len);
        bc2well.insert(make_pair(bc_bitset.to_ulong(), well));
        bclist.push_back(bc_str);
    }
    return bc_len;
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
 */
void dump_counts(string& filename, 
    robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_ms_umis,
    vector<string>& ms_names,
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts){

    FILE* f = fopen(filename.c_str(), "w");
    
    // Print header line
    fprintf(f, "bc");
    for (vector<string>::iterator n = ms_names.begin(); n != ms_names.end(); ++n){
        if (n->length() > 0){
            fprintf(f, "\t%s", n->c_str());
        }
    }
    fprintf(f, "\n");

    for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x =
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){

        vector<pair<int, string> > v;
        bc_ms_counts.emplace(x->first, v);

        string bc_str = bc2str(x->first);
        fprintf(f, "%s", bc_str.c_str());
        for (int i = 0; i < x->second.size(); ++i){
            if (ms_names[i].length() > 0){
                // This MULTIseq barcode was intended to be present in the data set
                if (x->second[i] == NULL){
                    fprintf(f, "\t0");
                    bc_ms_counts[x->first].push_back(make_pair(0,
                        ms_names[i]));
                }
                else{
                    bc_ms_counts[x->first].push_back(make_pair(-x->second[i]->count(),
                        ms_names[i]));
                    fprintf(f, "\t%d", x->second[i]->count());
                }
            }
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
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    vector<string>& labels){
    
    string line;
    ifstream inf(filename);
    bool first = true;
    bc cur_bc;

    vector<string> header;

    while (getline(inf, line)){
        istringstream splitter(line);
        string field;
        int idx = 0;
        while(getline(splitter, field, '\t')){
            if (first){
                header.push_back(field);
            }
            else{
                if (idx == 0){
                    str2bc(field.c_str(), cur_bc);
                    vector<pair<int, string> > v;
                    bc_ms_counts.emplace(cur_bc.to_ulong(), v);      
                }
                else{
                    int val = atoi(field.c_str());
                    if (val > 0){
                        bc_ms_counts[cur_bc.to_ulong()].push_back(make_pair(-val, header[idx]));
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

map<int, int> ddim;
map<int, pair<int, int> > ddoub;

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
 * For use with Brent 1-dim solver: infer proportion background
 * in a cell
 */
double ll_mix(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx1 = data_i.at("idx1");
    int idx2 = -1;
    double mixprop;
    if (data_i.count("idx2") > 0){
        idx2 = data_i.at("idx2");
        mixprop = data_d.at("mixprop");
    }
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
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgfrac = data_d.at(bufstr);
        if (idx2 != -1){
            if (i == idx1){
                p.push_back(param*bgfrac + (1.0 - param)*mixprop);
                x.push_back(x_i);
            }
            else if (i == idx2){
                p.push_back(param*bgfrac + (1.0-param)*(1.0 - mixprop));
                x.push_back(x_i);
            }
            else{
                p.push_back(param*bgfrac);
                x.push_back(x_i); 
            }
        }
        else{
            if (i == idx1){
                p.push_back((1.0 - param) + param*bgfrac);
                x.push_back(x_i);
            }
            else{
                p.push_back(param*bgfrac);
                x.push_back(x_i);
            }
        }
    }

    return dmultinom(x, p)/log2(exp(1));
}

/**
 * For use with Brent's 1-D solver: derivative function to 
 * find the MLE percent background in a cell
 */
double dll_mix(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx1 = data_i.at("idx1");
    int idx2 = -1;
    double mixprop;
    if (data_i.count("idx2") > 0){
        idx2 = data_i.at("idx2");
        mixprop = data_d.at("mixprop");
    }
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
        if (idx2 != -1){
            if (i == idx1){
                dL_df += (x_i*(bgfrac - mixprop))/((1-param)*mixprop + param*bgfrac);
            }
            else if (i == idx2){
                dL_df += (x_i*(mixprop + bgfrac - 1.0))/((1-param)*(1-mixprop) + param*bgfrac);
            }
            else{
                dL_df += x_i / param;
            }
        }
        else{
            if (i == idx1){
                dL_df += ((bgfrac-1.0)*x_i)/(param*(bgfrac-1.0) + 1.0);
            }
            else{
                dL_df += x_i/param;
            }
        }
    }
    
    return dL_df;
}

/**
 * Function for inferring background proportions of tag counts as
 * part of EM fitting callback function
 */
double bg_props(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    
    return params[idx];
}

/**
 * Derivative of function for inferring background proportions of tag counts
 * as part of EM fitting callback function
 */
void dbg_props(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    results[idx] += 1.0;
}

/**
 * Function for inferring singlet proportions as part of EM fitting callback
 * function
 */
double singlet_props(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int ndim = data_i.at("ndim");
    int idx = data_i.at("idx");
    int target = data_i.at("target");

    double frac_bg = params[ndim+target];
    double p_bg = params[idx];

    if (idx == target){
        return 1.0 - frac_bg + frac_bg*p_bg;        
    }
    else{
        return frac_bg*p_bg;
    }
}

/**
 * Derivative function for inferring singlet proportions as part of
 * EM fitting callback function
 */
void dsinglet_props(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int ndim = data_i.at("ndim");
    int idx = data_i.at("idx");
    int target = data_i.at("target");

    double frac_bg = params[ndim+target];
    double p_bg = params[idx];

    if (idx == target){
        results[idx] += frac_bg;
        results[ndim+target] += (p_bg - 1.0);
    }
    else{
        results[idx] += frac_bg;
        results[ndim+target] += p_bg;
    }
}

/**
 * Function for inferring doublet proportions in EM fitting
 * callback function (through optimization)
 */
double doublet_props(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int ndim = data_i.at("ndim");
    int idx = data_i.at("idx");
    int target = data_i.at("target");
    int target2 = data_i.at("target2");

    double frac_bg1 = params[ndim+target];
    double frac_bg2 = params[ndim+target2];
    double p_bg = params[idx];
    
    double reads1 = params[ndim*2 + target];
    double reads2 = params[ndim*2 + target2];
    
    double frac1 = reads1/(reads1+reads2);
    double frac2 = reads2/(reads1+reads2);

    if (idx == target){
        return frac1*(1.0 - frac_bg1 + frac_bg1*p_bg) + frac2*(frac_bg2*p_bg);    
    }
    else if (idx == target2){
        return frac1*(frac_bg1*p_bg) + frac2*(1.0 - frac_bg2 + frac_bg2*p_bg);
    }
    else{
        return frac1*frac_bg1*p_bg + frac2*frac_bg2*p_bg;
    }
}

/**
 * Derivative function for inferring doublet proportions in EM callback
 * function (through optimization)
 */
void ddoublet_props(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int ndim = data_i.at("ndim");
    int idx = data_i.at("idx");
    int target = data_i.at("target");
    int target2 = data_i.at("target2");

    double frac_bg1 = params[ndim+target];
    double frac_bg2 = params[ndim+target2];
    double p_bg = params[idx];
    
    double reads1 = params[ndim*2 + target];
    double reads2 = params[ndim*2 + target2];

    double frac1 = reads1/(reads1+reads2);
    double frac2 = reads2/(reads1+reads2);

    if (idx == target){
        results[ndim+target] += ((p_bg-1.0)*reads1)/(reads1+reads2);
        results[ndim+target2] += (p_bg*reads2)/(reads1+reads2);
        results[idx] += (frac_bg1*reads1 + frac_bg2*reads2)/(reads1+reads2);
        results[ndim*2 + target] += (reads2*(frac_bg1*(p_bg-1.0) - frac_bg2*p_bg + 1.0))/pow(reads1+reads2,2);
        results[ndim*2 + target2] += (reads1*(frac_bg1*-p_bg + frac_bg1 + frac_bg2*p_bg-1.0))/pow(reads1+reads2,2);
    }
    else if (idx == target2){
        results[ndim+target] += (p_bg*reads1)/(reads1+reads2);
        results[ndim+target2] += ((p_bg-1.0)*reads2)/(reads1+reads2);
        results[idx] += (frac_bg1*reads1 + frac_bg2*reads2)/(reads1+reads2);
        results[ndim*2 + target] += (reads2*(frac_bg1*p_bg + frac_bg2*-p_bg + frac_bg2 - 1.0))/pow(reads1+reads2, 2);
        results[ndim*2 + target2] += (reads1*(-frac_bg1*p_bg + frac_bg2*(p_bg-1.0) + 1.0))/pow(reads1+reads2, 2);
    }
    else{
        results[ndim+target] += (p_bg*reads1)/(reads1+reads2);
        results[ndim+target2] += (p_bg*reads2)/(reads1+reads2);
        results[idx] += (frac_bg1*reads1 + frac_bg2*reads2)/(reads1+reads2);
        results[ndim*2 + target] += (p_bg*reads2*(frac_bg1-frac_bg2))/pow(reads1+reads2,2);
        results[ndim*2 + target2] += (p_bg*reads1*(frac_bg2-frac_bg1))/pow(reads1+reads2,2);
    }
}

/**
 * After each round of EM fitting (M-step), reconciles all shared parameters: 
 * each distribution in the mixture model depends on proportion of counts of
 * each tag in ambient RNA (background distribution), an average fraction of 
 * foreground RNA per cell of each identity, and mean total tag counts per cell given
 * each possible identity.
 *
 * Solves for these parameters given current model parameters, and then updates model
 * parameters to reflect current solutions to these shared parameters.
 */
void callback_func(mixtureModel& mod, vector<double>& params){
    
    vector<int> idx1;
    vector<int> idx2;
    vector<double> frac;
     
    int ndim = mod.dists[0].params[0].size();
    double meansize = 0.0;
    int countsize = 0;
    for (int i = 0; i < ndim; ++i){
        params[ndim*2+i] = mod.dists[i].params[1][0];
        meansize += mod.dists[i].params[1][0];
        countsize++;
    }
    meansize /= (double)countsize;
    optimML::multivar_sys_solver solver(vector<double>{ } );
    
    // Add in proportions in background
    vector<double> pg;
    for (int i = 0; i < ndim; ++i){
        pg.push_back(params[i]);
    }
    solver.add_param_grp(pg);
    
    // Add in percent composed of background
    double meanpbg = 0.0;
    for (int i = ndim; i < ndim*2; ++i){
        meanpbg += params[i];
    }
    meanpbg /= (double)ndim;

    for (int i = ndim; i < ndim*2; ++i){
        solver.add_one_param(params[i]);
        solver.constrain_01(i);
    }
    
    // Add in read counts
    // Ensure all counts stay positive
    for (int i = ndim*2; i < ndim*3; ++i){
        solver.add_one_param(params[i]);
        solver.constrain_pos(i);
    }
    
    solver.add_data_fixed("ndim", ndim);

    for (map<int, pair<int, int> >::iterator d = ddoub.begin(); d != ddoub.end(); ++d){
        if (mod.weights[d->second.first] == 0.0 || mod.weights[d->second.second] == 0.0){
            mod.weights[d->first] = 0.0;
        }
    }
    double ws = 0.0;
    for (int j = 0; j < mod.n_components; ++j){
        ws += mod.weights[j];    
    }
    
    for (int j = 0; j < mod.n_components; ++j){
        mod.weights[j] /= ws;
    }
     
    for (int i = 0; i < mod.dists.size(); ++i){
        if (ddoub.count(i) > 0){
            int dim1 = ddim[ddoub[i].first];
            int dim2 = ddim[ddoub[i].second];

            double prop1 = mod.dists[i].params[0][dim1];
            double prop2 = mod.dists[i].params[0][dim2];
            
            if (mod.weights[i] > 0.0){
                for (int j = 0; j < ndim; ++j){
                    solver.add_equation(doublet_props, ddoublet_props, 
                        mod.dists[i].params[0][j], mod.weights[i]);
                    solver.add_data("idx", j);
                    solver.add_data("target", dim1);
                    solver.add_data("target2", dim2);
                }     
            }
        }
        else if (ddim.count(i) > 0 && mod.weights[i] > 0.0){
            int dim = ddim[i];
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                solver.add_equation(singlet_props, dsinglet_props, 
                    mod.dists[i].params[0][j], mod.weights[i]);
                solver.add_data("idx", j);
                solver.add_data("target", dim);
                solver.add_data("target2", -1);
            }
        }
        else if (ddim.count(i) == 0 && mod.weights[i] > 0.0){
            // Background dist.
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                solver.add_equation(bg_props, dbg_props, 
                    mod.dists[i].params[0][j], mod.weights[i]);
                solver.add_data("idx", j);
                solver.add_data("target", -1);
                solver.add_data("target2", -1);
            }
        }
    }
    solver.solve();
   
    /*
     * Uncomment to print current values of shared parameters on each iteration
     *
    for (int i = 0; i < solver.results.size(); ++i){
        if (i < ndim){
            fprintf(stderr, "BGPROP %d) %f\n", i, solver.results[i]);
        }
        else if (i >= ndim && i < ndim*2){
            fprintf(stderr, "%% BG %d) %f\n", i-ndim, solver.results[i]);
        }
        else{
            fprintf(stderr, "NREADS %d) %f\n", i-ndim*2, solver.results[i]);
        }
    }
    fprintf(stderr, "\n");
    */
    
    // Copy results back into shared_params
    for (int i = 0; i < solver.results.size(); ++i){
        params[i] = solver.results[i];
    }

    // Update model distribution parameters based on these current estimates of shared
    // parameters.

    for (map<int, pair<int, int> >::iterator dd = ddoub.begin(); dd != ddoub.end(); ++dd){
        int idx1 = ddim[dd->second.first];
        int idx2 = ddim[dd->second.second];
        
        double fracbg1 = solver.results[ndim+idx1];
        double fracbg2 = solver.results[ndim+idx2];

        double mixprop1 = solver.results[ndim*2+idx1]/
            (solver.results[ndim*2+idx1] + solver.results[ndim*2+idx2]);
        double mixprop2 = solver.results[ndim*2+idx2]/
            (solver.results[ndim*2+idx1] + solver.results[ndim*2+idx2]);

        for (int i = 0; i < ndim; ++i){
            if (i == idx1){
                mod.dists[dd->first].params[0][i] = 
                    mixprop1*(1.0 - fracbg1 + fracbg1*solver.results[i]) + 
                    mixprop2*(fracbg2*solver.results[i]);
            }
            else if (i == idx2){
                mod.dists[dd->first].params[0][i] = mixprop1*(fracbg1*solver.results[i]) + 
                    mixprop2*(1.0 - fracbg2 + fracbg2*solver.results[i]);
            }
            else{
                mod.dists[dd->first].params[0][i] = mixprop1*fracbg1*solver.results[i] + 
                    mixprop2*fracbg2*solver.results[i];
            }
        }
        
        mod.dists[dd->first].params[1][0] = solver.results[ndim*2+idx1] + 
            solver.results[ndim*2+idx2];
    }
    
    for (map<int, int>::iterator dd = ddim.begin(); dd != ddim.end(); ++dd){
        
        for (int i = 0; i < ndim; ++i){

            double fracbg = solver.results[ndim+dd->second];
            
            if (i == dd->second){
                mod.dists[dd->first].params[0][i] = 1.0 - fracbg + fracbg*solver.results[i];
            }
            else{
                mod.dists[dd->first].params[0][i] = fracbg*solver.results[i];
            }
        }
        
        mod.dists[dd->first].params[1][0] = solver.results[ndim*2+dd->second];
    }
    
    if (ddim.count(mod.dists.size()-1) == 0 && ddoub.count(mod.dists.size()-1) == 0){
        int bgidx = mod.dists.size()-1;
        if (mod.dists[bgidx].name == "BG"){
            for (int i = 0; i < ndim; ++i){
                mod.dists[bgidx].params[0][i] = solver.results[i];
            }
        }
    }
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
    double intercept = (ysum * xsum_sq - xsum * xysum) / (num * xsum_sq - xsum * xsum);
    return make_pair(slope, intercept);

}

double ll_mix_all(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int ndim = data_i.at("ndim");
    vector<double> x;
    vector<double> p;
    static char buf[50];
    static string bufstr;
    double fbg = params[params.size()-1];
    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        x.push_back(data_d.at(bufstr));
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgp = data_d.at(bufstr);
        double this_p = fbg * bgp + (1.0 - fbg)*params[i];
        p.push_back(this_p);
    }
    return dmultinom(x, p)/log2(exp(1));
}

void dll_mix_all(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
 
    int ndim = data_i.at("ndim");
    double tot = data_d.at("tot");

    static char buf[50];
    static string bufstr;
    double fbg = params[params.size()-1];
    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        double this_x = data_d.at(bufstr);
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgp = data_d.at(bufstr);
        double this_p = fbg * bgp + (1.0 - fbg)*params[i];
        double dLL_dp = this_x/this_p - tot;
        //results[results.size()-1] += bgp*dLL_dp;
        results[i] += (1-bgp)*dLL_dp;
    }
}

double infer_p_all(mixtureModel& mod,
    vector<string>& datakeys,
    vector<string>& bgpkeys,
    vector<vector<double> >& bgfrac,
    vector<double>& obsrow,
    int ndim,
    double pbg_exp,
    vector<double>& fgfracs){
    
    double tot = obsrow[ndim];
    vector<vector<double> > x1;
    for (int j = 0; j < ndim; ++j){
        x1.push_back(vector<double>{ obsrow[j] });
    } 
    
    vector<double> pg;
    for (int i = 0; i < ndim; ++i){
        pg.push_back(1.0 / (double)ndim);
    }    

    optimML::multivar_ml_solver mix(vector<double>{ }, ll_mix_all, dll_mix_all);
    mix.add_param_grp(pg);
    mix.add_one_param(pbg_exp);
    mix.constrain_01(ndim);
    //mix.set_delta(0.001);
    mix.add_data_fixed("ndim", ndim);
    mix.add_data_fixed("tot", tot);
    for (int k = 0; k < ndim; ++k){
        mix.add_data(datakeys[k], x1[k]);
        mix.add_data(bgpkeys[k], bgfrac[k]);
    }    
    mix.solve();
    double psum = 0.0;
    double fbg = mix.results[mix.results.size()-1];
    double psum2 = 0.0;
    for (int i = 0; i < ndim; ++i){
        fgfracs.push_back(mix.results[i]);
    }
    return fbg;
}

double infer_p(mixtureModel& mod,
    map<int, pair<int, int> >& ddoub,
    map<int, int>& ddim,
    int assn,
    vector<string>& datakeys,
    vector<string>& bgpkeys,
    vector<vector<double> >& bgfrac,
    vector<double>& obsrow,
    int ndim){
    
    double tot = obsrow[ndim];
    vector<vector<double> > x1;
    for (int j = 0; j < ndim; ++j){
        x1.push_back(vector<double>{ obsrow[j] });
    } 
    if (ddoub.count(assn) == 0){
        int dim = ddim[assn];
        
        optimML::brent_solver mix(ll_mix, dll_mix);
        mix.add_data_fixed("ndim", ndim);
        mix.constrain_01();
        vector<int> idx1{ dim };
        mix.add_data("idx1", idx1);
        for (int k = 0; k < ndim; ++k){
            mix.add_data(datakeys[k], x1[k]);
            mix.add_data(bgpkeys[k], bgfrac[k]);
        }    
        double res = mix.solve(0, 1);
        return res;
    }
    else{
        int dim1 = ddim[ddoub[assn].first];
        int dim2 = ddim[ddoub[assn].second];
        optimML::brent_solver mix2(ll_mix, dll_mix);
        mix2.add_data_fixed("ndim", ndim);
        mix2.constrain_01();
        vector<int> idxes1{ dim1 };
        vector<int> idxes2{ dim2 };
        mix2.add_data("idx1", idxes1);
        mix2.add_data("idx2", idxes2);
        vector<double> mixprop{ mod.shared_params[ndim*2 + dim1] / 
            (mod.shared_params[ndim*2 + dim1] + mod.shared_params[ndim*2 + dim2]) };
        mix2.add_data("mixprop", mixprop);
        for (int k = 0; k < ndim; ++k){
            if (mod.weights[k] > 0.0){
                mix2.add_data(datakeys[k], x1[k]);
                mix2.add_data(bgpkeys[k], bgfrac[k]);
            }
        }
        double res2 = mix2.solve(0,1);
        return res2;
    }
}

double get_cell_bg(vector<double>& obsrow,
    set<int>& assns,
    vector<double>& model_means,
    vector<double>& bgfracs,
    int ndim){
    
    // Infer background prop
    optimML::brent_solver mix(ll_mix_multi, dll_mix_multi);
    
    double thistot = 0.0;
    for (set<int>::iterator a = assns.begin(); a != assns.end(); ++a){
        thistot += model_means[*a];
        //onsizes[*a];
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
    if (!mix.root_found){
        return this_offtarget;
        //return -1;
    } 
    else{
        double n_bg = p_inferred * obsrow[obsrow.size()-1];
        return n_bg;
    }
}

bool assign_cell2(vector<double>& bgprops,
    vector<double>& sizes,
    double bg_mean,
    double bg_var,
    vector<double>& obsrow,
    vector<double>& obs2row,
    int ndim,
    pair<double, double>& a_b,
    bool add_bg,
    set<int>& result,
    double& llr,
    double& bgcount){
    
    /*
    double bgmeanlog = -1;
    double bgsdlog = -1;
    if (bg_mean != -1 && bg_sd != -1){
        double bgv = bg_sd*bg_sd;
        bgmeanlog = log(bg_mean) - log(bgv/(bg_mean*bg_mean) + 1)/2.0;
        bgsdlog = sqrt(log(bgv/(bg_mean*bg_mean) + 1.0));
    }
    */
    
    // Get expected percent background counts from overall count
    double tot = obsrow[obsrow.size()-1];
    double pbg = a_b.first * log(tot) + a_b.second;
    pbg = exp(pbg)/(exp(pbg) + 1);

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

    for (int i = 0; i < dim_enrich_sort.size(); ++i){
        int j = dim_enrich_sort[i].second;
        if (sizes[j] <= 0){
            continue;
        }

        included.insert(j);
        
        //mean += mod.dists[j].params[1][0];
        //sd += mod.dists[j].params[1][1];
        //tot_counts += mod.shared_params[ndim*2 + j];
        
        tot_counts += sizes[j];
        // Assume exponential variance (and use scaling factor to prevent overflow)
        var += pow(sizes[j], 2)/1e6;
        //var += vars[j];
        
        tot_fg += obsrow[j] * (1-bgprops[j]);

        vector<double> p;
        for (int k = 0; k < ndim; ++k){
            //if (mod.weights[k] > 0){
            if (sizes[k] > 0){
                double p_this;
                if (included.find(k) != included.end()){
                    p_this = bgprops[k] * pbg + (1.0 - pbg) * (sizes[k]/tot_counts);
                    //p_this = mod.shared_params[k] * pbg + (1.0 - pbg) * (mod.shared_params[ndim*2 + k]/tot_counts);
                }
                else{
                    p_this = bgprops[k] * pbg;
                    //p_this = mod.shared_params[k] * pbg;
                }
                p.push_back(p_this);
            }
        }
        double ll = dmultinom(obsrowclean, p);
        
        if (bg_mean != -1 && bg_var != -1){
            double tot_exp = tot_counts + bg_mean;
            //double tot_sd = sqrt(var) * sqrt(1e-6);
            var += bg_var/1e6;
            
            // Use log-normal distribution
            double m2 = pow(tot_exp, 2)/1e6;
            double sigma = sqrt(log(var/m2 + 1.0));
            double mu = log(tot_exp) - log(var/m2 + 1.0)/2.0;
            ll += dnorm(log(obsrow[obsrow.size()-1]), mu, sigma);

            //double tot_exp_v = var + bg_sd*bg_sd;
            //double tot_exp_log = log(tot_exp) - log(tot_exp_v/(tot_exp*tot_exp) + 1)/2.0;
            //double tot_sd_log = sqrt(log(tot_exp_v/(tot_exp*tot_exp) + 1.0));
            //ll += dnorm(log(obsrow[obsrow.size()-1]), tot_exp_log, tot_sd_log); 
            
            //ll += dexp(obsrow[obsrow.size()-1], 1.0/tot_exp);

            //double cbg = get_cell_bg(obsrow, included, sizes, bgprops, ndim);
            //ll += dnorm(log(cbg), bg_mean, bg_sd);

            //ll += dnbinom(cbg, bg_mean, bg_sd);
                
            //ll += dnorm(tot_fg, tot_counts, sqrt(var));
           
            //pair<double, double> muphi = nbinom_moments(tot_counts + bg_mean, var + bg_sd*bg_sd);
            //ll += dnbinom(obsrow[obsrow.size()-1], muphi.first, muphi.second);
            
            //ll += dnorm(obsrow[obsrow.size()-1], tot_counts + bg_mean, sqrt(var + bg_sd*bg_sd));
        }
        lls.push_back(make_pair(-ll, i));
    }
    
    if (add_bg){
        vector<double> params_bg;
        for (int j = 0; j < ndim; ++j){
            //if (mod.weights[j] > 0){
            if (sizes[j] > 0){
                params_bg.push_back(pbg * bgprops[j]);
                //params_bg.push_back(pbg * mod.shared_params[j]);
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

            //ll += dnbinom(cbg, bg_mean, bg_sd);
            //ll += dnorm(log(cbg), bg_mean, bg_sd);
            //ll += dnorm(log(cbg), bgmeanlog, bgsdlog);
            
            //ll += dexp(cbg, 1.0/bg_mean);

            //pair<double, double> muphi = nbinom_moments(bg_mean, bg_sd*bg_sd);
            //ll += dnbinom(cbg, muphi.first, muphi.second);
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

bool assign_cell(mixtureModel& mod,
    double bg_mean,
    double bg_sd,
    vector<double>& obsrow,
    vector<double>& obs2row,
    int ndim,
    pair<double, double>& a_b,
    map<int, pair<int, int> >& ddoub,
    map<int, int>& ddim,
    bool add_bg,
    set<int>& result,
    double& llr){
    
    double tot = obsrow[obsrow.size()-1];
    double pbg = a_b.first * log(tot) + a_b.second;
    pbg = exp(pbg)/(exp(pbg) + 1);

    vector<pair<double, int> > lls;
    vector<string> names;
    vector<int> nums;
    
    vector<double> obsrowclean;

    // Sort in decreasing order of enrichment of p over expectation in background
    vector<pair<double, int> > dim_enrich_sort;
    for (int j = 0; j < ndim; ++j){
        if (mod.weights[j] > 0){
            double p = obsrow[j]/tot - mod.shared_params[j];
            dim_enrich_sort.push_back(make_pair(-p, j));
            obsrowclean.push_back(obsrow[j]);
        }
    }
    sort(dim_enrich_sort.begin(), dim_enrich_sort.end());
    double tot_counts = 0;
    set<int> included;
    double mean = 0.0;
    double sd = 0.0;
    vector<double> llsize;

    for (int i = 0; i < dim_enrich_sort.size(); ++i){
        int j = dim_enrich_sort[i].second;
        included.insert(j);
        mean += mod.dists[j].params[1][0];
        sd += mod.dists[j].params[1][1];
        tot_counts += mod.shared_params[ndim*2 + j];
        vector<double> p;
        for (int k = 0; k < ndim; ++k){
            if (mod.weights[k] > 0){
                double p_this;
                if (included.find(k) != included.end()){
                    p_this = mod.shared_params[k] * pbg + (1.0 - pbg) * (mod.shared_params[ndim*2 + k]/tot_counts);
                }
                else{
                    p_this = mod.shared_params[k] * pbg;
                }
                p.push_back(p_this);
            }
        }
        double ll = dmultinom(obsrowclean, p);
        pair<double, double> nbpar = nbinom_moments(mean, sd*sd);
        double ll2 = dnbinom(obsrow[obsrow.size()-1], nbpar.first, nbpar.second);
        llsize.push_back(ll2);
        lls.push_back(make_pair(-ll, i));
    }

    vector<double> params_bg;
    for (int j = 0; j < ndim; ++j){
        if (mod.weights[j] > 0){
            params_bg.push_back(pbg * mod.shared_params[j]);
        }
    }
    double ll = dmultinom(obsrowclean, params_bg);
    pair<double, double> nbpar = nbinom_moments(bg_mean, bg_sd*bg_sd);
    double llsize_bg = dnbinom(obsrow[obsrow.size()-1], nbpar.first, nbpar.second);
    //ll += dnorm(obsrow[obsrow.size()-1], bg_mean, bg_sd);
    double ll_bg = ll;
    lls.push_back(make_pair(-ll, -1));

    sort(lls.begin(), lls.end());
    llr = -lls[0].first - -lls[1].first;
    int chosen = lls[0].second;
    if (chosen == -1){
        fprintf(stderr, "elim %ld, %f: %f %f | %f %f\n", result.size(), obsrow[obsrow.size()-1],
            -lls[1].first, llsize[lls[1].second], ll_bg, llsize_bg);
        return false;    
    }
    else{
        for (int i = 0; i <= chosen; ++i){
            result.insert(dim_enrich_sort[i].second);
        }
        if (result.size() > 2){
            fprintf(stderr, "%ld, %f: %f %f | %f %f\n", result.size(), obsrow[obsrow.size()-1], 
                -lls[0].first, llsize[chosen],
               ll_bg, llsize_bg); 
        }
        return true;
    }

    /*
    for (int j = 0; j < mod.n_components; ++j){
        // Skip 0-weight components
        if (mod.weights[j] > 0.0){
            int dim1 = -1;
            int dim2 = -1;
            if (ddoub.count(j) > 0){
                dim1 = ddim[ddoub[j].first];
                dim2 = ddim[ddoub[j].second];
            }
            else{
                dim1 = ddim[j];
            }
            vector<double> params;
            for (int k = 0; k < ndim; ++k){
                if (dim2 == -1){
                    if (k == dim1){
                        params.push_back(pbg*mod.shared_params[k] + (1.0 - pbg));
                    }
                    else{
                        params.push_back(pbg*mod.shared_params[k]);
                    }
                }
                else{
                    double m1 = mod.shared_params[ndim*2 + dim1];
                    double m2 = mod.shared_params[ndim*2 + dim2];
                    double m = m1/(m1+m2);
                    if (k == dim1){
                        params.push_back(pbg*mod.shared_params[k] + (1.0-pbg)*m);
                    }
                    else if (k == dim2){
                        params.push_back(pbg*mod.shared_params[k] + (1.0-pbg)*(1.0-m));
                    }
                    else{
                        params.push_back(pbg*mod.shared_params[k]);
                    }
                }
            }
            double ll_this = dmultinom(obs2row, params);

            //ll_this += log2(mod.weights[j]);
            lls.push_back(make_pair(-ll_this, j));
        }
    }
    
    if (add_bg){
        // Include a component for background-only
        vector<double> params_blank;
        for (int j = 0; j < ndim; ++j){
            params_blank.push_back(mod.shared_params[j]);
        }
        double ll_this = dmultinom(obs2row, params_blank);
        lls.push_back(make_pair(-ll_this, -1));
    }

    sort(lls.begin(), lls.end());
    ll = -lls[0].first;
    if (lls[0].second == -1){
        llr = 0.0;
        return -1;
    }
    else{
        llr = -lls[0].first - -lls[1].first;
        if (llr <= 0){
            return -1;
        }
        return lls[0].second;
    }
    */
}

double percentile(vector<double>& vec, double quant){
    sort(vec.begin(), vec.end());
    double tot = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        tot += vec[i];
    }
    double runtot = 0.0;
    for (int i = 0; i < vec.size(); ++i){
        runtot += vec[i];
        if (runtot >= tot*quant){
            return vec[i];
        }
    }
    return vec[vec.size()-1];
}

/**
 * Finds the point of intersection between two gaussians.
 */
pair<double, double> int2gauss(double mu1, double mu2, double sigma1, double sigma2){
    if (mu2 < mu1){
        // Flip.
        double tmp = mu2;
        double tmp2 = sigma2;
        mu2 = mu1;
        sigma2 = sigma1;
        mu1 = tmp;
        sigma1 = tmp2;
    }
    double sig21 = sigma1*sigma1;
    double sig22 = sigma2*sigma2;
    double mu21 = mu1*mu1;
    double mu22 = mu2*mu2;
    
    double a = (-1.0/sig21 + 1.0/sig22)/2;
    double b = -mu2/sig22 + mu1/sig21;
    double c = (mu22/sig22 - mu21/sig21 + log(sig22/sig21))/2;
    
    // Solve ax^2 + bx + c
    double sol1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    double sol2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
    return make_pair(sol1, sol2);
}

double ov2gauss(double mu1, double mu2, double sigma1, double sigma2){
    
    if (mu2 < mu1){
        // Flip
        double tmp1 = mu1;
        mu1 = mu2;
        mu2 = tmp1;
        double tmp2 = sigma1;
        sigma1 = sigma2;
        sigma2 = tmp2;
    }
    pair<double, double> intpts = int2gauss(mu1, mu2, sigma1, sigma2);
    double ov1 = (1.0 - pnorm(intpts.first, mu1, sigma1) + pnorm(intpts.first, mu2, sigma2));
    double ov2 = (1.0 - pnorm(intpts.second, mu1, sigma1) + pnorm(intpts.second, mu2, sigma2));
    if (ov1 >= 1.0 || ov1 <= 0.0){
        return ov2;
    }
    else{
        return ov1;
    }
}


void assign_ids2(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    robin_hood::unordered_map<unsigned long, double>& cell_bg,
    vector<string>& labels,
    bool filt){
    
    vector<vector<double> > obscol;
    map<string, int> label2idx;
    vector<double> meanscol;

    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
        vector<double> v;
        obscol.push_back(v);
        meanscol.push_back(0.0);
    }
    
    vector<vector<double> > obs;
    vector<unsigned long> bcs;
    vector<double> tots;
    
    // We'll dump totals into the "obs" vectors. "obs2" will not include totals.
    vector<vector<double> > obs2;

    // Prepare data
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        vector<double> row;
        for (int i = 0; i < labels.size(); ++i){
            row.push_back(0.0);
        }
        double tot = 0.0;
        for (vector<pair<int, string> >::iterator y = x->second.begin(); y != x->second.end(); ++y){
            row[label2idx[y->second]] += (double)-y->first;
            tot += (double)-y->first;
            
            //double val = log((double)-y->first + 1);
            double val = (double)-y->first + 1;
            obscol[label2idx[y->second]].push_back(val);
            meanscol[label2idx[y->second]] += val;
        }
        // Skip rows with 0 counts
        if (tot > 0.0){
            tots.push_back(tot);
            obs2.push_back(row);
            row.push_back(tot);
            obs.push_back(row);
            bcs.push_back(x->first);
        }
    }
    
    // Get mean & SD sizes
    vector<double> model_means;
    
    vector<double> lomeans;
    vector<double> lophis;
    vector<double> himeans;

    for (int i = 0; i < labels.size(); ++i){
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
        fprintf(stderr, "%s\n", labels[i].c_str());
        fprintf(stderr, "  %f %f\n", mod.dists[0].params[0][0], 1.0/mod.dists[1].params[0][0]);
        lomeans.push_back(mod.dists[0].params[0][0]);
        lophis.push_back(mod.dists[0].params[0][1]);
        himeans.push_back(1.0/mod.dists[1].params[0][0]);
    }
    
    pair<double, double> lomv = welford(lomeans);
    pair<double, double> himv = welford(himeans);
    fprintf(stderr, "== MEAN %f %f ==\n", lomv.first, himv.first);
    fprintf(stderr, "== SD   %f %f ==\n", sqrt(lomv.second), sqrt(himv.second));
    for (int i = 0; i < labels.size(); ++i){
        fprintf(stderr, "%f\n", himeans[i]);
        double d1 = dnorm(himeans[i], lomv.first, sqrt(lomv.second));
        double d2 = dnorm(himeans[i], himv.first, sqrt(himv.second));
        fprintf(stderr, "%f %f\n", d1, d2);
        string filtstr = "";
        if (d2 - d1 < 1){
            fprintf(stderr, "REMOVE %s\n", labels[i].c_str());
            model_means.push_back(-1);
            lomeans[i] = -1;
            lophis[i] = -1;
            himeans[i] = -1;
        }
        else{
            model_means.push_back(himeans[i]);

        }
    }

    double p_thresh = 0.5;
    
    vector<double> bgtots;
    vector<double> bgcounts;
    
    for (int i = 0; i < labels.size(); ++i){
        bgtots.push_back(0.0);
        bgcounts.push_back(0.0);
        vector<double> v;
    }

    vector<double> logtots;
    vector<double> logitpercs;
    vector<unsigned long> logtots_bc;
   
    double bgcount = 0.0;
    
    vector<double> bgcount_all;
    
    map<unsigned long, set<int> > assnorig;
 
    vector<vector<double> > counts_all;
    vector<vector<double> > weights_all;
    for (int i = 0; i < labels.size(); ++i){
        vector<double> v;
        counts_all.push_back(v);
        weights_all.push_back(v);
    }

    for (int i = 0; i < obs.size(); ++i){
        map<int, double> probs;
        
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

            probs.insert(make_pair(j, p_high));
            if (p_high > p_thresh){
                counts_all[j].push_back(obs[i][obs[i].size()-1]);
                chosen.insert(j);
            }
            else{
                countbg += obs[i][j];
                bgtots[j] += obs[i][j];
            }
        }

        assnorig.insert(make_pair(bcs[i], chosen));

        bgcount += countbg;
        bgcount_all.push_back(countbg);
        
        if (count_tot > 0 && countbg > 0 && countbg < count_tot){
            logtots.push_back(log(count_tot));
            logitpercs.push_back(log(countbg/count_tot) - log(1-countbg/count_tot));
            logtots_bc.push_back(bcs[i]);
        }
    }

    // Get background fracs
    for (int i = 0; i < labels.size(); ++i){
        bgtots[i] /= bgcount;
        fprintf(stderr, "BG %s: %f\n", labels[i].c_str(), bgtots[i]);
    }
    
    for (int i = 0; i < labels.size(); ++i){
        fprintf(stderr, "MODEL %s: %f\n", labels[i].c_str(), model_means[i]);
    }

    pair<double, double> slope_int = linreg(logtots, logitpercs);
    fprintf(stderr, "A,B = %f %f\n", slope_int.first, slope_int.second);
    
    bgcount_all.clear();
    
    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr;
        double bgcount; 
        assign_cell2(bgtots, model_means, -1, -1, obs[i], obs2[i],   
            labels.size(), slope_int, true, assns, llr, bgcount);
        if (bgcount != -1){
            bgcount_all.push_back(round(bgcount) + 1);    
            cell_bg.emplace(bcs[i], bgcount / obs[i][obs[i].size()-1]);
        }
    }
    
    
    pair<double, double> muvarbg = welford(bgcount_all);
    
    mixtureModel bgmod;

    if (filt){ 
        /* 
        // Test for bimodality of log-transformed background size dist
        vector<double> bgcount_log;
        for (int i = 0; i < bgcount_all.size(); ++i){
            bgcount_log.push_back(log(bgcount_all[i] + 1.0));
        }
        pair<double, double> bglog_mv = welford(bgcount_log);
        double bglog_s = sqrt(bglog_mv.second);
        double skew = 0.0;
        double kurt = 0.0;
        double n = (double)bgcount_all.size();
        double w = 1.0/(n-1);
        for (int i = 0; i < bgcount_all.size(); ++i){
            double xms = (bgcount_log[i] - bglog_mv.first)/bglog_s;
            double xms3 = xms*xms*xms;
            double xms4 = xms3*xms;
            skew += w*xms3;
            kurt += w*xms4;
        }
        kurt -= 3.0;

        // Calculate bimodality coefficient
        // https://www.frontiersin.org/journals/psychology/articles/10.3389/fpsyg.2013.00700/full
        double frac = ((n-1)*(n-1))/((n-2)*(n-3));
        double bim_coef = (skew*skew + 1.0)/(kurt + 3*frac);
        fprintf(stderr, "BIMODALITY COEF %f\n", bim_coef);

        if (bim_coef <= 0.6){
            // Don't bother filtering.
            filt = false;
        }
        else{
        */
        if (true){
            pair<double, double> muvar_all = welford(bgcount_all);
            fprintf(stderr, "MEAN %f SD %f\n", muvar_all.first, sqrt(muvar_all.second));
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
            bgmod.init(bgdists, w);
            vector<vector<double> > bgobs;
            for (int i = 0; i < bgcount_all.size(); ++i){
                bgobs.push_back(vector<double>{ bgcount_all[i] });
            }
            bgmod.fit(bgobs);
            bgmod.print();
            fprintf(stderr, "%f %f\n", bgmod.dists[0].params[0][0], 1.0/bgmod.dists[1].params[0][0]);
            
            // Don't filter if we'll lose too much data
            // or if means are too close together.

            double mean1 = bgmod.dists[0].params[0][0];
            double mean2 = 1.0/bgmod.dists[1].params[0][0];
            double p = pnbinom(mean2, bgmod.dists[0].params[0][0], bgmod.dists[0].params[0][1]);
            fprintf(stderr, "p = %f\n", p);
            if (pnbinom(mean2, bgmod.dists[0].params[0][0], bgmod.dists[0].params[0][1]) < 0.9){
                filt = false;
            }
        }
    }

    int count_filt = 0;

    for (int i = 0; i < obs.size(); ++i){
        set<int> assns;
        double llr; 
        double bgcount;
        bool assigned = assign_cell2(bgtots, model_means, muvarbg.first, muvarbg.second, 
            obs[i], obs2[i],   
            labels.size(), slope_int, true, assns, llr, bgcount);
        
        if (assigned && filt){
            vector<double> row{ bgcount };
            double ll1 = bgmod.dists[0].loglik(row);
            double ll2 = bgmod.dists[1].loglik(row);
            if (ll2 > ll1){
                assigned = false;
                count_filt++;
            }    
        }
        if (assigned){
            assn.emplace(bcs[i], assns);
            assn_llr.emplace(bcs[i], llr);
        }
    }
    
    fprintf(stderr, "Filtered %d of %ld assignments\n", count_filt, obs.size());

}

/**
 * Main function to determine the identities of cells. Performs EM fitting, solves
 * for a set of common parameters at each step, then (given first-round assignments)
 * infers the relationship between a cell's total tag counts and its percent 
 * background/ambient tags, and then uses all solved parameters and this relationship
 * to infer a final identity per cell.
 */
void assign_ids(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn,
    robin_hood::unordered_map<unsigned long, bool>& assn_doublet,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    vector<string>& labels,
    double filt_prob,
    int sampsize){
    
    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }
    
    vector<vector<double> > obs;
    vector<vector<double> > obs_subsamp;
    vector<unsigned long> bcs;
    vector<double> tots;
    
    // We'll dump totals into the "obs" vectors. "obs2" will not include totals.
    vector<vector<double> > obs2;

    double samp_prob = -1;
    if (sampsize > 0){
        samp_prob = (double)sampsize/(double)bc_ms_counts.size();
    }
    srand(time(NULL));

    // Prepare data
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        vector<double> row;
        for (int i = 0; i < labels.size(); ++i){
            row.push_back(0.0);
        }
        double tot = 0.0;
        for (vector<pair<int, string> >::iterator y = x->second.begin(); y != x->second.end(); ++y){
            row[label2idx[y->second]] += (double)-y->first;
            tot += (double)-y->first;
        }
        tots.push_back(tot);
        obs2.push_back(row);
        row.push_back(tot);
        obs.push_back(row);
        bcs.push_back(x->first);
        if (samp_prob <= 0 || ((double) rand() / (RAND_MAX)) < samp_prob){
            obs_subsamp.push_back(row);
        }
    }
    
    // Now prepare distributions for EM fitting

    pair<double, double> totmuvar = welford(tots);
    double totmean = totmuvar.first;
    double totsd = sqrt(totmuvar.second);
    double totsd_doub = sqrt(totmuvar.second*2);

    double off_target = 0.1;
    map<int, int> dist2dim;
    map<int, pair<int, int> > dist2doublet;
    vector<mixtureDist> dists;
    vector<double> shared_params;
    
    for (int i = 0; i < labels.size(); ++i){
        // Fraction in background 
        shared_params.push_back(1.0 / (double)labels.size());
    }
    for (int i = 0; i < labels.size(); ++i){
        // Proportion composed of background
        shared_params.push_back(off_target);
    }

    for (int i = 0; i < labels.size(); ++i){
        // Mean size (read count)
        shared_params.push_back(totmean);

        vector<double> params;
        for (int j = 0; j < labels.size(); ++j){
            if (j == i){
                params.push_back(1.0 - off_target);
            }
            else{
                params.push_back(off_target / (double)(labels.size()-1));
            }
        }
        
        // Fit a mixture model both over the proportions of reads per cell belonging to
        // each label, and to the total counts per cell (i.e. doublets should have 
        // higher counts on average)

        mixtureDist dist(vector<string>{ "multinomial", "normal" }, 
            vector<vector<double> >{ params, vector<double>{ totmean, totsd } });
        dist.set_num_inputs(0, labels.size());

        dist.name = labels[i];
        int didx = dists.size();
        ddim.insert(make_pair(didx, i));
        dists.push_back(dist);
    }
    
    // Create distributions representing doublets
    map<pair<int, int>, int> doub2dist;

    for (int i = 0; i < labels.size(); ++i){
        for (int j = i + 1; j < labels.size(); ++j){
            vector<double> params2;
            for (int k = 0; k < labels.size(); ++k){
                if (k == i || k == j){
                    params2.push_back((1.0-off_target)/2.0);
                }
                else{
                    params2.push_back(off_target / (double)(labels.size()-2));
                }
            }
            mixtureDist dist2(vector<string>{ "multinomial", "normal" }, 
                vector<vector<double> >{ params2, vector<double>{ 2*totmean, totsd_doub }});
            dist2.set_num_inputs(0, labels.size());

            string name;
            if (labels[i] < labels[j]){
                name = labels[i] + "+" + labels[j];
            }
            else{
                name = labels[j] + "+" + labels[i];
            }
            dist2.name = name;
            int didx = dists.size();
            ddoub.insert(make_pair(didx, make_pair(i, j)));
            doub2dist.insert(make_pair(make_pair(i, j), didx));
            dists.push_back(dist2);
        }
    }
    
    // Include a distribution just representing "background" cells?
    bool add_bg = false;

    if (add_bg){
        vector<double> params_bg;
        for (int i = 0; i < labels.size(); ++i){
            params_bg.push_back(1.0 / (double)labels.size());
        }
        mixtureDist dist(vector<string>{ "multinomial", "normal" },
            vector<vector<double> >{ params_bg, vector<double>{ totmean, totsd } });
        dist.set_num_inputs(labels.size());
        dist.name = "BG";
        dists.push_back(dist);
    }
    
    // Fit the model, with a function to reconcile shared parameters across all
    // distributions after each "M" step

    mixtureModel mod(dists);
    mod.set_callback(callback_func, shared_params);
    mod.fit(obs_subsamp);
    
    fprintf(stderr, "=== Ambient tag count profile ===\n"); 
    for (int i = 0; i < labels.size(); ++i){
        if (mod.weights[i] > 0){
            fprintf(stderr, "%s\t%.3f\n", labels[i].c_str(), mod.shared_params[i]);
        }
    }
    fprintf(stderr, "=== Mean percent ambient counts per cell type ===\n");
    for (int i = 0; i < labels.size(); ++i){
        if (mod.weights[i] > 0){
            fprintf(stderr, "%s\t%.3f\n", labels[i].c_str(), mod.shared_params[labels.size() + i]);
        }
    }
    fprintf(stderr, "=== Mean total tag counts per cell type ===\n");
    for (int i = 0; i < labels.size(); ++i){
        if (mod.weights[i] > 0){
            fprintf(stderr, "%s\t%d\n", labels[i].c_str(), 
                (int)round(mod.shared_params[labels.size()*2 + i]));
        }
    }

    // Next step: fit a function to predict the percent of a cell's barcode
    // counts consisting of background noise, given the total reads counted
    // for that cell

    // First, we will infer (via multinomial MLE) the likeliest percent background
    // per cell, given the current model.

    // Then we will use linear regression to fit logit(percent background) to 
    // log(total counts) across all cells

    vector<string> datakeys;
    vector<string> bgpkeys;
    char buf[50];
    string bufstr;
    for (int i = 0; i < labels.size(); ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        datakeys.push_back(bufstr);
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        bgpkeys.push_back(bufstr);
    }
    
    // Remove 0-weight dists
    double perctot = 0.0;
    for (int j = 0; j < labels.size(); ++j){
        if (mod.weights[j] > 0.0){
            perctot += mod.shared_params[j];
        }
    }
    for (int j = 0; j < labels.size(); ++j){
        mod.shared_params[j] /= perctot;
    }
    
    vector<vector<double> > bgfrac;
    for (int j = 0; j < labels.size(); ++j){
        bgfrac.push_back(vector<double>{ mod.shared_params[j] });
    }
    vector<double> totcounts;
    vector<double> pbgs;
    
    map<unsigned long, double> bc2p;

    for (int i = 0; i < obs.size(); ++i){
        double p = infer_p(mod, ddoub, ddim, mod.assignments[i],
            datakeys, bgpkeys, bgfrac, obs[i], labels.size());
        totcounts.push_back(log(obs[i][obs[i].size()-1]));
        pbgs.push_back(log(p) - log(1.0-p));      
    } 
   
    // Now re-visit all cells and identify them, using a tweaked model that estimates
    // percent background using this function, and which leaves out the count component
    // (we will be using multinomial probabilities only)
        
    // Find best fit relationship
    pair<double, double> a_b = linreg(totcounts, pbgs);
    for (int i = 0 ; i < totcounts.size(); ++i){
        double pred = a_b.first * totcounts[i] + a_b.second;
    } 
    fprintf(stderr, "=== Total/ambient tag count relationship ===\n");
    fprintf(stderr, "b = (%.4f*c^%.4f)/(%.4f*c^%.4f + 1)\n",
        exp(a_b.second), a_b.first, exp(a_b.second), a_b.first);
    fprintf(stderr, "  b = percent of tag counts from ambient\n");
    fprintf(stderr, "  c = total tag counts per cell\n");    
    
    /*
    fprintf(stdout, "barcode\tbg");
    for (int i = 0; i < labels.size(); ++i){
        fprintf(stdout, "\t%s", labels[i].c_str());
    }
    fprintf(stdout, "\n");
    for (int i = 0; i < obs.size(); ++i){
        double pred = a_b.first * log(obs[i][obs[i].size()-1]) + a_b.second;
        double percbg = exp(pred)/(exp(pred) + 1.0);
        string bcstr = bc2str(bcs[i]);
        
        vector<pair<double, int> > counts;
        double tot = obs[i][obs[i].size()-1];

        for (int j = 0; j < labels.size(); ++j){
            if (mod.weights[j] > 0.0){
                double bgcount = mod.shared_params[j] * percbg * tot;
                double count = obs[i][j] - bgcount;
                if (count > 0){
                    counts.push_back(make_pair(-count, j));
                }
            }
        }
        sort(counts.begin(), counts.end());
        double llprev = 0.0;
        int truenum = -1;
        for (int num = 1; num < counts.size(); ++num){
            double meantrue = 0.0;
            double counttrue = 0.0;
            for (int x = 0; x < num; ++x){
                meantrue += -counts[x].first;
                counttrue++;
            }
            meantrue /= counttrue;
            double ll = 0.0;
            for (int x = 0; x < num; ++x){
                ll += dpois(-counts[x].first, meantrue); 
            }
            double meanbg = 0.0;
            double countbg = 0.0;
            for (int x = num; x < counts.size(); ++x){
                meanbg += -counts[x].first;
                countbg++;
            }
            meanbg /= countbg;
            for (int x = num; x < counts.size(); ++x){
                ll += dpois(-counts[x].first, meanbg);
            }
            if (llprev != 0 && ll < llprev){
                break;
            }
            else{
                truenum = num;
                llprev = ll;
            }
        }
        if (counts.size() == 1){
            truenum = 1;
        }
        set<string> names;
        for (int x = 0; x < truenum; ++x){
            names.insert(labels[counts[x].second]);
        }
        string namestr = "";
        for (set<string>::iterator n = names.begin(); n != names.end(); ++n){
            if (namestr != ""){
                namestr += "+";
            }
            namestr += *n;
        }
        if (truenum == -1){
            fprintf(stderr, "percbg %f\n", percbg);
            fprintf(stderr, "raw\n");
            for (int j = 0; j < obs[i].size()-1; ++j){
                fprintf(stderr, "%f\n", obs[i][j]);
            }
            fprintf(stderr, "COUNTS\n");
            for (int i = 0; i < counts.size(); ++i){
                fprintf(stderr, "%f %d\n ", counts[i].first, counts[i].second);
            }
            exit(0);
        }
        fprintf(stdout, "%s\t%d\t%s\n", bcstr.c_str(), truenum, namestr.c_str());

    }    
    exit(0);
    */
    
    /*
    for (int i = 0; i < obs.size(); ++i){
        double pred = a_b.first * log(obs[i][obs[i].size()-1]) + a_b.second;
        pred = exp(pred)/(exp(pred) + 1.0);
        string bcstr = bc2str(bcs[i]);
        fprintf(stdout, "%s", bcstr.c_str());
        vector<double> fgfracs;
        double bgf = infer_p_all(mod, datakeys, bgpkeys, bgfrac, obs[i], labels.size(), pred, fgfracs);
        fprintf(stdout, "\t%f", bgf);
        for (int j = 0; j < fgfracs.size(); ++j){
            fprintf(stdout, "\t%f", fgfracs[j]);
        }
        fprintf(stdout, "\n");

    }
    exit(0);
    */
    
    vector<double> bgsizes;
    for (int i = 0; i < obs.size(); ++i){
        double tot = obs[i][obs[i].size()-1];
        double pbg = a_b.first * log(tot) + a_b.second;
        pbg = exp(pbg)/(exp(pbg) + 1);
        bgsizes.push_back(pbg*tot);
    }
    pair<double, double> bgmuvar = welford(bgsizes);
    double bg_mean = bgmuvar.first;
    double bg_sd = sqrt(bgmuvar.second);

    vector<vector<double> > filt_dat;
    vector<unsigned long> filt_dat_bc;
    
    fprintf(stderr, "BG: %f %f\n", bg_mean, bg_sd);

    for (int i = 0; i < obs.size(); ++i){
        double llr = 0.0;
        double ll = 0.0;
        int assn_this;

        set<int> assned;
        string bcx = bc2str(bcs[i]);
        fprintf(stderr, "%s\n", bcx.c_str());
        bool not_bg = assign_cell(mod, bg_mean, bg_sd, obs[i], obs2[i], labels.size(),
            a_b, ddoub, ddim, add_bg, assned, llr);
        string name = "";
        set<string> names;
        for (set<int>::iterator a = assned.begin(); a != assned.end(); ++a){
            names.insert(labels[*a]);
        }
        for (set<string>::iterator n = names.begin(); n != names.end(); ++n){
            if (name != ""){
                name += "+";
            }
            name += *n;
        }
        char s_d = 'S';
        if (assned.size() == 2){
            s_d = 'D';
        }
        else if (assned.size() > 2){
            s_d = 'M';
        }
        string bcstr = bc2str(bcs[i]);
        fprintf(stdout, "%s\t%s\t%c\t%f\n", bcstr.c_str(), name.c_str(), s_d, llr);

        continue;
        if (assn_this != -1){
            
            double p = infer_p(mod, ddoub, ddim, assn_this,
                datakeys, bgpkeys, bgfrac, obs[i], labels.size());
            
            double pred = a_b.first * log(obs[i][obs[i].size()-1]) + a_b.second;
            pred = exp(pred)/(exp(pred) + 1.0);
            
            vector<double> fdrow;
            fdrow.push_back(log(obs[i][obs[i].size()-1]));
            fdrow.push_back(log(pred)-log(p));
            filt_dat.push_back(fdrow);
            filt_dat_bc.push_back(bcs[i]);
            
            string bcstr = bc2str(bcs[i]);
            //fprintf(stdout, "%s\t%f\t%f\t%f\n", bcstr.c_str(), obs[i][obs[i].size()-1], p, pred);

            //totcounts.push_back(log(obs[i][obs[i].size()-1]));
            //pbgs.push_back(log(p) - log(1-p));
            //fprintf(stdout, "%f\n", p);
            
            assn.emplace(bcs[i], mod.dists[assn_this].name);
            assn_llr.emplace(bcs[i], llr);
            if (ddoub.count(assn_this) > 0){
                assn_doublet.emplace(bcs[i], true);
            }
            else{
                assn_doublet.emplace(bcs[i], false);
            }
        }
    }

    if (filt_prob > 0 && filt_prob < 1){
        vector<mixtureDist> dists_filt;
        vector<string> names{ "normal", "normal" };
        vector<double> pgood{ 0.0, 0.25 };
        vector<double> pbad{ -0.5, 0.25 };
        double meansize = 0;
        double meansize_count = 0;
        for (int i = 0; i < labels.size(); ++i){
            if (mod.weights[i] > 0){
                meansize += mod.shared_params[labels.size()*2 + i];
                meansize_count++;
            }
        }
        meansize /= meansize_count;
        vector<double> sgood{ log(meansize), log(meansize) };
        vector<double> sbad{ log(meansize + 1), log(meansize) };
        //mixtureDist dgood(names, vector<vector<double> >{ sgood, pgood }, vector<double>{ 0.5, 0.5 });
        //mixtureDist dbad(names, vector<vector<double> >{ sbad, pbad }, vector<double>{ 0.5, 0.5 });
        mixtureDist dgood("2dgauss", vector<double>{ log(meansize), 0.0, log(meansize), 0.25, 0.01 });
        mixtureDist dbad("2dgauss", vector<double>{ log(meansize) + 1, -0.5, log(meansize), 0.25, 0.01 }); 
        dists_filt.push_back(dgood);
        dists_filt.push_back(dbad);
        mixtureModel mod_filt(dists_filt);
        mod_filt.print();
        mod_filt.fit(filt_dat);
        mod_filt.print();
        
        int cfilt = 0;
        int idx = 0;
        
        for (vector<short>::iterator a = mod_filt.assignments.begin(); a != mod_filt.assignments.end(); ++a){
            if (*a == 1 && mod_filt.responsibility_matrix[idx][*a] > filt_prob){
                cfilt++;
                assn.erase(filt_dat_bc[idx]);
                assn_llr.erase(filt_dat_bc[idx]);
            }
            ++idx;
        }
        
        fprintf(stderr, "Filtered out %d of %ld cells (%.2f%%)\n", cfilt, obs.size(),
            100.0*(double)cfilt/(double)obs.size());
    }
}

/**
 * Spill cell -> identity assignments to disk (in .assignments file)
 */
void write_assignments(string filename, 
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    vector<string>& labels,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    string sep,
    const string& batch_id,
    bool sgrna){
    
    FILE* f_out = fopen(filename.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, set<int> >::iterator a = assn.begin(); 
        a != assn.end(); ++a){

        if (a->second.size() == 0){
            // Filtered/background. Omit.
            continue;
        }
        else{
            char type = 'S';
            if (a->second.size() == 2){
                type = 'D';
            }
            else if (a->second.size() > 2){
                type = 'M';
            }
            string name = "";
            set<string> names;
            for (set<int>::iterator i = a->second.begin(); i != a->second.end(); ++i){
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
            string bc_str = bc2str(a->first);
            if (batch_id != ""){
                bc_str += "-" + batch_id;
            }
            if (sgrna){
                fprintf(f_out, "%s\t%s\t%ld\t%f\n", bc_str.c_str(), name.c_str(), a->second.size(), 
                    assn_llr[a->first]);
            }
            else{
                fprintf(f_out, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type, assn_llr[a->first]); 
            }
        }
    }
    fclose(f_out);
}

void write_assignments_table(string filename,
    robin_hood::unordered_map<unsigned long, set<int> >& assn,
    vector<string>& labels,
    const string& batch_id){
    
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

    for (robin_hood::unordered_map<unsigned long, set<int> >::iterator a = assn.begin(); 
        a != assn.end(); ++a){
        
        string bcstr = bc2str(a->first);
        if (batch_id != ""){
            bcstr += "-" + batch_id;
        }
        fprintf(f_out, "%s", bcstr.c_str());
        for (vector<pair<string, int> >::iterator l = labelsort.begin(); l != labelsort.end(); ++l){
            if (a->second.find(l->second) != a->second.end()){
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
    
    fprintf(stderr, "check missing MS wells\n");    
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
    if (minrank_missing < toprank_present){
        fprintf(stderr, "WARNING: one or more un-named wells has higher counts than one or \
more named/expected wells. You may have mis-specified well-to-ID mappings in provided \
metadata. Please see the output well counts file to evaluate.\n");
    }

}

int main(int argc, char *argv[]) {    

    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"mapping", required_argument, 0, 'm'},
       {"whitelist", required_argument, 0, 'w'},
       {"barcodes", required_argument, 0, 'b'},
       {"cell_barcodes", required_argument, 0, 'B'},
       {"mismatches", required_argument, 0, 'M'},
       {"filt", no_argument, 0, 'f'},
       {"comma", no_argument, 0, 'c'},
       {"batch_id", required_argument, 0, 'i'},
       {"sgRNA", no_argument, 0, 'g'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string output_prefix = "";
    vector<string> read1fn;
    vector<string> read2fn;
    string mapfn = "";
    string wlfn = "";
    string cell_barcodesfn;
    int mismatches = 2;
    bool filt = false;
    string sep = "+";
    string batch_id = "";
    bool sgrna = false;

    fprintf(stderr, "%s\n", STRINGIZE(PROJ_ROOT));
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
    while((ch = getopt_long(argc, argv, "o:i:1:2:m:M:w:b:B:gfch", 
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
            case 'i':
                batch_id = optarg;
                break;
            case '1':
                read1fn.push_back(optarg);
               break;
            case '2':
                read2fn.push_back(optarg);
                break;
            case 'm':
                mapfn = optarg;
                break;
            case 'M':
                mismatches = atoi(optarg);
                break;
            case 'w':
                wlfn = optarg;
                break;
            case 'b':
                barcodesfn = optarg;
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
    
    if (sgrna && sep == "+"){
        sep = ",";
    }

    // Data structure to store counts
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > > bc_ms_counts;
    
    // How many possible MULTI-seq labels are there?
    int n_labels;
    vector<string> labels;

    string countsfilename = output_prefix + ".counts";
    if (file_exists(countsfilename)){
        fprintf(stderr, "Loading previously-computed counts (delete file to avoid this \
next time)...\n");
        n_labels = load_counts(countsfilename, bc_ms_counts, labels);
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
        
        // Load MULTIseq barcode -> well ID mapping
        map<unsigned long, string> bc2well;
        vector<string> bclist;
        int ms_bc_len = load_well2bc(barcodesfn, bc2well, bclist);
        
        // Load cell barcodes, if provided
        set<unsigned long> cell_barcodes;
        bool has_cell_barcodes;
        if (cell_barcodesfn != ""){
            has_cell_barcodes = true;
            parse_barcode_file(cell_barcodesfn, cell_barcodes);
        }    

        // Load MULTIseq well -> unique identifier mapping
        map<string, string> well2id;
        if (mapfn == ""){
            fprintf(stderr, "WARNING: no MULTIseq well to identifier mapping provided.\n");
            fprintf(stderr, "Using well names as identifiers.\n");    
            for (map<unsigned long, string>::iterator bw = bc2well.begin(); 
                bw != bc2well.end(); ++bw){
                if (well2id.count(bw->second) == 0){
                    well2id.insert(make_pair(bw->second, bw->second));
                }
            }
        }
        else{
            load_well_mapping(mapfn, well2id);
            // Make sure wells specified here have associated valid barcodes.
            set<string> wellnames;
            for (map<unsigned long, string>::iterator bw = bc2well.begin(); bw != 
                bc2well.end(); ++bw){
                wellnames.insert(bw->second);
            }
            for (map<string, string>::iterator wi = well2id.begin(); wi != well2id.end();
                ++wi){
                if (wellnames.find(wi->first) == wellnames.end()){
                    fprintf(stderr, "ERROR: well -> ID mapping file contains a well named %s\n",
                        wi->first.c_str());
                    fprintf(stderr, "This well ID is not present in the MULTIseq barcode -> \
well mappings\n");
                    exit(1);
                }
            }
            // If the user has omitted one or more wells that turn up often in the data, we will
            // alert the user about this later.
        }
        n_labels = well2id.size();

        // Initiate MULTIseq barcode whitelist
        seq_fuzzy_match msmatch(bclist, mismatches, true, true);
        
        vector<unsigned long> ms_bc;
        bc cur_bc;
        for (int i = 0; i < bclist.size(); ++i){
            str2bc(bclist[i].c_str(), cur_bc);
            ms_bc.push_back(cur_bc.to_ulong());
        }
        /*
        // Determine appropriate k-mer length for fuzzy matching
        // k must be < (bc_length + 1)/2
        int ms_k;
        if (ms_bc_len % 2 == 0){
            // Even barcode length.
            ms_k = ms_bc_len/2;
        }
        else{
            ms_k = (ms_bc_len-1)/2;
        }
        */
        //bc_whitelist wl_ms(bclist, ms_bc_len, ms_k);
        
        map<unsigned long, int> ms_bc2idx;
        vector<string> ms_names;
        vector<string> ms_wells;
        int ix = 0;
        for (map<unsigned long, string>::iterator bw = bc2well.begin(); 
            bw != bc2well.end(); ++bw){
            
            ms_bc2idx.insert(make_pair(bw->first, ix));
            ++ix;
            ms_wells.push_back(bw->second);
            if (well2id.count(bw->second) > 0){
                ms_names.push_back(well2id[bw->second]);
                labels.push_back(well2id[bw->second]);
            }
            else{
                ms_names.push_back("");
            }
        }
        
        // Data structure to store MULTIseq barcode counts per cell barcode
        // counts come from UMIs
        robin_hood::unordered_map<unsigned long, vector<umi_set*> > bc_ms_umis;
        
        // Set up object that will scan each pair of read files
        bc_scanner scanner;
        scanner.init_multiseq_v3(wlfn);
        
        char ms_bc_buf[ms_bc_len+1];

        for (int i = 0; i < read1fn.size(); ++i){
            
            // Initiate object to read through FASTQs and find cell barcodes.
            scanner.add_reads(read1fn[i], read2fn[i]);
            
            while (scanner.next()){
                
                // Skip reads if we have a cell barcode whitelist and the current
                // barcode isn't in it
                if (has_cell_barcodes && 
                    cell_barcodes.find(scanner.barcode) == cell_barcodes.end()){
                    continue;
                }
                
                // At this point, there's a valid cell barcode.
                // scanner.barcode_read holds R1
                // scanner.read_f holds R2
                // scanner.umi holds UMI 
                if (scanner.has_umi){
                    // MULTIseq barcodes are the first 8 bp (or however long supplied barcodes are)
                    // of R2, forward orientation
                    unsigned long ms_ul;
                    
                    strncpy(&ms_bc_buf[0], scanner.read_f, ms_bc_len);
                    ms_bc_buf[ms_bc_len] = '\0';            
                    int idx = msmatch.match(&ms_bc_buf[0]);
                    if (idx != -1){
                        ms_ul = ms_bc[idx];
                    //if (wl_ms.lookup(scanner.read_f, false, ms_ul)){
                        // Store UMI
                        umi this_umi(scanner.umi, scanner.umi_len);    
                        if (bc_ms_umis.count(scanner.barcode) == 0){
                            vector<umi_set*> v(bclist.size(), NULL);
                            bc_ms_umis.emplace(scanner.barcode, v);
                            //map<unsigned long, umi_set> m;
                            //bc_ms_umis.emplace(scanner.barcode, m);
                        }
                        umi_set* ptr = bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]];
                        if (ptr == NULL){
                            // Initialize.
                            bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]] = new umi_set(scanner.umi_len);
                            ptr = bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]];
                        }
                        ptr->add(this_umi); 
                    }
                }
            }
        }

        // Write counts to disk
        string countsfn = output_prefix + ".counts";
        dump_counts(countsfn, bc_ms_umis, ms_names, bc_ms_counts);
       
        // Check to see that the desired MULTIseq barcodes are the most common ones.
        // Warn the user if unexpected ones are more common. 
        string wellcountsfn = output_prefix + ".wells";
        check_missing_ms_wells(bc_ms_umis, ms_wells, well2id, wellcountsfn);        

        // Free stuff
        for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x = 
            bc_ms_umis.begin(); x != bc_ms_umis.end(); ){
            for (int i = 0; i < x->second.size(); ++i){
                if (x->second[i] != NULL){
                    delete x->second[i];
                }
            }
            bc_ms_umis.erase(x++);
        }
    }

    // Now we can do the actual demultiplexing
    robin_hood::unordered_map<unsigned long, set<int> > assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    robin_hood::unordered_map<unsigned long, double> cell_bg; 
    assign_ids2(bc_ms_counts, assn, assn_llr, cell_bg, labels, filt);

    string assnfilename = output_prefix + ".assignments";
    write_assignments(assnfilename, assn, labels, assn_llr, sep, batch_id, sgrna);
    
    if (sgrna){
        string tablefilename = output_prefix + ".table";
        write_assignments_table(tablefilename, assn, labels, batch_id);
    }

    string bgfilename = output_prefix + ".bg";
    write_bg(bgfilename, cell_bg, batch_id);

    return 0;  
}
