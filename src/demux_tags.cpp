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
    fprintf(stderr, "demux_multiseq [OPTIONS]\n");
    fprintf(stderr, "Given reads containing MULTIseq data, counts the occurrences of each\n");
    fprintf(stderr, "  MULTIseq barcode per cell barcode and assigns each cell barcode\n");
    fprintf(stderr, "  an identity based on MULTIseq barcode meaning.\n");
    fprintf(stderr, "A single run will produce a .counts file and an .assignments file. If\n");
    fprintf(stderr, "  you run again with the same output_prefix, the counts from the prior\n");
    fprintf(stderr, "  run will be loaded.\n");
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
    
    fprintf(stderr, "   --barcodes -b A path to a file listing MULTIseq well IDs and corresponding\n");
    fprintf(stderr, "       barcodes, tab separated. If you do not provide one, the default file (in\n");
    fprintf(stderr, "       the data directory) will be used.\n");
    fprintf(stderr, "   --mismatches -M The number of allowable mismatches to a MULTIseq barcode to accept\n");
    fprintf(stderr, "       it. Default = -1 (take best match overall, if there are no ties)\n");
    fprintf(stderr, "   --doublet_rate -D What is the prior expected doublet rate?\n");
    fprintf(stderr, "       (OPTIONAL; default = 0.1). Must be a decimal between 0 and 1,\n");
    fprintf(stderr, "       exclusive.\n");
    fprintf(stderr, "   --output_unsassigned -u Print entries for unassigned barcodes in assignments file\n");
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

int assign_cell(mixtureModel& mod,
    vector<double>& obsrow,
    vector<double>& obs2row,
    int ndim,
    pair<double, double>& a_b,
    map<int, pair<int, int> >& ddoub,
    map<int, int>& ddim,
    bool add_bg,
    double& ll,
    double& llr){
    
    double tot = obsrow[obsrow.size()-1];
    double pbg = a_b.first * log(tot) + a_b.second;
    pbg = exp(pbg)/(exp(pbg) + 1);

    vector<pair<double, int> > lls;
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
    vector<string>& labels){
    
    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
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
        }
        tots.push_back(tot);
        obs2.push_back(row);
        row.push_back(tot);
        obs.push_back(row);
        bcs.push_back(x->first);
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
    mod.fit(obs);
    
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
    double sumsq = 0.0;
    for (int i = 0 ; i < totcounts.size(); ++i){
        double pred = a_b.first * totcounts[i] + a_b.second;
        sumsq += pow(pred-pbgs[i], 2);
    } 
    fprintf(stderr, "=== Total/ambient tag count relationship ===\n");
    fprintf(stderr, "b = (%.4f*c^%.4f)/(%.4f*c^%.4f + 1)\n",
        exp(a_b.second), a_b.first, exp(a_b.second), a_b.first);
    fprintf(stderr, "  b = percent of tag counts from ambient\n");
    fprintf(stderr, "  c = total tag counts per cell\n");    
    fprintf(stderr, "SSE = %f\n", sumsq);
    
    vector<vector<double> > filt_dat;
    vector<unsigned long> filt_dat_bc;

    for (int i = 0; i < obs.size(); ++i){
        double llr = 0.0;
        double ll = 0.0;
        int assn_this = assign_cell(mod, obs[i], obs2[i], labels.size(),
            a_b, ddoub, ddim, add_bg, ll, llr);
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
        if (*a == 1 && mod_filt.responsibility_matrix[idx][*a] > 0.75){
            cfilt++;
            assn.erase(filt_dat_bc[idx]);
            assn_llr.erase(filt_dat_bc[idx]);
        }
        ++idx;
    }
    
    fprintf(stderr, "%d filt\n", cfilt);
}

/**
 * Spill cell -> identity assignments to disk (in .assignments file)
 */
void write_assignments(string filename, 
    robin_hood::unordered_map<unsigned long, string>& assn,
    robin_hood::unordered_map<unsigned long, bool>& assn_doublet,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    FILE* f_out = fopen(filename.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn.begin(); 
        a != assn.end(); ++a){
        char type = 'S';
        if (assn_doublet.count(a->first) > 0 && assn_doublet[a->first]){
            type = 'D';
        }
        string name = a->second;
        string bc_str = bc2str(a->first);
        fprintf(f_out, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type, assn_llr[a->first]); 
    }
    fclose(f_out);
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
    /* 
    ifstream inf("mvr.samp.txt");
    double n1;
    double n2;
    vector<vector<double> > obstest;
    while (inf >> n1 >> n2){
        vector<double> row{ n1, n2 };
        obstest.push_back(row);
    } 
    mixtureDist d1("2dgauss", vector<double>{ 200, 0.1, 50, 0.25, 0 });
    mixtureDist d2("2dgauss", vector<double>{ 600, 0.5, 100, 0.1, 0 });
    vector<mixtureDist> dists;
    dists.push_back(d1);
    dists.push_back(d2);
    mixtureModel m(dists);
    m.fit(obstest);
    m.print();
    exit(0);
    */

    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"mapping", required_argument, 0, 'm'},
       {"whitelist", required_argument, 0, 'w'},
       {"barcodes", required_argument, 0, 'b'},
       {"cell_barcodes", required_argument, 0, 'B'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"mismatches", required_argument, 0, 'M'},
       {"output_unassigned", no_argument, 0, 'u'},
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
    bool output_unassigned = false;

    fprintf(stderr, "%s\n", STRINGIZE(PROJ_ROOT));
    string pr = STRINGIZE(PROJ_ROOT);
    if (pr[pr.length()-1] != '/'){
        pr += "/";
    }
    string barcodesfn = pr + "data/multiseq_indices.txt";
    double doublet_rate = 0.1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:1:2:m:M:w:b:B:D:uh", 
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
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'u':
                output_unassigned = true;
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
    robin_hood::unordered_map<unsigned long, string> assn;
    robin_hood::unordered_map<unsigned long, bool> assn_doublet;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    assign_ids(bc_ms_counts, assn, assn_doublet, assn_llr, labels);

    string assnfilename = output_prefix + ".assignments";
    write_assignments(assnfilename, assn, assn_doublet, assn_llr);

    return 0;  
}
