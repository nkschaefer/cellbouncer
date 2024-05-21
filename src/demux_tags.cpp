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

double ll_multinom_mix(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    vector<double> params_untx;
    double grptot = 0.0;
    for (int i = 0; i < params.size()-1; ++i){
        double xt = expit(params[i]);
        grptot += xt;
        params_untx.push_back(xt);
    }
    for (int i = 0; i < params_untx.size(); ++i){
        params_untx[i] /= grptot;
    }

    double prop_noise = params[params.size()-1];
    int dim_idx = data_i.at("dim_idx");
    vector<double> params_this;
    vector<double> obs_this;
    char keybuf[100];
    string keystr;
    for (int i = 0; i < params.size()-1; ++i){
        if (i == dim_idx){
            params_this.push_back((1.0-prop_noise)*1.0 + prop_noise*params_untx[i]);
        }
        else{
            params_this.push_back((1.0-prop_noise)*0.0 + prop_noise*params_untx[i]);
        }
        sprintf(&keybuf[0], "n_%d", i);
        keystr = keybuf;
        int x = data_i.at(keystr);
        obs_this.push_back((double)x);
    }
    return dmultinom(obs_this, params_this) / log2(exp(1));
}

void dll_multinom_mix(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    vector<double> params_untx;
    double grptot = 0.0;
    for (int i = 0; i < params.size()-1; ++i){
        double xt = expit(params[i]);
        grptot += xt;
        params_untx.push_back(xt);
    }
    for (int i = 0; i < params_untx.size(); ++i){
        params_untx[i] /= grptot;
    }

    double prop_noise = params[params.size()-1];
    int prop_noise_idx = params.size()-1;
    int dim_idx = data_i.at("dim_idx");
    
    // First compute dLL / dp_i for each i
    vector<double> dLL_dp;
    vector<double> p_i_all;

    char keybuf[100];
    string keystr;
    for (int i = 0; i < params.size()-1; ++i){
        
        // Retrieve count
        sprintf(&keybuf[0], "n_%d", i);
        keystr = keybuf;
        int x = data_i.at(keystr);

        // Compute p_i
        double p_i;
        if (i == dim_idx){
            p_i = (1.0 - prop_noise)*1.0 + prop_noise * params_untx[i];
        }
        else{
            p_i = prop_noise * params_untx[i];
        }
        p_i_all.push_back(p_i);
        
        double dLL_dpi = (double)x / p_i;
        
        // Incorporate constraint
        double e_negx1 = exp(-params[i]);
        double e_negx1_p1_2 = pow(e_negx1 + 1, 2);
        double e_negx1_p1_3 = e_negx1_p1_2 * (e_negx1 + 1);
        double chain_rule = e_negx1 / (e_negx1_p1_2 * grptot) - e_negx1 / (e_negx1_p1_3 * grptot * grptot);
        dLL_dpi *= chain_rule;
        dLL_dp.push_back(dLL_dpi);
    }
    
    for (int i = 0 ; i < dLL_dp.size(); ++i){
        results[i] += dLL_dp[i] * prop_noise;
        if (i == dim_idx){
            // on target
            results[results.size()-1] += dLL_dp[i] * (params_untx[i] - 1.0);
        }
        else{
            // off target
            results[results.size()-1] += dLL_dp[i] * params_untx[i];
        }
    }
}

map<int, pair<int, int> > dist2doublet;
map<int, int> dist2dim;
vector<vector<double> > obs;

double callback_func_eq(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    // NOTE: having to un-transform all variables on each function call is obviously
    // computationally wasteful. But since there shouldn't be all that many mixture
    // model components, hopefully this won't be a huge problem.

    vector<double> params_untx;
    double grptot = 0.0;
    for (int i = 0; i < params.size()-1; ++i){
        double xt = expit(params[i]);
        grptot += xt;
        params_untx.push_back(xt);
    }
    for (int i = 0; i < params_untx.size(); ++i){
        params_untx[i] /= grptot;
    }
    
    double prop_noise = params[params.size()-1];

    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    int param_idx = data_i.at("param_idx");

    if (idx2 != -1){
        if (param_idx == idx1){
            return prop_noise * params_untx[param_idx] + (1.0 - prop_noise)*0.5;
        }
        else if (param_idx == idx2){
            return prop_noise * params_untx[param_idx] + (1.0 - prop_noise)*0.5;
        }
        else{
            return prop_noise * params_untx[param_idx];
        }
    }
    else{
        if (param_idx == idx1){
            return prop_noise * params_untx[param_idx] + (1.0 - prop_noise);
        }
        else{
            return prop_noise * params_untx[param_idx];
        }    
    }
}

/*
void callback_func_eq_d(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    vector<double> params_untx;
    double grptot = 0.0;
    for (int i = 0; i < params.size()-1; ++i){
        double xt = expit(params[i]);
        grptot += xt;
        params_untx.push_back(xt);
    }
    for (int i = 0; i < params_untx.size(); ++i){
        params_untx[i] /= grptot;
    } 

    double prop_noise = params[params.size()-1];
    int prop_noise_idx = params.size()-1;

    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    int param_idx = data_i.at("param_idx");
    
    if (idx2 != -1){
        if (param_idx == idx1 || param_idx == idx2){
            results[param_idx] += prop_noise;
            results[prop_noise_idx] += (params_untx[param_idx] - 0.5);
        }
        else{
            results[param_idx] += prop_noise;
            results[prop_noise_idx] += params_untx[param_idx];
        }
    }
    else{
        if (param_idx == idx1){
            results[param_idx] += prop_noise;       
            results[prop_noise_idx] += (params_untx[param_idx] - 1.0);
        }
        else{
            results[param_idx] += prop_noise;
            results[prop_noise_idx] += params_untx[param_idx];
        }
    }
    
    double e_negx1 = exp(-params[param_idx]);
    double e_negx1_p1_2 = pow(e_negx1 + 1, 2);
    double e_negx1_p1_3 = e_negx1_p1_2 * (e_negx1 + 1);
    double chain_rule = e_negx1 / (e_negx1_p1_2 * grptot) - e_negx1 / (e_negx1_p1_3 * grptot * grptot);
    results[param_idx] *= chain_rule;
}

void callback_func(const vector<mixtureDist*>& dists,
    const vector<double>& weights,
    vector<double>& shared_params){
    
    double perc_bg = shared_params[shared_params.size()-1];
    
    optimML::multivar_sys_solver solver(shared_params);

    for (int i = 0; i < dists.size(); ++i){
        if (weights[i] > 0.0){
            int idx1_this;
            int idx2_this;
            if (dist2doublet.count(i) > 0){
                idx1_this = dist2dim[dist2doublet[i].first];
                idx2_this = dist2dim[dist2doublet[i].second];
            }
            else{
                idx1_this = dist2dim[i];
                idx2_this = -1;
            }

            for (int j = 0; j < dists[i]->params[0].size(); ++j){
                solver.add_equation(callback_func_eq, callback_func_eq_d, dists[i]->params[0][j]);
                solver.add_data("idx1", idx1_this);
                solver.add_data("idx2", idx2_this);
                solver.add_data("param_idx", j);
            }
        }
    } 
    
    solver.solve();
    
    shared_params[shared_params.size()-1] = solver.results[shared_params.size()-1];
    fprintf(stderr, "Contam frac %f\n", shared_params[shared_params.size()-1]);

    double grptot = 0.0;
    for (int i = 0; i < shared_params.size()-1; ++i){
        // Keep logit transformed
        shared_params[i] = solver.results[i]; 
        grptot += expit(solver.results[i]);
    }
    
    fprintf(stderr, "RESULTS\n");
    for (int i = 0; i < solver.results.size(); ++i){
        double val = solver.results[i];
        if (i < solver.results.size()-1){
            val = expit(solver.results[i]) / grptot;
        }
        fprintf(stderr, "%d) %f\n", i, val);
    }
    exit(1);
    double frac_noise = shared_params[shared_params.size()-1];

    // Update params in mixture model
    for (int i = 0; i < dists.size(); ++i){
        if (dist2doublet.count(i) > 0){
            int dim1 = dist2dim[dist2doublet[i].first];
            int dim2 = dist2dim[dist2doublet[i].second];
            for (int j = 0; j < dists[i]->params[0].size(); ++j){
                if (j == dim1 || j == dim2){
                    dists[i]->params[0][j] = (1.0 - frac_noise)*0.5 + frac_noise*expit(solver.results[j])/grptot;
                }
                else{
                    dists[i]->params[0][j] = frac_noise*expit(solver.results[j])/grptot;
                }
            }
        }
        else{
            int dim = dist2dim[i];
            for (int j = 0; j < dists[i]->params[0].size(); ++j){
                if (j == dim){
                    dists[i]->params[0][j] = (1.0 - frac_noise) + frac_noise*shared_params[j];
                }
                else{
                    dists[i]->params[0][j] = frac_noise * shared_params[j];
                }
            }
        }
    } 
    
}
*/

double tags_ll(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int n_dim = data_i.at("ndim");
    int target1 = data_i.at("target1");
    int target2 = data_i.at("target2");

    double frac_bg = params[0];
    
    static vector<double> x;
    static vector<double> p;
    if (x.size() == 0){
        for (int i = 0; i < n_dim; ++i){
            x.push_back(0.0);
            p.push_back(0.0);
        }
    }
    static char buf[50];
    static string bufstr;
    for (int i = 0; i < n_dim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        int x_this = data_i.at(bufstr);
        x[i] = (double)x_this;
        double p_this = params[i+1];
        
        if (i == target1){
            if (target2 != -1){
                p_this = frac_bg * p_this + (1.0 - frac_bg)*0.5;
            }
            else{
                p_this = frac_bg * p_this + (1.0 - frac_bg)*1.0;
            }
        }
        else if (i == target2){
            p_this = frac_bg * p_this + (1.0 - frac_bg)*0.5;
        }
        p[i] = p_this;
        
        if (p[i] <= 0 || p[i] >= 1){
            fprintf(stderr, "WTF %d %d %d\n", i, target1, target2);
            fprintf(stderr, "ndim %d\n", n_dim);
            fprintf(stderr, "frac_bg %f p0[i] %f\n", frac_bg, params[i+1]);
        }
    }
    return dmultinom(x, p) * log(2);
}

void tags_dll(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
 
    int n_dim = data_i.at("ndim");
    int target1 = data_i.at("target1");
    int target2 = data_i.at("target2");

    double frac_bg = params[0];
   
    static char buf[50];
    static string bufstr;
    
    double tot = 0.0;
    for (int i = 0; i < n_dim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        double x_this = (double)data_i.at(bufstr);
        tot += x_this;
    }

    for (int i = 0; i < n_dim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        double x_this = (double)data_i.at(bufstr);
        double p_this = params[i+1];
        
        double dLL_dp = 0.0;
        double dp_dp0 = 0.0;
        double dp_df = 0.0;

        if (i == target1){
            if (target2 != -1){
                p_this = frac_bg * p_this + (1.0 - frac_bg)*0.5;
                dLL_dp = x_this / p_this - tot;
                dp_df = params[i+1] - 1.0;
                dp_dp0 = frac_bg;
            }
            else{
                p_this = frac_bg * p_this + (1.0 - frac_bg)*1.0;
                dLL_dp = x_this / p_this - tot;
                dp_dp0 = frac_bg;
                dp_df = params[i+1] - 0.5;
            }
        }
        else if (i == target2){
            p_this = frac_bg * p_this + (1.0 - frac_bg)*0.5;
            dLL_dp = x_this / p_this - tot;
            dp_dp0 = frac_bg;
            dp_df = params[i+1] - 0.5;
        }
        else{
            dLL_dp = x_this / p_this - tot;
            dp_dp0 = 1.0;
            dp_df = 0.0;
        }
        results[0] += (dLL_dp*dp_df);
        results[i+1] += (dLL_dp*dp_dp0);    
    }
}

double callback_eq_target(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    return (params[idx]*frac_bg + (1.0 - frac_bg)*1.0); 
}

void callback_eq_target_d(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    results[idx] += frac_bg;
    results[0] += (params[idx] - 1.0);
}

double callback_eq_doublet(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    return (params[idx]*frac_bg + (1.0 - frac_bg)*0.5);
}

void callback_eq_doublet_d(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    results[idx] += frac_bg;
    results[0] += (params[idx] - 0.5);
}

double callback_eq_offtarget(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    return (params[idx]*frac_bg);
}

void callback_eq_offtarget_d(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    double frac_bg = params[0];
    results[idx] += frac_bg;
    results[0] += params[idx];
}

double callback_eq_noise(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    return params[idx];
}

void callback_eq_noise_d(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    results[0] = 0.0;
    results[idx] = 1.0;
}

void callback_func2(mixtureModel& mod, vector<double>& params){
    
    optimML::multivar_sys_solver solver(vector<double>{ params[0] });
    solver.constrain_01(0);
    vector<double> pg;
    for (int i = 1; i < params.size(); ++i){
        pg.push_back(params[i]);
    }
    solver.add_param_grp(pg);
    for (int i = 0; i < mod.dists.size(); ++i){
        for (int dim_idx = 0; dim_idx < mod.dists[i].params[0].size(); ++dim_idx){
            if (dist2doublet.count(i) > 0){
                if (dim_idx == dist2dim[dist2doublet[i].first] || 
                    dim_idx == dist2dim[dist2doublet[i].second]){
                    solver.add_equation(callback_eq_doublet, callback_eq_doublet_d,
                        mod.dists[i].params[0][dim_idx], mod.weights[i]);
                    solver.add_data("idx", dim_idx);
                }
                else{
                    solver.add_equation(callback_eq_offtarget, callback_eq_offtarget_d,
                        mod.dists[i].params[0][dim_idx], mod.weights[i]);
                    solver.add_data("idx", dim_idx);
                }    
            }
            else{
                if (dist2dim.count(i) > 0){
                    if (dist2dim[i] == dim_idx){
                        solver.add_equation(callback_eq_target, callback_eq_target_d,
                            mod.dists[i].params[0][dim_idx], mod.weights[i]);
                        solver.add_data("idx", dim_idx);
                    }
                    else{
                        solver.add_equation(callback_eq_offtarget, callback_eq_offtarget_d,
                            mod.dists[i].params[0][dim_idx], mod.weights[i]);
                        solver.add_data("idx", dim_idx);
                    }
                }
                else{
                    solver.add_equation(callback_eq_noise, callback_eq_noise_d,
                        mod.dists[i].params[0][dim_idx], mod.weights[i]);
                    solver.add_data("idx", dim_idx);
                }
            }
        }
    }
    solver.solve();
    double frac_bg = solver.results[0];

    for (int i = 0; i < mod.dists.size(); ++i){
        for (int dim_idx = 0; dim_idx < mod.dists[i].params[0].size(); ++dim_idx){
            if (dist2doublet.count(i) > 0){
                if (dist2dim[dist2doublet[i].first] == dim_idx ||
                    dist2dim[dist2doublet[i].second] == dim_idx){
                    mod.dists[i].params[0][dim_idx] = frac_bg * solver.results[dim_idx+1] + 
                        (1.0-frac_bg)*0.5;                
                }
                else{
                    mod.dists[i].params[0][dim_idx] = frac_bg * solver.results[dim_idx+1];
                }
            }
            else if (dist2dim.count(i) > 0){
                if (dist2dim[i] == dim_idx){
                    mod.dists[i].params[0][dim_idx] = frac_bg * solver.results[dim_idx+1] + 
                        (1.0-frac_bg)*1.0;
                }
                else{
                    mod.dists[i].params[0][dim_idx] = frac_bg * solver.results[dim_idx+1];
                }
            }
            else{
                mod.dists[i].params[0][dim_idx] = solver.results[dim_idx+1];
            }
        }
    }
}

void callback_func(mixtureModel& mod, vector<double>& params){
    
    // params[0] is fraction background

    optimML::multivar_ml_solver solver(vector<double>{ params[0] }, tags_ll, tags_dll);
     
    vector<double> weights;
    vector<vector<int> > ns;
    vector<int> target1;
    vector<int> target2; 
    
    int ndim = obs[0].size();
    vector<double> paramgrp;
    for (int i = 0; i < ndim; ++i){
        vector<int> v;
        ns.push_back(v);
        paramgrp.push_back(params[i+1]);
    }
    solver.add_param_grp(paramgrp);

    for (int i = 0; i < obs.size(); ++i){
        for (int j = 0; j < mod.dists.size(); ++j){
            double weight = mod.responsibility_matrix[i][j];
            weights.push_back(weight);
            for (int k = 0; k < ndim; ++k){
                ns[k].push_back(obs[i][k]);
            }
            if (dist2doublet.count(j) > 0){
                target1.push_back(dist2dim[dist2doublet[j].first]);
                target2.push_back(dist2dim[dist2doublet[j].second]);
            }
            else{
                target1.push_back(dist2dim[j]);
                target2.push_back(-1);
            }
        }
    }
    
    solver.add_data("target1", target1);
    solver.add_data("target2", target2);
    char buf[50];
    for (int j = 0; j < ndim; ++j){
        sprintf(&buf[0], "n_%d", j);
        string bufstr = buf;
        solver.add_data(bufstr, ns[j]);
    }
    solver.add_data_fixed("ndim", ndim);
    solver.add_weights(weights);
    solver.constrain_01(0);
    fprintf(stderr, "solving\n");
    solver.solve();
    fprintf(stderr, "done\n");
    for (int i = 0; i < solver.results.size(); ++i){
        fprintf(stderr, "RESULT %d) %f\n", i, solver.results[i]);
    }
    exit(0);
}

void assign_ids_new2(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    vector<string> labels,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned){
    
    int n_labels = labels.size();
    int ndim = n_labels;

    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }

    vector<unsigned long> bcs;

    vector<mixtureDist> dists;
    
    double off_target = 0.1;
    double covhi = 100;
    double covlow = 1;
    
    double perc_noise = 0.1;

    vector<double> background_dist;
    for (int i = 0; i < n_labels; ++i){
        background_dist.push_back(1.0 / (double)n_labels);
    }
    
    mixtureDist dist_bg("multinomial", background_dist);
    dist_bg.name = "Unassigned";
    dist_bg.set_num_inputs(0, n_labels);
    dists.push_back(dist_bg);

    double perc_background = 0.1;
    
    for (int i = 0; i < n_labels; ++i){
        int di = dists.size();
        dist2dim.insert(make_pair(di, i));
        vector<double> params;
        for (int j = 0; j < n_labels; ++j){
            if (j == i){
                //params.push_back(1.0 - off_target);
                params.push_back(perc_background * background_dist[j] + (1.0-perc_background)*1.0);
            }
            else{
                params.push_back(perc_background * background_dist[j]);
                //params.push_back(off_target / (double)(n_labels-1));
            }
        }
        /*
        for (int j = 0; j < n_labels; ++j){
            params[j] = (1.0-perc_background)*params[j] + perc_background*background_dist[j];
        }
        */
        mixtureDist dist("multinomial", params);
        dist.name = labels[i];
        dist.set_num_inputs(0, n_labels);
        dists.push_back(dist);
    } 

    // Handle doublets
    for (int i = 0; i < n_labels-1; ++i){
        for (int j = i + 1; j < n_labels; ++j){
            
            int di = dists.size();
            dist2doublet.insert(make_pair(di, make_pair(i, j)));

            string name;
            if (dists[i].name < dists[j].name){
                name = dists[i].name + "+" + dists[j].name;
            }
            else{
                name = dists[j].name + "+" + dists[i].name;
            }
            vector<double> params;
            for (int k = 0; k < n_labels; ++k){
                params.push_back(0.5 * dists[i].params[0][k] + 0.5 * dists[j].params[0][k]);
            }
            mixtureDist dist("multinomial", params);
            dist.name = name;
            dist.set_num_inputs(0, n_labels);
            dists.push_back(dist);
        } 
    }
    
    mixtureModel mod(dists);
    vector<double> shared_params;
    shared_params.push_back(perc_background);
    for (int i = 0; i < background_dist.size(); ++i){
        shared_params.push_back(background_dist[i]);
    }
    mod.set_callback(callback_func2, shared_params);
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        int tot = 0;
        vector<double> row;
        for (vector<string>::iterator l = labels.begin(); l != labels.end(); ++l){
            row.push_back(0.0);            
        }
        for (int i = 0; i < x->second.size(); ++i){
            row[label2idx[x->second[i].second]] += (double)-x->second[i].first;
            tot += -x->second[i].first;
        }
        if (tot > 0){
            obs.push_back(row);
            bcs.push_back(x->first);
        }
    }
    mod.fit(obs);
    
    for (int i = 0; i < obs.size(); ++i){
        vector<pair<double, int> > lls;
        for (int j = 0; j < mod.dists.size(); ++j){
            double ll = mod.dists[j].loglik(obs[i]);
            lls.push_back(make_pair(-ll, j)); 
        }
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;
        int j = lls[0].second;
        if (llr > 0){
            unsigned long bc = bcs[i];
            if (dist2doublet.count(j) > 0){
                assn1.emplace(bc, mod.dists[dist2doublet[j].first].name);
                assn2.emplace(bc, mod.dists[dist2doublet[j].second].name);
                assn_llr.emplace(bc, llr);
            }
            else if (dist2dim.count(j) > 0){
                assn1.emplace(bc, mod.dists[j].name);
                assn_llr.emplace(bc, llr);
            }
        }
    }

}

void assign_ids3(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    int n_labels,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned){

    vector<vector<double> > obs;
    vector<unsigned long> bcs;

    // Define 3 distributions:
    // Single-label, high count
    // Double-label, high count
    // Noise (3+ labels), low count
    vector<mixtureDist> dists;
    
    double off_target = 0.1;
    vector<double> params_single = { 1.0-off_target, off_target/2, off_target/2 };
    vector<double> params_double = { (1.0-off_target)/2, (1.0-off_target)/2, off_target };
    vector<double> params_noise = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
    
    /* 
    vector<string> names = {"multinomial", "poisson"};
    mixtureDist singlet(names, vector<vector<double> >{ params_single, vector<double>{ 100 } });
    singlet.name = "singlet";
    mixtureDist doublet(names, vector<vector<double> >{ params_double, vector<double>{ 100 } });
    doublet.name = "doublet";
    mixtureDist noise(names, vector<vector<double> >{ params_noise, vector<double>{ 1 } });
    noise.name = "noise";
    singlet.set_num_inputs(0, 3);
    doublet.set_num_inputs(0, 3);
    noise.set_num_inputs(0, 3); 
    dists.push_back(singlet);
    dists.push_back(doublet);
    dists.push_back(noise);
    */
    mixtureDist singlet("multinomial", params_single);
    mixtureDist doublet("multinomial", params_double);
    mixtureDist noise("multinomial", params_noise);
    singlet.name = "singlet";
    doublet.name = "doublet";
    noise.name = "noise";
    singlet.set_num_inputs(0, 3);
    doublet.set_num_inputs(0, 3);
    noise.set_num_inputs(0, 3);
    dists.push_back(singlet);
    dists.push_back(doublet);
    dists.push_back(noise);

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        sort(x->second.begin(), x->second.end());
        if (x->second[0].first == 0){
            // No counts
            continue;
        }    

        int tot = 0;
        int past2 = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
            if (i > 1){
                past2 += -x->second[i].first;
            }
        }
        vector<double> row = { (double)-x->second[0].first, (double)-x->second[1].first, (double)past2 };
        fprintf(stderr, "ROW %d %d %d\n", -x->second[0].first, -x->second[1].first, past2);
        obs.push_back(row);
        bcs.push_back(x->first);
    }
    
    mixtureModel mod(dists);
    mod.fit(obs);
    mod.print();
    exit(0);

}

double median(vector<double> vals){
    sort(vals.begin(), vals.end());
    if (vals.size() == 0){
        return 0;
    }
    else if (vals.size() == 1){
        return vals[0];
    }
    if (vals.size() % 2 == 0){
        double v1 = vals[vals.size()/2];
        double v2 = vals[vals.size()/2-1];
        return (0.5*v1 + 0.5*v2);
    }
    else{
        return vals[(vals.size()-1)/2];
    }
}

double median(vector<int>& vals){
    sort(vals.begin(), vals.end());
    if (vals.size() == 0){
        return 0;
    }
    else if (vals.size() == 1){
        return vals[0];
    }
    if (vals.size() % 2 == 0){
        double v1 = (double)vals[vals.size()/2];
        double v2 = (double)vals[vals.size()/2-1];
        return (0.5*v1 + 0.5*v2);
    }
    else{
        return vals[(vals.size()-1)/2];
    }
}

map<int, int> ddim;
map<int, pair<int, int> > ddoub;

double doub_prop(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int i1 = data_i.at("idx1");
    int i2 = data_i.at("idx2");
    fprintf(stderr, "rhs %f\n", params[i1]/(params[i1] + params[i2]));
    return params[i1] / (params[i1] + params[i2]);
}

void ddoub_prop(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int i1 = data_i.at("idx1");
    int i2 = data_i.at("idx2");
    results[i1] += params[i2] / (pow(params[i1] + params[i2], 2));
    results[i2] -= params[i1] / (pow(params[i1] + params[i2], 2));
    fprintf(stderr, "deriv %d) %f\n", i1, results[i1]);
    fprintf(stderr, "deriv %d) %f\n", i2, results[i2]);
}

double ll_mix_bg_binom(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");

    double bp = data_d.at("bg");
    
    return dbinom(n, k, bp*param);
   
}

double dll_mix_bg_binom(double param,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double bp = data_d.at("bg");

    return (k - bp*n*param)/(param - bp*param*param);
}

double ll_mix_new(double param,
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

double dll_mix_new(double param,
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

double ll_mix(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    int idx1 = data_i.at("idx1");
    int idx2 = -1;
    if (data_i.count("idx2") > 0){
        idx2 = data_i.at("idx2");
    }
    int ndim = data_i.at("ndim");
    vector<double> x;
    vector<double> p;
    
    /* 
    fprintf(stderr, "idx1 %d idx2 %d ndim %d\n", idx1, idx2, ndim);
    fprintf(stderr, "params:");
    for (int i = 0; i < params.size(); ++i){
        fprintf(stderr, " %f", params[i]);
    } 
    fprintf(stderr, "\n");
    */

    char buf[50];
    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        x.push_back(data_d.at(bufstr));
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgfrac = data_d.at(bufstr);
        if (idx2 != -1){
            if (i == idx1){
                //p.push_back(params[0] * 1.0 + (1.0 - params[0] - params[1])*bgfrac);
                p.push_back(params[0] * 1.0 + params[2]*bgfrac);
            }
            else if (i == idx2){
                //p.push_back(params[1] * 1.0 + (1.0 - params[0] - params[1])*bgfrac);
                p.push_back(params[1] * 1.0 + params[2]*bgfrac);
            }
            else{
                //p.push_back((1.0 - params[0] - params[1])*bgfrac);
                p.push_back(params[2]*bgfrac);
            }
        }
        else{
            if (i == idx1){
                p.push_back(params[0] * 1.0 + (1.0 - params[0])*bgfrac);
            }
            else{
                p.push_back((1.0 - params[0]) * bgfrac);
            }
        }
    }

    return dmultinom(x, p)/log2(exp(1));
}

void dll_mix(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx1 = data_i.at("idx1");
    int idx2 = -1;
    if (data_i.count("idx2") > 0){
        idx2 = data_i.at("idx2");
    }
    int ndim = data_i.at("ndim");

    char buf[50];

    vector<double> x_all;
    vector<double> p_bg;

    double xsum_offtarget = 0.0;
    double psum_offtarget = 0.0;

    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        double x_i = data_d.at(bufstr);
        sprintf(&buf[0], "bg_%d", i);
        bufstr = buf;
        double bgfrac = data_d.at(bufstr);
        x_all.push_back(x_i);
        p_bg.push_back(bgfrac);
        if (idx2 != -1){
            if (i != idx1 && i != idx2){
                xsum_offtarget += x_i;
                psum_offtarget += bgfrac;
            }
        }
    }

    double lambda;
    if (idx2 != -1){
        lambda = xsum_offtarget / psum_offtarget / params[2];
    }
    
    for (int i = 0; i < ndim; ++i){
        double bgfrac = p_bg[i];
        double x = x_all[i];
        if (idx2 != -1){
            if (i == idx1){
                results[0] += x/(params[0] + params[2]*bgfrac);
                results[2] += (bgfrac*x)/(params[0] + params[2]*bgfrac);
                //results[0] += ((bgfrac - 1)*x)/(params[0]*(bgfrac-1) + (params[1] - 1)*bgfrac);
                //results[1] += (bgfrac*x)/(params[0]*(bgfrac-1) + (params[1]-1)*bgfrac);
            }
            else if (i == idx2){
                results[1] += x/(params[1] + params[2]*bgfrac);
                results[2] += (bgfrac*x)/(params[1] + params[2]*bgfrac);
                //results[0] += (x*bgfrac)/((params[0]-1)*bgfrac + params[1]*(bgfrac-1));
                //results[1] += ((bgfrac-1)*x)/((params[0]-1)*bgfrac + params[1]*(bgfrac-1));
            }
            else{
                results[2] += x/bgfrac;
                //results[0] += x/(params[0] + params[1] - 1);
                //results[1] += x/(params[0] + params[1] - 1);
            }
        }
        else{
            if (i == idx1){
                results[0] += (bgfrac*x - x)/(params[0]*bgfrac - params[0] - bgfrac);
            }
            else{
                results[0] += x/(params[0] - 1.0);
            }
        }
    }
    if (idx2 != -1){
        results[0] -= lambda;
        results[1] -= lambda;
        results[2] -= lambda;
    }
}

double bg_props(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int idx = data_i.at("idx");
    
    return params[idx];
}

void dbg_props(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    int idx = data_i.at("idx");
    results[idx] += 1.0;
}

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

void cb_new_take2(mixtureModel& mod, vector<double>& params){
    
    vector<int> idx1;
    vector<int> idx2;
    vector<double> frac;
     
    /* 
    optimML::multivar_sys_solver solver(params);
    for (int i = 0; i < params.size(); ++i){
        solver.constrain_01(i);
    }
    */
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

    //optimML::multivar_sys_solver solver( vector<double>{ });
    //solver.add_param_grp(params);
    
    for (map<int, pair<int, int> >::iterator d = ddoub.begin(); d != ddoub.end(); ++d){
        if (mod.weights[d->second.first] == 0.0 || mod.weights[d->second.second] == 0.0){
            mod.weights[d->first] = 0.0;
        }
    }
    double ws = 0.0;
    for (int j = 0; j < mod.n_components; ++j){
        ws += mod.weights[j];    
    }
    fprintf(stderr, "weightsum %f\n", ws);
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
                    solver.add_equation(doublet_props, ddoublet_props, mod.dists[i].params[0][j], mod.weights[i]);
                    solver.add_data("idx", j);
                    solver.add_data("target", dim1);
                    solver.add_data("target2", dim2);
                }     
            }
        }
        else if (ddim.count(i) > 0 && mod.weights[i] > 0.0){
            int dim = ddim[i];
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                solver.add_equation(singlet_props, dsinglet_props, mod.dists[i].params[0][j], mod.weights[i]);
                solver.add_data("idx", j);
                solver.add_data("target", dim);
                solver.add_data("target2", -1);
            }
        }
        else if (ddim.count(i) == 0 && mod.weights[i] > 0.0){
            // Background dist.
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                solver.add_equation(bg_props, dbg_props, mod.dists[i].params[0][j], mod.weights[i]);
                solver.add_data("idx", j);
                solver.add_data("target", -1);
                solver.add_data("target2", -1);
            }
        }
    }
    solver.solve();
    
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
        params[i] = solver.results[i];
    }
    fprintf(stderr, "\n");
    
    for (map<int, pair<int, int> >::iterator dd = ddoub.begin(); dd != ddoub.end(); ++dd){
        int idx1 = ddim[dd->second.first];
        int idx2 = ddim[dd->second.second];
        
        double fracbg1 = solver.results[ndim+idx1];
        double fracbg2 = solver.results[ndim+idx2];

        double mixprop1 = solver.results[ndim*2+idx1]/(solver.results[ndim*2+idx1] + solver.results[ndim*2+idx2]);
        double mixprop2 = solver.results[ndim*2+idx2]/(solver.results[ndim*2+idx1] + solver.results[ndim*2+idx2]);

        for (int i = 0; i < ndim; ++i){
            if (i == idx1){
                mod.dists[dd->first].params[0][i] = mixprop1*(1.0 - fracbg1 + fracbg1*solver.results[i]) + 
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
        
        mod.dists[dd->first].params[1][0] = solver.results[ndim*2+idx1] + solver.results[ndim*2+idx2];
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

void cb_new(mixtureModel& mod, vector<double>& params){
    
    vector<int> idx1;
    vector<int> idx2;
    vector<double> frac;
     
    /* 
    optimML::multivar_sys_solver solver(params);
    for (int i = 0; i < params.size(); ++i){
        solver.constrain_01(i);
    }
    */
    params.clear();
    for (int i = 0; i < mod.dists[0].params[0].size(); ++i){
        params.push_back(mod.dists[i].params[1][0]);
    }
    optimML::multivar_sys_solver solver(params);
    for (int i = 0; i < mod.dists[0].params[0].size(); ++i){
        solver.constrain_pos(i);
    }

    //optimML::multivar_sys_solver solver( vector<double>{ });
    //solver.add_param_grp(params);
    
    vector<double> bg_props;
    for (int i = 0; i < mod.dists[0].params[0].size(); ++i){
        bg_props.push_back(0.0);
    }
    /*
    for (map<int, pair<int, int> >::iterator d = ddoub.begin(); d != ddoub.end(); ++d){
        if (mod.weights[ddim[d->second.first]] == 0.0 || mod.weights[ddim[d->second.second]] == 0.0){
            mod.weights[d->first] = 0.0;
        }
    } 
    */
    double ws = 0.0;
    for (int i = 0; i < mod.n_components; ++i){
        ws += mod.weights[i];
    }
    fprintf(stderr, "ws %f\n", ws);
    for (int i = 0; i < mod.n_components; ++i){
        mod.weights[i] /= ws;
    }

    for (int i = 0; i < mod.dists.size(); ++i){
        if (ddoub.count(i) > 0){
            int dim1 = ddim[ddoub[i].first];
            int dim2 = ddim[ddoub[i].second];

            double prop1 = mod.dists[i].params[0][dim1];
            double prop2 = mod.dists[i].params[0][dim2];
            double fracbg = 1.0 - prop1-prop2;
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                if (j != dim1 && j != dim2){
                    bg_props[j] += mod.weights[i] * (mod.dists[i].params[0][j] / fracbg);
                }
            }
            
            if (mod.weights[i] > 0.0){
                solver.add_equation(doub_prop, ddoub_prop, prop1/(prop1+prop2));
                solver.add_data("idx1", dim1);
                solver.add_data("idx2", dim2);
            }
        }
        else if (ddim.count(i) > 0 && mod.weights[i] > 0.0){
            int dim = ddim[i];
            double fracbg = 1.0 - mod.dists[i].params[0][dim];
            for (int j = 0; j < mod.dists[i].params[0].size(); ++j){
                if (j != dim){
                    bg_props[j] += mod.weights[i] * (mod.dists[i].params[0][j] / fracbg);
                }
            }
        }
    }
    solver.solve();
    
    for (map<int, pair<int, int> >::iterator dd = ddoub.begin(); dd != ddoub.end(); ++dd){
        int idx1 = ddim[dd->second.first];
        int idx2 = ddim[dd->second.second];

        double fracbg = 1.0 - mod.dists[dd->first].params[0][idx1] - mod.dists[dd->first].params[0][idx2];

        for (int i = 0; i < mod.dists[dd->first].params[0].size(); ++i){
            if (i != idx1 && i != idx2){
                mod.dists[dd->first].params[0][i] = fracbg*bg_props[i];
            }
        }

        double newfrac1 = solver.results[idx1] / (solver.results[idx1] + solver.results[idx2]);
        double newfrac2 = 1.0 - newfrac1;
        newfrac1 *= (1.0 - fracbg);
        newfrac2 *= (1.0 - fracbg);
        mod.dists[dd->first].params[0][idx1] = newfrac1;
        mod.dists[dd->first].params[0][idx2] = newfrac2;
        
        mod.dists[dd->first].params[1][0] = params[idx1] + params[idx2];
    }
    for (map<int, int>::iterator dd = ddim.begin(); dd != ddim.end(); ++dd){
        double fracbg = 1.0 - mod.dists[dd->second].params[0][dd->second];
        for (int i = 0; i < mod.dists[dd->second].params[0].size(); ++i){
            if (i != dd->second){
                mod.dists[dd->second].params[0][i] = fracbg*bg_props[i];
            }
        }
        mod.dists[dd->second].params[1][0] = params[dd->second];

    }

    for (int i = 0; i < bg_props.size(); ++i){
        fprintf(stderr, "BGPROP %d) %f\n", i, bg_props[i]);
    }
    fprintf(stderr, "\n");
    for (int i = 0; i < solver.results.size(); ++i){
        fprintf(stderr, "NREADS %d) %f\n", i, solver.results[i]);
    }


}

double size_pbg(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double size = data_d.at("size");
    
    return params[0]* size + params[1];
     
    //return params[0] * size*size + params[1] * size + params[1];

    //return params[0] * size + params[1];
}

void dsize_pbg(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){

    double size = data_d.at("size");
    
    /*
    results[0] += 2*size;
    results[1] += size;
    results[2] += 1.0;
    */

    results[0] += size;
    results[1] += 1.0;
}

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

double ll_withtot(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double tot = data_d.at("tot");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    int ndim = data_i.at("ndim");

    double a = params[params.size()-2];
    double b = params[params.size()-1];
    
    double f = 1.0 - 1/(exp(a)*pow(tot, b) + 1.0);
    
    vector<double> x;
    vector<double> p;
    
    static char buf[50];
    static string bufstr;
    for (int idx = 0; idx < ndim; ++idx){
        sprintf(&buf[0], "n_%d", idx);
        bufstr = buf;
        double count = data_d.at(bufstr);

        double p_bg = params[idx];
        if (idx2 == -1){
            if (idx1 == idx){
                p.push_back(p_bg * f + (1.0 - f));
            }
            else{
                p.push_back(p_bg * f);
            }
        }
        else{
            double num1 = params[ndim + idx1];
            double num2 = params[ndim + idx2];

            if (idx == idx1 || idx == idx2){
                double m = num1/(num1+num2);
                if (idx == idx2){
                    m = 1.0 - m;
                }
                p.push_back(p_bg * f + (1.0 - f)*m);
            }
            else{
                p.push_back(p_bg * f);
            }
        }
        x.push_back(count);
    }
    return dmultinom(x, p)/log2(exp(1));
}

void dll_withtot(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
 
    double tot = data_d.at("tot");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    int ndim = data_i.at("ndim");
    
    double a = params[params.size()-2];
    double b = params[params.size()-1];
    
    double f = 1.0 - 1/(exp(a)*pow(tot, b) + 1.0);
    
    double df_da = (exp(a)*pow(tot, b))/(pow(exp(a)*pow(tot, b) + 1, 2));
    double df_db = df_da * log(tot);

    double lambda = tot;
    
    static char buf[50];
    static string bufstr;
    for (int idx = 0; idx < ndim; ++idx){
        double p;
        double p_bg = params[idx];
        sprintf(&buf[0], "n_%d", idx);
        bufstr = buf;
        double count = data_d.at(bufstr);
        if (idx2 == -1){
            if (idx1 == idx){
                p = (p_bg * f + (1.0 - f));
                double dLL_dp = (count/p - lambda);
                results[idx] += dLL_dp*f;
                double dp_df = p_bg - 1.0;
                results[results.size()-2] += dLL_dp*dp_df*df_da;
                results[results.size()-1] += dLL_dp*dp_df*df_db;
            }
            else{
                p = (p_bg * f);
                double dLL_dp = (count/p - lambda)*f;
                results[idx] += dLL_dp*f;
                double dp_df = p_bg;
                results[results.size()-2] += dLL_dp*dp_df*df_da;
                results[results.size()-1] += dLL_dp*dp_df*df_db;
            }
        }
        else{
            double num1 = params[ndim + idx1];
            double num2 = params[ndim + idx2];

            if (idx == idx1 || idx == idx2){
                double m = num1/(num1+num2);
                if (idx == idx2){
                    m = 1.0 - m;
                }
                p = (p_bg * f + (1.0 - f)*m);
                double dLL_dp = (count/p - lambda)*f;
                results[idx] += dLL_dp*f;
                double dp_df = (p_bg - m);
                results[results.size()-2] += dLL_dp*dp_df*df_da;
                results[results.size()-1] += dLL_dp*dp_df*df_db;
                double denom = pow(num1+num2, 2);
                if (idx == idx1){
                    results[ndim+idx1] += (-num2*(f-1))/denom;
                    results[ndim+idx2] += (num1*(f-1))/denom;
                }
                else{
                    results[ndim+idx2] += (-num2*(f-1))/denom;
                    results[ndim+idx1] += (num1*(f-1))/denom;
                }
            }
            else{
                p = (p_bg * f);
                double dLL_dp = (count/p - lambda)*f;
                results[idx] += dLL_dp*f;
                double dp_df = p_bg;
                results[results.size()-2] += dLL_dp*dp_df*df_da;
                results[results.size()-2] += dLL_dp*dp_df*df_db;
            }
        }
    }
}

void assign_ids_yetagain(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned,
    vector<string>& labels){
    
    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }
    
    vector<vector<double> > obs;
    vector<unsigned long> bcs;

    vector<double> tots;
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = bc_ms_counts.begin(); 
        x != bc_ms_counts.end(); ++x){
        
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
        row.push_back(tot);
        obs.push_back(row);
        bcs.push_back(x->first);
    }
    
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
        
        mixtureDist dist(vector<string>{ "multinomial", "normal" }, vector<vector<double> >{ params, vector<double>{ totmean, totsd } });
        dist.set_num_inputs(0, labels.size());

        dist.name = labels[i];
        int didx = dists.size();
        ddim.insert(make_pair(didx, i));
        dists.push_back(dist);
    }
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
            mixtureDist dist2(vector<string>{ "multinomial", "normal" }, vector<vector<double> >{ params2, vector<double>{ 2*totmean, totsd_doub }});
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
    
    vector<vector<vector<double> > > obsdist;
    vector<double> weights;
    vector<int> idx1;
    vector<int> idx2;
    vector<double> tots_fit;
    vector<vector<double> > countsdat;

    for (int i = 0; i < labels.size(); ++i){
        vector<double> v;
        countsdat.push_back(v);
    }

    double prob_cutoff = 0.001;

    for (int i = 0; i < obs.size(); ++i){
        vector<double> lls;
        double llmax = 0.0;
        for (int j = 0; j < dists.size(); ++j){
            double ll = dists[j].loglik(obs[i]);    
            if (llmax == 0.0 || ll > llmax){
                llmax = ll;
            }
            lls.push_back(ll);
        }
        double lltot = 0.0;
        for (int j = 0; j < lls.size(); ++j){
            lltot += (lls[j] - llmax);
        }
        for (int j = 0; j < lls.size(); ++j){
            double prob = pow(2, (lls[j] - llmax) - lltot);
            
            if (prob > prob_cutoff){        
                weights.push_back(prob);
                if (ddoub.count(j) > 0){
                    idx1.push_back(ddim[ddoub[j].first]);
                    idx2.push_back(ddim[ddoub[j].second]);
                }
                else{
                    idx1.push_back(j);
                    idx2.push_back(-1);
                }
                tots_fit.push_back(tots[i]); 
                for (int k = 0; k < labels.size(); ++k){
                    countsdat[k].push_back(obs[i][k]);
                }                    
            }
        }

    }
    
    int ndim = (int)labels.size();
    vector<double> params;
    // Fraction in background
    for (int i = 0; i < ndim; ++i){
        params.push_back(1.0/(double)ndim);
    }
    optimML::multivar_ml_solver solver(vector<double>{ }, ll_withtot, dll_withtot);
    solver.add_param_grp(params);
    // Mean size
    for (int i = 0; i < ndim; ++i){
        solver.add_one_param(totmean);
        solver.constrain_pos(ndim+i);
    }
    // Regression coeffs
    solver.add_one_param(-1.0);
    solver.add_one_param(1.0);
   
    solver.add_data("idx1", idx1);
    solver.add_data("idx2", idx2);
    solver.add_data("tot", tots_fit);
    char buf[50];
    for (int i = 0; i < labels.size(); ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        solver.add_data(bufstr, countsdat[i]);
    }
    solver.add_data_fixed("ndim", (int)labels.size());
    solver.solve();
    fprintf(stderr, "DONE???\n");
    for (int i = 0; i < solver.results.size(); ++i){
        fprintf(stderr, "%d) %f\n", i, solver.results[i]);
    }
    exit(0);
}

void assign_ids_brandnew(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned,
    vector<string>& labels){
    
    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }
    
    vector<vector<double> > obs;
    vector<unsigned long> bcs;

    vector<double> tots;
    
    vector<vector<double> > obs2;

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = bc_ms_counts.begin(); 
        x != bc_ms_counts.end(); ++x){
        
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
        
        //mixtureDist dist("multinomial", params);
        //dist.set_num_inputs(labels.size());
        mixtureDist dist(vector<string>{ "multinomial", "normal" }, vector<vector<double> >{ params, vector<double>{ totmean, totsd } });
        dist.set_num_inputs(0, labels.size());

        dist.name = labels[i];
        int didx = dists.size();
        ddim.insert(make_pair(didx, i));
        dists.push_back(dist);
    }
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
            //mixtureDist dist2("multinomial", params2);
            //dist2.set_num_inputs(labels.size());
            mixtureDist dist2(vector<string>{ "multinomial", "normal" }, vector<vector<double> >{ params2, vector<double>{ 2*totmean, totsd_doub }});
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

    mixtureModel mod(dists);
    mod.set_callback(cb_new_take2, shared_params);
    mod.fit(obs);
    //cb_new_take2(mod, shared_params);
    mod.print();
    
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
        double tot = obs[i][labels.size()];
        vector<vector<double> > x1;
        for (int j = 0; j < labels.size(); ++j){
            x1.push_back(vector<double>{ obs[i][j] });
        } 
        int assn = mod.assignments[i];
        if (ddoub.count(assn) == 0){
            int dim = ddim[assn];
            
            optimML::brent_solver mix(ll_mix_new, dll_mix_new);
            mix.add_data_fixed("ndim", (int)labels.size());
            mix.constrain_01();
            vector<int> idx1{ dim };
            mix.add_data("idx1", idx1);
            for (int k = 0; k < labels.size(); ++k){
                mix.add_data(datakeys[k], x1[k]);
                mix.add_data(bgpkeys[k], bgfrac[k]);
            }    
            double res = mix.solve(0, 1);
            totcounts.push_back(log(tot));
            pbgs.push_back(log(res)-log(1-res));

        }
        else{
            int dim1 = ddim[ddoub[assn].first];
            int dim2 = ddim[ddoub[assn].second];
            optimML::brent_solver mix2(ll_mix_new, dll_mix_new);
            mix2.add_data_fixed("ndim", (int)labels.size());
            mix2.constrain_01();
            vector<int> idxes1{ dim1 };
            vector<int> idxes2{ dim2 };
            mix2.add_data("idx1", idxes1);
            mix2.add_data("idx2", idxes2);
            int ndim = labels.size();
            vector<double> mixprop{ mod.shared_params[ndim*2 + dim1] / 
                (mod.shared_params[ndim*2 + dim1] + mod.shared_params[ndim*2 + dim2]) };
            mix2.add_data("mixprop", mixprop);
            for (int k = 0; k < labels.size(); ++k){
                if (mod.weights[k] > 0.0){
                    mix2.add_data(datakeys[k], x1[k]);
                    mix2.add_data(bgpkeys[k], bgfrac[k]);
                }
            }
            double res2 = mix2.solve(0,1);
            totcounts.push_back(log(tot));
            pbgs.push_back(log(res2) - log(1-res2));
        }

        continue; 



        vector<vector<double> > x;
        vector<vector<double> > bgp;
        for (int j = 0; j < labels.size(); ++j){
            x.push_back(vector<double>{ obs[i][j] });
            bgp.push_back(vector<double>{ mod.shared_params[j] });
        }
        
        vector<pair<double, int> > lls;
        map<int, double> mixp;

        string bcstr = bc2str(bcs[i]);
        bool print = bcstr == "AAAGTCCAGATCCCAT";
        //bool print = bcstr == "CGCGTGACATCTGCGG";

        for (int j = 0; j < labels.size(); ++j){
            
            if (mod.weights[j] == 0.0){
                continue;
            }
            
            /*
            optimML::brent_solver mixbg(ll_mix_bg_binom, dll_mix_bg_binom);
            mixbg.constrain_01();
            vector<double> n;
            vector<double> x;
            vector<double> bg;
            for (int k = 0; k < labels.size(); ++k){
                if (mod.weights[k] > 0.0 && k != j){
                    n.push_back(tot);
                    x.push_back(obs[i][k]);
                    bg.push_back(mod.shared_params[k]);
                }
            } 
            mixbg.add_data("n", n);
            mixbg.add_data("k", x);
            mixbg.add_data("bg", bg);
            double p = mixbg.solve(0,1);
            vector<double> params;
            vector<double> counts;
            for (int k = 0; k < labels.size(); ++k){
                if (mod.weights[k] > 0.0){
                    counts.push_back(obs[i][k]);
                    if (k == j){
                        params.push_back(mod.shared_params[k]*p + (1.0-p));
                    }
                    else{
                        params.push_back(mod.shared_params[k]*p);
                    }
                }
            }
            double ll = dmultinom(counts, params);
            lls.push_back(make_pair(-ll, j));
            mixp.push_back(p);
            if (print){
                fprintf(stderr, "%s\t%f\t%f\n", labels[j].c_str(), p, ll);
            }
            continue;
            */
             
            // Test if this one is the true identity
            optimML::brent_solver mix(ll_mix_new, dll_mix_new);
            mix.add_data_fixed("ndim", (int)labels.size());
            mix.constrain_01();
            vector<int> idx1{ j };
            mix.add_data("idx1", idx1);
            for (int k = 0; k < labels.size(); ++k){
                if (mod.weights[k] > 0.0){
                    mix.add_data(datakeys[k], x[k]);
                    mix.add_data(bgpkeys[k], bgp[k]);
                }               
            }    
            double res = mix.solve(0, 1);
            
            double ll = mix.log_likelihood + dnorm(obs[i][obs[i].size()-1], 
                dists[j].params[1][0], dists[j].params[1][1]);
            ll = mix.log_likelihood;
            //ll += dnorm(obs[i][obs[i].size()-1], dists[j].params[1][0],
            //    dists[j].params[1][1]);
            //ll += dpois(obs[i][obs[i].size()-1], dists[j].params[1][0]);
            if (print){
                fprintf(stderr, "%s\t%f\t%f\n", labels[j].c_str(), res, ll);
            }
            lls.push_back(make_pair(-ll, j));
            mixp.insert(make_pair(j,res));
            
        }
        sort(lls.begin(), lls.end());
        int idx1 = lls[0].second;
        
        vector<pair<double, int> > lls2;
        map<int, double> mixp2;

        for (int j = 0; j < labels.size(); ++j){
            if (mod.weights[j] > 0.0 && j != idx1){
                /*    
                optimML::brent_solver mixbg2(ll_mix_bg_binom, dll_mix_bg_binom);
                mixbg2.constrain_01();
                vector<double> n;
                vector<double> x;
                vector<double> bg;
                for (int k = 0; k < labels.size(); ++k){
                    if ( mod.weights[k] > 0.0 && k != idx1 && k != j){
                        x.push_back(obs[i][k]);
                        n.push_back(tot);
                        bg.push_back(mod.shared_params[k]);
                    }
                }
                mixbg2.add_data("n", n);
                mixbg2.add_data("k", x);
                mixbg2.add_data("bg", bg);
                double res = mixbg2.solve(0,1);
                vector<double> params;
                vector<double> counts;

                int nd = labels.size();
                double mxp = mod.shared_params[nd*2 + idx1] / (mod.shared_params[nd*2 + idx1] + mod.shared_params[nd*2 + j]);
                for (int k = 0; k < labels.size(); ++k){
                    if (mod.weights[k] > 0.0){
                        counts.push_back(obs[i][k]);
                        if (k == idx1){
                            params.push_back(res*mod.shared_params[k] + (1.0-res)*mxp);
                        }
                        else if (k == j){
                            params.push_back(res*mod.shared_params[k] + (1.0-res)*(1.0-mxp));
                        }
                        else{
                            params.push_back(res*mod.shared_params[k]);
                        }
                    }
                }
                
                double ll = dmultinom(counts, params);
                lls2.push_back(make_pair(-ll, j));
                mixp2.push_back(res);
                if (print){
                    fprintf(stderr, "%s+%s\t%f\t%f\n", labels[idx1].c_str(), labels[j].c_str(),
                        res, ll);
                }
                continue;
                */
                  
                optimML::brent_solver mix2(ll_mix_new, dll_mix_new);
                mix2.add_data_fixed("ndim", (int)labels.size());
                mix2.constrain_01();
                vector<int> idxes1{ idx1 };
                vector<int> idxes2{ j };
                mix2.add_data("idx1", idxes1);
                mix2.add_data("idx2", idxes2);
                int ndim = labels.size();
                vector<double> mixprop{ mod.shared_params[ndim*2 + idx1] / 
                    (mod.shared_params[ndim*2 + idx1] + mod.shared_params[ndim*2 + j]) };
                mix2.add_data("mixprop", mixprop);
                for (int k = 0; k < labels.size(); ++k){
                    if (mod.weights[k] > 0.0){
                        mix2.add_data(datakeys[k], x[k]);
                        mix2.add_data(bgpkeys[k], bgp[k]);
                    }
                }
                double res2 = mix2.solve(0,1);
                
                pair<int, int> key;
                if (idx1 < j){
                    key = make_pair(idx1, j);
                }
                else{
                    key = make_pair(j, idx1);
                }
                double ll = mix2.log_likelihood + dnorm(obs[i][obs[i].size()-1], 
                    mod.dists[doub2dist[key]].params[1][0],
                    mod.dists[doub2dist[key]].params[1][1]);
                ll = mix2.log_likelihood;
                //ll += dpois(obs[i][obs[i].size()-1], mod.dists[doub2dist[key]].params[1][0]);
                //ll += dnorm(obs[i][obs[i].size()-1], mod.dists[doub2dist[key]].params[1][0],
                //    mod.dists[doub2dist[key]].params[1][1]);
                if (print){
                    fprintf(stderr, "%s+%s\t%f\t%f\n", labels[idx1].c_str(), labels[j].c_str(), res2, ll);
                }
                lls2.push_back(make_pair(-ll, j));
                mixp2.insert(make_pair(j, res2));
            }
        }
        sort(lls2.begin(), lls2.end());
        double llr;
        string name;
        char s_d;
        double pprint;

        if (-lls2[0].first > -lls[0].first){
            s_d = 'D';
            // Doublet.
            pprint = mixp2[lls2[0].second];
            llr = -lls2[0].first - -lls2[1].first;
            if (-lls[0].first > -lls2[1].first){
                llr = -lls2[0].first - -lls[0].first;
            }
            if (labels[lls2[0].second] < labels[idx1]){
                name = labels[lls2[0].second] + "+" + labels[idx1];
            }
            else{
                name = labels[idx1] + "+" + labels[lls2[0].second];
            }
        }
        else{
            pprint = mixp[lls[0].second];
            s_d = 'S';
            llr = -lls[0].first - -lls[1].first;
            if (-lls2[0].first > -lls[1].first){
                llr = -lls[0].first - -lls2[0].first;
            }
            name = labels[idx1];
        }
        fprintf(stdout, "%s\t%s\t%c\t%f\t%f\n", bcstr.c_str(), name.c_str(), s_d, llr, pprint);
    } 


    pair<double, double> a_b = linreg(totcounts, pbgs);
    fprintf(stderr, "%f %f\n", a_b.first, a_b.second);
    /*
    for (int i = 0; i < tots.size(); ++i){
        double pred = 1.0 - 1/(exp(a_b.first) * pow(exp(totcounts[i]), a_b.second) + 1.0);
        pred = log(pred) - log(1-pred);
        fprintf(stdout, "%f\t%f\t%f\n", totcounts[i], pbgs[i], pred);
    }
    */
    for (int i = 0; i < obs.size(); ++i){
        double tot = obs[i][obs[i].size()-1];
        double pbg = a_b.first * log(tot) + a_b.second;
        pbg = exp(pbg)/(exp(pbg) + 1);
        vector<pair<double, string> > lls;
        for (int j = 0; j < mod.n_components; ++j){
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
            for (int k = 0; k < labels.size(); ++k){
                if (dim2 == -1){
                    if (k == dim1){
                        params.push_back(pbg*mod.shared_params[k] + (1.0 - pbg));
                    }
                    else{
                        params.push_back(pbg*mod.shared_params[k]);
                    }
                }
                else{
                    double m1 = mod.shared_params[labels.size()*2 + dim1];
                    double m2 = mod.shared_params[labels.size()*2 + dim2];
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
            double ll = dmultinom(obs2[i], params);
            lls.push_back(make_pair(-ll, mod.dists[j].name));
        }
        vector<double> params_blank;
        for (int j = 0; j < labels.size(); ++j){
            params_blank.push_back(mod.shared_params[j]);
        }
        double ll = dmultinom(obs2[i], params_blank);
        lls.push_back(make_pair(-ll, ""));
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0 && lls[0].second != ""){
            char s_d = 'S';
            for (int x = 0; x < lls[0].second.size(); ++x){
                if (lls[0].second[x] == '+'){
                    s_d = 'D';
                    break;
                }
            }
            string bcstr = bc2str(bcs[i]);
            fprintf(stdout, "%s\t%s\t%c\t%f\n", bcstr.c_str(), lls[0].second.c_str(), s_d, llr);
        }
    }

    exit(0);

    /*
    for (int i = 0; i < obs.size(); ++i){
        string bcstr = bc2str(bcs[i]);
        fprintf(stdout, "%s", bcstr.c_str());
        double tot = obs[i][labels.size()];
        for (int j = 0; j < labels.size(); ++j){
            double x = obs[i][j] / tot / mod.shared_params[j];
            fprintf(stdout, "\t%f", x);
        }
        fprintf(stdout, "\n");
    } 
    exit(0);
    */

    // Find count / % background relationship
    //optimML::multivar_sys_solver pbg( vector<double>{ 0.1, 1.0 } );
    
    /* 
    map<int, vector<double> > logcounts;
    map<int, vector<double> > logitrhs;
    map<int, vector<double> > inferredp;
    map<int, vector<double> > doublet_f;
    map<int, vector<unsigned long> > logcount_bc;
    
    map<int, pair<double, double> > regression_coeffs;
    */

    vector<double> logcounts;
    vector<double> logitrhs;
    //vector<double> p_weights;
    vector<double> inferredp;
    vector<unsigned long> logcount_bc;
    
    for (int i = 0; i < obs.size(); ++i){
        
        int jmax = -1;
        double maxp = 0.0;
        for (int j = 0; j < mod.n_components; ++j){
            if (jmax == -1 || mod.responsibility_matrix[i][j] > maxp){
                jmax = j;
                maxp = mod.responsibility_matrix[i][j];
            }
        }
        if (ddoub.count(jmax) == 0 && ddim.count(jmax) > 0 && mod.weights[jmax] > 0 && mod.dists[jmax].name != "BG"){
            /*
            if (logcounts.count(jmax) == 0){
                vector<double> v;
                logcounts.insert(make_pair(jmax, v));
                logitrhs.insert(make_pair(jmax, v));
                inferredp.insert(make_pair(jmax, v));
                vector<unsigned long> v2;
                logcount_bc.insert(make_pair(jmax, v2));
            }
            */

            // Singlet. Infer mixing proportion.
            int dim1 = ddim[jmax];
            int ndim = labels.size();

            // Solve for mixprop
            optimML::multivar_ml_solver prop1(vector<double>{ 1.0 - mod.shared_params[ndim+dim1] }, 
                ll_mix, dll_mix);
            prop1.constrain_01(0);
            prop1.add_data_fixed("ndim", (int)labels.size());
            
            vector<int> prop1_idx1{ dim1 };
            vector<int> prop1_idx2{ -1 };                
            prop1.add_data("idx1", prop1_idx1);
            prop1.add_data("idx2", prop1_idx2);
            vector<vector<double> > x;
            vector<vector<double> > bgp;

            for (int j = 0; j < labels.size(); ++j){
                vector<double> x_i { obs[i][j] };
                x.push_back(x_i);
                vector<double> bgp_i { mod.shared_params[j ] };
                bgp.push_back(bgp_i);
            }

            char buf[50];
            for (int j = 0; j < labels.size(); ++j){
                sprintf(&buf[0], "n_%d", j);
                string bufstr = buf;
                prop1.add_data(bufstr, x[j]);
                sprintf(&buf[0], "bg_%d", j);
                bufstr = buf;
                prop1.add_data(bufstr, bgp[j]);
            }
            prop1.solve();
            double count = obs[i][dim1];
            
            //string xx = bc2str(bcs[i]);
            //fprintf(stdout, "%s\t%s\t%f\t%f\t%f\n", xx.c_str(), labels[dim1].c_str(), count, 1.0-prop1.results[0], mod.shared_params[ndim+dim1]);

            //pbg.add_equation(size_pbg, dsize_pbg, log(prop1.results[0]) - log(1.0 - prop1.results[0]), maxp);
            //pbg.add_data("size", log(count));
            //p_weights.push_back(maxp);
            
            logcounts.push_back(log(count));
            double fracbg = 1.0-prop1.results[0];
            logitrhs.push_back(log(fracbg)-log(1.0-fracbg));
            inferredp.push_back(1.0-prop1.results[0]);
            logcount_bc.push_back(bcs[i]);
            
            /*
            logcounts[jmax].push_back(log(count));
            double fracbg = 1.0-prop1.results[0];
            logitrhs[jmax].push_back(log(fracbg) - log(1.0-fracbg));
            //logrhs.push_back(log(1.0-prop1.results[0])); 
            inferredp[jmax].push_back(1.0-prop1.results[0]);
            logcount_bc[jmax].push_back(bcs[i]);
            */
        }
        /*
        else if (ddoub.count(jmax) >  0 && mod.dists[jmax].name != "BG" && mod.weights[jmax] > 0){
            if (logcounts.count(jmax) == 0){
                vector<double> v;
                logcounts.insert(make_pair(jmax, v));
                logitrhs.insert(make_pair(jmax, v));
                inferredp.insert(make_pair(jmax, v));
                vector<unsigned long> v2;
                logcount_bc.insert(make_pair(jmax, v2));
                doublet_f.insert(make_pair(jmax, v));
            }
            
            int dim1 = ddim[ddoub[jmax].first];
            int dim2 = ddim[ddoub[jmax].second];

            int ndim = labels.size();

            // Solve for mixprop and bgfrac
            double mp = mod.shared_params[ndim*2+dim1] / (mod.shared_params[ndim*2+dim1] + mod.shared_params[ndim*2+dim2]);
            double bgf1 = mod.shared_params[ndim+dim1];
            double bgf2 = mod.shared_params[ndim+dim2];
            double bgf = mp*bgf1 + (1.0-mp)*bgf2;
            optimML::multivar_ml_solver prop(vector<double>{}, ll_mix, dll_mix);
            //vector<double> pg{ (1.0-bgf)*mp, (1.0-bgf)*(1.0-mp), bgf};
            vector<double> pg{ 1.0/3.0, 1.0/3.0, 1.0/3.0 };
            fprintf(stderr, "PG %f %f %f\n", pg[0], pg[1], pg[2]);
            prop.add_param_grp(pg);
            prop.add_data_fixed("ndim", (int)labels.size());
            
            vector<int> prop1_idx1{ dim1 };
            vector<int> prop1_idx2{ dim2 };                
            prop.add_data("idx1", prop1_idx1);
            prop.add_data("idx2", prop1_idx2);
            vector<vector<double> > x;
            vector<vector<double> > bgp;

            for (int j = 0; j < labels.size(); ++j){
                vector<double> x_i { obs[i][j] };
                x.push_back(x_i);
                vector<double> bgp_i { mod.shared_params[j ] };
                bgp.push_back(bgp_i);
            }

            char buf[50];
            for (int j = 0; j < labels.size(); ++j){
                sprintf(&buf[0], "n_%d", j);
                string bufstr = buf;
                prop.add_data(bufstr, x[j]);
                fprintf(stderr, "  %s %f\n", bufstr.c_str(), x[j][0]);
                sprintf(&buf[0], "bg_%d", j);
                bufstr = buf;
                prop.add_data(bufstr, bgp[j]);
                fprintf(stderr, "  %s %f\n", bufstr.c_str(), bgp[j][0]);
            }
            prop.solve();
            double f = prop.results[0] / (prop.results[0] + prop.results[1]);
            double bgfrac = prop.results[2];
            logcounts[jmax].push_back(log(obs[i][dim1] + obs[i][dim2]));
            logitrhs[jmax].push_back(log(bgfrac) - log(1.0-bgfrac));
            inferredp[jmax].push_back(bgfrac);
            logcount_bc[jmax].push_back(bcs[i]);
            doublet_f[jmax].push_back(f);
        }
        */
    }
    
    pair<double, double> res = linreg(logcounts, logitrhs);
    fprintf(stderr, "RES %f %f\n", res.first, res.second);

    /*
    for (map<int, vector<double> >::iterator p = inferredp.begin(); p != inferredp.end(); ++p){
        fprintf(stderr, "NUM %ld\n", logcounts[p->first].size());
        pair<double, double> res = linreg(logcounts[p->first], logitrhs[p->first]);
        fprintf(stderr, "%s\t%f\t%f\n", mod.dists[p->first].name.c_str(), res.first, res.second);
        regression_coeffs.insert(make_pair(p->first, res));
    } 
    exit(1);
    */
    /*
    for (map<int, vector<double> >::iterator df = doublet_f.begin(); df != doublet_f.end(); ++df){
        for (int i = 0; i < df->second.size(); ++i){
            string bcst = bc2str(logcount_bc[df->first][i]);
            fprintf(stdout, "%s\t%s\t%f\t%f\n", bcst.c_str(), mod.dists[df->first].name.c_str(),
                logcounts[df->first][i], inferredp[df->first][i]);
        }
    }
    exit(0);
    
    pair<double, double> res = linreg(logcounts, logrhs);
    fprintf(stderr, "res %f %f\n", res.first, res.second);
    */ 

    // Extract all global parameters
    vector<double> bg_prof;
    for (int i = 0; i < labels.size(); ++i){
        bg_prof.push_back(mod.shared_params[i]);
    }
    vector<string> types{ "multinomial", "normal" };
    for (int i = 0; i < obs.size(); ++i){
        vector<mixtureDist> dists;
        vector<pair<double, string> > lls;
        vector<string> distnames;
        int cur_idx = 0;

        for (int j = 0; j < labels.size(); ++j){
            // Infer percent background
            double count = obs[i][j];
            if (count > 0.0 && mod.weights[cur_idx] > 0.0){
                double perc_bg = 1.0 - 1.0/(exp(res.second) * pow(count, res.first) + 1.0);
                
                fprintf(stderr, "count %f pbg %f\n", count, perc_bg);
                vector<double> params;
                for (int k = 0; k < labels.size(); ++k){
                    if (k == j){
                        params.push_back((1.0 - perc_bg) + perc_bg*bg_prof[k]);
                    }
                    else{
                        params.push_back(perc_bg * bg_prof[k]);
                    }
                }
                //mixtureDist dist(types, vector<vector<double> >{ params, vector<double>{ mod.dists[j].params[1][0],
                //   mod.dists[j].params[1][1] } });
                mixtureDist dist(vector<string>{ "multinomial" }, params);
                dist.set_num_inputs(0, labels.size());
                dist.name = labels[j];
                distnames.push_back(dist.name);
                //double ll = dist.loglik(obs2[i]);
                double ll = dmultinom(obs2[i], params);
                lls.push_back(make_pair(-ll, dist.name));
            }
            cur_idx++;
        }    

        for (int j = 0; j < labels.size()-1; ++j){
            double perc_bg_j = 1.0 - 1.0/(exp(res.second) * pow(obs[i][j], res.first) + 1.0);
            for (int k = j + 1; k < labels.size(); ++k){
                if (obs[i][j] > 0.0 && obs[i][k] > 0.0 && mod.weights[cur_idx] > 0.0){
                    double perc_bg_k = 1.0 - 1.0/(exp(res.second) * pow(obs[i][k], res.first) + 1.0);
                    
                    double perc_bg = 1.0 - 1.0/(exp(res.second) * pow(obs[i][j] + obs[i][k], res.first) + 1.0);
                    
                    //double frac = (obs[i][j])/(obs[i][j] + obs[i][k]);
                    int ndim = labels.size();
                    double frac = mod.shared_params[ndim*2+j] / (mod.shared_params[ndim*2+j] + mod.shared_params[ndim*2+k]);
                    //frac = obs[i][j]/(obs[i][j] + obs[i][k]);

                    vector<double> params;
                    for (int l = 0; l < labels.size(); ++l){
                        if (l == j){
                            params.push_back(perc_bg*bg_prof[l] + frac*(1.0-perc_bg));

                            //params.push_back(frac*perc_bg_j*bg_prof[l] + frac*(1.0-perc_bg_j) + 
                            //    (1.0-frac)*perc_bg_k*bg_prof[l]);
                        }
                        else if (l == k){
                            params.push_back(perc_bg*bg_prof[l] + (1.0-frac)*(1.0-perc_bg));

                            //params.push_back(frac*perc_bg_j*bg_prof[l] + (1.0-frac)*perc_bg_k*bg_prof[l] + 
                            //    (1.0-frac)*(1.0-perc_bg_k));
                        }
                        else{
                            params.push_back(frac*perc_bg*bg_prof[l]);

                            //params.push_back(frac*perc_bg_j*bg_prof[l] + (1.0-frac)*perc_bg_k*bg_prof[l]);
                        }
                    }
                    
                    int dist_idx = lls.size();
                    //mixtureDist dist(types, vector<vector<double> >{ params, 
                    //    vector<double>{ mod.dists[cur_idx].params[1][0],
                    //    mod.dists[cur_idx].params[1][1] } });
                    mixtureDist dist(vector<string>{ "multinomial" },
                        params);
                    dist.set_num_inputs(0, labels.size());
                    if (labels[j] < labels[k]){
                        dist.name = labels[j] + "+" + labels[k];
                    }
                    else{
                        dist.name = labels[k] + "+" + labels[j];
                    }
                    
                    string bs = bc2str(bcs[i]);
                    if (false && bs == "AGGTGTTAGTAGTCAA"){
                    //if (false && dist.name == "CMO304+CMO307" && bs == "AAACCCAAGTCCTGCG"){
                        fprintf(stderr, "here\n");
                        for (int j = 0; j < obs[i].size()-1; ++j){
                            fprintf(stderr, "obs[%d] = %f\n", j, obs[i][j]);
                        }
                        fprintf(stderr, "frac %f\n", frac);
                        fprintf(stderr, "bpg j %f k %f\n", perc_bg_j, perc_bg_k);
                        fprintf(stderr, "allparams\n");
                        for (int j = 0; j < params.size(); ++j){
                            fprintf(stderr, "%d) %f\n", j, params[j]);
                        }
                        fprintf(stderr, "sizep %f %f\n", mod.dists[cur_idx].params[1][0], mod.dists[cur_idx].params[1][1]);
                        fprintf(stderr, "LLS\n");
                        for (int j = 0; j < lls.size(); ++j){
                            fprintf(stderr, "%s: %f\n", lls[j].second.c_str(), -lls[j].first);
                        }
                        exit(1);
                    }

                    distnames.push_back(dist.name);
                    //double ll = dist.loglik(obs2[i]);
                    double ll = dmultinom(obs2[i], params);
                    lls.push_back(make_pair(-ll, dist.name));
                }
                cur_idx++;
            }
        }
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0.0){
            string bcstr = bc2str(bcs[i]);
            char s_d = 'S';
            for (int x = 0; x < lls[0].second.size(); ++x){
                if (lls[0].second[x] == '+'){
                    s_d = 'D';
                    break;
                }
            }
            fprintf(stdout, "%s\t%s\t%c\t%f\n", bcstr.c_str(), lls[0].second.c_str(),
                s_d, llr);
        }
    }
    
    exit(0);
    /*
    pbg.solve();

    fprintf(stderr, "%f %f\n", pbg.results[0], pbg.results[1]);
    */
    

    /*
    for (int i = 0; i < obs.size(); ++i){
        //double pred = exp(pbg.results[0]*log(test1[i]) + pbg.results[1]);
        double pred = exp(pbg.results[0] * pow(logcounts[i], 2) + pbg.results[1] * logcounts[i] + pbg.results[2]);
        //fprintf(stdout, "%f\t%f\t%f\n", exp(logcounts[i]), pred, exp(logrhs[i]));
    }
    */
    exit(0);

   // fprintf(stderr, "%f %f\n", pbg.results[0], pbg.results[1]);
    
    exit(0);

    for (int i = 0; i < obs.size(); ++i){
        vector<pair<double, int> > lls;
        
        for (int j = 0 ; j < mod.n_components; ++j){
            if (mod.weights[j] > 0.0){
                double ll = mod.dists[j].loglik(obs[i]);
                lls.push_back(make_pair(-ll, j));
            }
        }
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;

        if (llr > 0.0){
            string name1;
            string name2 = "";
            if (ddoub.count(lls[0].second) > 0){

                string bcstr = bc2str(bcs[i]);

                int dim1 = ddim[ddoub[lls[0].second].first];
                int dim2 = ddim[ddoub[lls[0].second].second];
                
                int dist1 = ddoub[lls[0].second].first;
                int dist2 = ddoub[lls[0].second].second;
                
                int ndim = labels.size(); 
                double mp = mod.shared_params[ndim*2 + dim1] / 
                    (mod.shared_params[ndim*2 + dim1] + mod.shared_params[ndim*2 + dim2]);
                
                optimML::multivar_ml_solver prop1(vector<double>{ 1.0 - mod.shared_params[ndim+dim1] }, ll_mix, dll_mix);
                prop1.constrain_01(0);
                optimML::multivar_ml_solver prop2(vector<double>{ 1.0 - mod.shared_params[ndim+dim2] }, ll_mix, dll_mix);
                prop2.constrain_01(0);
                optimML::multivar_ml_solver prop3(vector<double>{ }, ll_mix, dll_mix);
                double f1 = mp*(1.0 - mod.shared_params[ndim+dim1]);
                double f2 = (1.0-mp)*(1.0 - mod.shared_params[ndim+dim2]);
                double f3 = 1.0 - f1 - f2;
                vector<double> pg{ f1, f2, f3};
                prop3.add_param_grp(pg);
                /*
                optimML::brent_solver brent1(ll_mix_new, dll_mix_new);
                optimML::brent_solver brent2(ll_mix_new, dll_mix_new);
                optimML::brent_solver brent3(ll_mix_new, dll_mix_new);
                
                brent1.constrain_01();
                brent2.constrain_01();
                brent3.constrain_01();
                */

                prop1.add_data_fixed("ndim", (int)labels.size());
                prop2.add_data_fixed("ndim", (int)labels.size());
                prop3.add_data_fixed("ndim", (int)labels.size());

                vector<int> prop1_idx1{ dim1 };
                vector<int> prop1_idx2{ -1 };
                vector<int> prop2_idx1{ dim2 };
                vector<int> prop2_idx2{ -1 };
                vector<int> prop3_idx1{ dim1 };
                vector<int> prop3_idx2{ dim2 };
                 
                /*
                vector<double> brent1_mixprop{ -1 };
                vector<double> brent2_mixprop{ -1 };
                vector<double> brent3_mixprop{ mp };
                */
                prop1.add_data("idx1", prop1_idx1);
                prop1.add_data("idx2", prop1_idx2);
                //brent1.add_data("mixprop", brent1_mixprop);
                prop2.add_data("idx1", prop2_idx1);
                prop2.add_data("idx2", prop2_idx2);
                //brent2.add_data("mixprop", brent2_mixprop);
                prop3.add_data("idx1", prop3_idx1);
                prop3.add_data("idx2", prop3_idx2);
                //brent3.add_data("mixprop", brent3_mixprop);
                
                vector<vector<double> > x;
                vector<vector<double> > bgp;

                for (int j = 0; j < labels.size(); ++j){
                    vector<double> x_i { obs[i][j] };
                    x.push_back(x_i);
                    vector<double> bgp_i { mod.shared_params[ndim + j ] };
                    bgp.push_back(bgp_i);
                }

                char buf[50];
                for (int j = 0; j < labels.size(); ++j){
                    sprintf(&buf[0], "n_%d", j);
                    string bufstr = buf;
                    prop1.add_data(bufstr, x[j]);
                    prop2.add_data(bufstr, x[j]);
                    prop3.add_data(bufstr, x[j]);
                    sprintf(&buf[0], "bg_%d", j);
                    bufstr = buf;
                    prop1.add_data(bufstr, bgp[j]);
                    prop2.add_data(bufstr, bgp[j]);
                    prop3.add_data(bufstr, bgp[j]);
                }
                
                prop1.solve();
                prop2.solve();
                prop3.solve();
                
                /*
                double brent1_f = brent1.solve(0,1);
                double brent2_f = brent2.solve(0,1);
                double brent3_f = brent3.solve(0,1);
                
                double brent1_size_ll = dnorm(obs[i][obs[i].size()-1], mod.dists[dist1].params[1][0],
                    mod.dists[dist1].params[1][1]);
                double brent2_size_ll = dnorm(obs[i][obs[i].size()-1], mod.dists[dist2].params[1][0],
                    mod.dists[dist2].params[1][1]);
                double brent3_size_ll = dnorm(obs[i][obs[i].size()-1], mod.dists[lls[0].second].params[1][0],
                    mod.dists[lls[0].second].params[1][1]);
                */
                    fprintf(stderr, "%s\n", bcstr.c_str());
                    fprintf(stderr, "%s: %f %f\n", labels[dim1].c_str(), prop1.results[0], 
                        prop1.log_likelihood);
                    fprintf(stderr, "%s: %f %f\n", labels[dim2].c_str(), prop2.results[0],
                        prop2.log_likelihood);
                    fprintf(stderr, "%s+%s: %f %f %f | %f\n", labels[dim1].c_str(), labels[dim2].c_str(),
                        prop3.results[0], prop3.results[1], prop3.results[2], prop3.log_likelihood);
                if (bcstr == "AAAGGTAAGTAAGGGA"){
                                        
                    /*
                    fprintf(stderr, "%s: %f %f %f\n", labels[dim1].c_str(), brent1_f, brent1.log_likelihood,
                        brent1_size_ll);
                    fprintf(stderr, "%s: %f %f %f\n", labels[dim2].c_str(), brent2_f, brent2.log_likelihood,
                        brent2_size_ll);
                    fprintf(stderr, "%s+%s: %f %f %f\n", labels[dim1].c_str(), labels[dim2].c_str(), brent3_f,
                        brent3.log_likelihood, brent3_size_ll);
                    */
                    //exit(1);
                }
               
                name1 = labels[ddim[ddoub[lls[0].second].first]];
                name2 = labels[ddim[ddoub[lls[0].second].second]];
                if (name2 < name1){
                    string tmp = name1;
                    name1 = name2;
                    name2 = tmp;
                }
            } 
            else if (ddim.count(lls[0].second) > 0){
                int dim1 = ddim[lls[0].second];
                int ndim = labels.size();

                // Solve for mixprop
                optimML::multivar_ml_solver prop1(vector<double>{ 1.0 - mod.shared_params[ndim+dim1] }, ll_mix, dll_mix);
                prop1.constrain_01(0);
                prop1.add_data_fixed("ndim", (int)labels.size());
                
                vector<int> prop1_idx1{ dim1 };
                vector<int> prop1_idx2{ -1 };                
                prop1.add_data("idx1", prop1_idx1);
                prop1.add_data("idx2", prop1_idx2);
                vector<vector<double> > x;
                vector<vector<double> > bgp;

                for (int j = 0; j < labels.size(); ++j){
                    vector<double> x_i { obs[i][j] };
                    x.push_back(x_i);
                    vector<double> bgp_i { mod.shared_params[ndim + j ] };
                    bgp.push_back(bgp_i);
                }

                char buf[50];
                for (int j = 0; j < labels.size(); ++j){
                    sprintf(&buf[0], "n_%d", j);
                    string bufstr = buf;
                    prop1.add_data(bufstr, x[j]);
                    sprintf(&buf[0], "bg_%d", j);
                    bufstr = buf;
                    prop1.add_data(bufstr, bgp[j]);
                }
                
                prop1.solve();
                double count = obs[i][dim1];
                string bcstr = bc2str(bcs[i]);
                fprintf(stdout, "%s\t%s\t%f\t%f\n", bcstr.c_str(), labels[dim1].c_str(), count, 1.0-prop1.results[0]);

                name1 = labels[ddim[lls[0].second]];
            }
            else{
                // Background
                continue;
            }
            assn_llr.emplace(bcs[i], llr);
            assn1.emplace(bcs[i], name1);
            if (name2 != ""){
                assn2.emplace(bcs[i], name2);
            }
        }
        else{
            fprintf(stderr, "LLR == 0!\n");
            string bc_str = bc2str(bcs[i]);
            fprintf(stderr, "%s", bc_str.c_str());
            for (int j = 0; j < obs[i].size(); ++j){
                fprintf(stderr, "\t%.0f", obs[i][j]);
            }
            fprintf(stderr, "\n");
        }
    }
}

double ll_new(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
     
    int ndim = data_i.at("ndim");
    
    int true_idx = data_i.at("true_idx");
    
    double frac_bg = params[true_idx];

    char buf[50];
    string bufstr;
    vector<double> x;
    vector<double> p;
    double psum = 0.0;
    double tot = 0.0;
    
    //double lambda = params[params.size()-1];

    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        int x_i = data_i.at(bufstr);
        double p_0 = params[ndim+i];
        psum += p_0;

        if (i == true_idx){
            p.push_back(frac_bg * p_0 + (1.0 - frac_bg)*1.0);
            //p.push_back(1.0-frac_bg);
        }
        else{
            p.push_back(frac_bg * p_0);
            //p.push_back(frac_bg * (p_0 / (1.0-params[ndim+true_idx])));
        }
        x.push_back((double)x_i);
        tot += x[x.size()-1];
    }
    double lambda = tot*frac_bg;

    return dmultinom(x, p) / log2(exp(1)) - lambda*(1.0 - psum);

    /*
    double ll = 0.0;
    for (int i = 0; i < x.size(); ++i){
        ll += dbinom(tot, x[i], p[i]);
    }
    return ll / log2(exp(1));
    */
}

void dll_new2(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
 
    int ndim = data_i.at("ndim");
    int true_idx = data_i.at("true_idx");
    double frac_bg = params[true_idx];
    char buf[50];
    string bufstr;
    vector<double> x;
    vector<double> p;
    double tot = 0.0;
    double psum = 0.0;

    for (int i = 0; i < ndim; ++i){
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        int x_i = data_i.at(bufstr);
        x.push_back((double)x_i);
        tot += (double)x_i;
        psum += params[ndim+i];
    }
    //double lambda = params[params.size()-1];
    double lambda = tot*frac_bg;
    for (int i = 0; i < ndim; ++i){
        double p0 = params[ndim+i];
        if (i == true_idx){
            results[true_idx] += ((p0 - 1)*x[i])/((p0 - 1)*frac_bg + 1);
            results[ndim+i] += (frac_bg*x[i])/((p0-1)*frac_bg + 1) - lambda;
        }
        else{
            results[true_idx] += x[i] / frac_bg;
            results[ndim+i] += x[i]/p0 - lambda;
        }
    }
    //results[results.size()-1] += 1.0 - psum;
}

void assign_ids_new_2step(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned,
    vector<string>& labels){
    
    vector<vector<int> > counts;
    vector<int> true_idx;
    for (int i = 0; i < labels.size(); ++i){
        vector<int> v;
        counts.push_back(v);
    }
    
    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }
    
    map<int, int> id_tots;
    map<int, int> id_counts;

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
         
        // If we don't have a positive count for the second BC, add a pseudocount
        bool pseudo2 = false;
        // If we don't have a positive count for the third BC, add a pseudocount
        bool pseudo3 = false;

        sort(x->second.begin(), x->second.end());
        if (x->second.size() == 0){
            // No counts
            continue;
        }
        else if (x->second.size() == 1){
            pseudo2 = true;
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
            x->second.push_back(make_pair(-1, ""));
        }
        else if (x->second.size() == 2){
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
        }
        
        int tot = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
        }
        double ll1 = 0;
        double ll2 = 0;
        double ll3 = 0;

        for (int ngroups = 1; ngroups <= 3; ++ngroups){
            double mean_true = 0;
            double tot_true = 0;
            for (int i = 0; i < ngroups; ++i){
                mean_true += -(double)x->second[i].first;
                tot_true++;
            }
            if (tot_true == 0){
                tot_true = 1;
            }
            mean_true /= tot_true;

            double mean_err = 0;
            double tot_err = 0;
            for (int i = ngroups; i < x->second.size(); ++i){
                mean_err += -(double)x->second[i].first;
                tot_err++;
            }
            if (tot_err == 0){
                tot_err = 1;
            }
            mean_err /= tot_err;
            if (mean_err == 0){
                // Prevent issues using Poisson dist
                mean_err = 1;
            }
            double ll = 0;
            for (int i = 0; i < ngroups; ++i){
                ll += dpois(-x->second[i].first, mean_true);
            }
            for (int i = ngroups; i < x->second.size(); ++i){
                ll += dpois(-x->second[i].first, mean_err);
            }
            
            if (ngroups == 1){
                ll1 = ll;
            }
            else if (ngroups == 2){
                ll2 = ll;
            }
            else if (ngroups == 3){
                ll3 = ll;
            }
        }
        
        if (pseudo2){
            tot--;
        }
        if (pseudo3){
            tot--;
        }
        
        if (ll3 > ll1 && ll3 > ll2){
            // 3 or more == noise. Don't make an assignment.        
        }
        else if (ll2 > ll1 && ll2 > ll3){
            if (!pseudo2){
                // Doublet.
            }
            else{
                // Noise
            }
        }
        else if (ll1 > ll2 && ll1 > ll3){
            // Singlet
            bool first = true;
            int counts_idx;
            
            int id_idx = -1;
            int countsum = 0;

            for (vector<pair<int, string> >::iterator y = x->second.begin(); y != x->second.end(); ++y){
                if (first){
                    id_idx = label2idx[y->second];
                    true_idx.push_back(label2idx[y->second]);
                    counts_idx = counts[0].size();
                    for (int z = 0; z < labels.size(); ++z){
                        counts[z].push_back(0);
                    }
                    first = false;
                }
                
                if (y->second != ""){
                    counts[label2idx[y->second]][counts_idx] = -y->first;
                    countsum += -y->first;
                }
            }

            if (id_tots.count(id_idx) == 0){
                id_tots.insert(make_pair(id_idx, 0));
                id_counts.insert(make_pair(id_idx, 0));
            }
            id_tots[id_idx] += countsum;
            id_counts[id_idx]++;
        }
    }
    
    fprintf(stderr, "%ld\n", true_idx.size());
    vector<double> fracs_bg;
    for (int i = 0; i < labels.size(); ++i){
        fprintf(stderr, "counts %ld\n", counts[i].size());
        fracs_bg.push_back(0.05);
    }
    optimML::multivar_ml_solver solver( fracs_bg, ll_new, dll_new2 );
    for (int i = 0; i < labels.size(); ++i){
        solver.constrain_01(i);
    }
    vector<double> grp;
    for (int i = 0; i < labels.size(); ++i){
        grp.push_back(1.0 / (double)labels.size());
    }
    solver.add_param_grp(grp);
    solver.add_data("true_idx", true_idx);
    char buf[50];
    for (int i = 0; i < labels.size(); ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        solver.add_data(bufstr, counts[i]);
    }
    solver.add_data_fixed("ndim", (int)labels.size());
    solver.solve();
    fprintf(stderr, "RESULTS\n");
    for (int i = 0; i < solver.results.size(); ++i){
        string lab;
        if (i < labels.size()){
            fprintf(stderr, "frac_bg ");
            lab = labels[i];
        }
        else{
            lab = labels[i-labels.size()];
        }
        fprintf(stderr, "%s) %f\n", lab.c_str(), solver.results[i]);
    }
    
    map<pair<int, int>, pair<double, double> > allfrac;
    map<pair<int, int>, double> allcount;

    // Make assignments
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        vector<double> countvec;
        double count_tot = 0.0;
        for (int i = 0; i < labels.size(); ++i){
            countvec.push_back(0.0);
        }
        for (vector<pair<int, string> >::iterator y = x->second.begin(); y != x->second.end();
            ++y){
            if (y->second != "" && label2idx.count(y->second) > 0){
                countvec[label2idx[y->second]] = -(double)y->first;
                count_tot += -(double)y->first;
            }
        }
        int top_idx = label2idx[x->second.begin()->second];

        optimML::multivar_ml_solver mixsolver(vector<double>{ 0.5 }, ll_mix, dll_mix);
        mixsolver.constrain_01(0);
        mixsolver.add_data_fixed("idx1", top_idx);
        
        vector<vector<double> > xs;
        vector<vector<double> > bgfracs;

        for (int i = 0; i < labels.size(); ++i){
            vector<double> v1{ solver.results[labels.size() + i] };
            bgfracs.push_back(v1);
            vector<double> v2{ countvec[i] };
            xs.push_back(v2);
        }
        
        char buf[50];
        string bufstr;
        for (int i = 0; i < labels.size(); ++i){
            sprintf(&buf[0], "n_%d", i);
            bufstr = buf;
            mixsolver.add_data(bufstr, xs[i]);
            sprintf(&buf[0], "bg_%d", i);
            bufstr = buf;
            mixsolver.add_data(bufstr, bgfracs[i]);
        }

        mixsolver.add_data_fixed("ndim", (int)labels.size());
        mixsolver.solve();
        
        vector<pair<double, int> > lls;
        
        map<int, pair<double, double> > idx_frac;

        fprintf(stderr, "COUNTS:\n");
        for (int i = 0; i < labels.size(); ++i){
            fprintf(stderr, "%s) %f\n", labels[i].c_str(), countvec[i]);
        }
        fprintf(stderr, "FRAC FG %f\n", mixsolver.results[0]);
        fprintf(stderr, "LL %f\n", mixsolver.log_likelihood);
        double meansize = (double)id_tots[top_idx]/(double)id_counts[top_idx];
        fprintf(stderr, "count %f mean %f LL %f\n", count_tot, meansize,
            dpois(count_tot, meansize)/log2(exp(1)));
        double ll_pois_sing = dpois(count_tot, meansize)/log2(exp(1));
        //lls.push_back(make_pair(-(mixsolver.log_likelihood + ll_pois_sing), -1));
        lls.push_back(make_pair(-mixsolver.log_likelihood, -1));
        idx_frac.insert(make_pair(-1, make_pair(mixsolver.results[0], 0.0)));

        // Test possible doublet combinations
        for (int j = 0; j < labels.size(); ++j){
            if (j != top_idx && countvec[j] > 0.0 && id_counts.count(j) > 0 && id_counts[j] > 0){
                
                optimML::multivar_ml_solver mixsolver_doublet( vector<double>{}, ll_mix, dll_mix);
                vector<double> pg{ 1.0/3.0, 1.0/3.0, 1.0/3.0 };
                mixsolver_doublet.add_param_grp(pg);
                mixsolver_doublet.add_data_fixed("idx1", top_idx);
                mixsolver_doublet.add_data_fixed("idx2", j);
                mixsolver_doublet.add_data_fixed("ndim", (int)labels.size());
                for (int z = 0; z < labels.size(); ++z){
                    sprintf(&buf[0], "n_%d", z);
                    bufstr = buf;
                    mixsolver_doublet.add_data(bufstr, xs[z]);
                    sprintf(&buf[0], "bg_%d", z);
                    bufstr = buf;
                    mixsolver_doublet.add_data(bufstr, bgfracs[z]);
                }
                mixsolver_doublet.solve();
                fprintf(stderr, "combination %s + %s\n", labels[top_idx].c_str(), labels[j].c_str());
                fprintf(stderr, "   log lik %f\n", mixsolver_doublet.log_likelihood);
                fprintf(stderr, "   mix fracs %f %f %f\n", 
                    mixsolver_doublet.results[0], 
                    mixsolver_doublet.results[1],
                    1.0 - mixsolver_doublet.results[0] - mixsolver_doublet.results[1]);
                double meansize_doub = meansize + (double)id_tots[j]/(double)id_counts[j];
                fprintf(stderr, "   count %f mean %f LL %f\n", count_tot, meansize_doub,
                    dpois(count_tot, meansize_doub)/log2(exp(1)));
                double ll_pois_doub = dpois(count_tot, meansize_doub)/log2(exp(1));
                //lls.push_back(make_pair(-(mixsolver_doublet.log_likelihood + ll_pois_doub), j));
                lls.push_back(make_pair(-mixsolver_doublet.log_likelihood, j));
                idx_frac.insert(make_pair(j, make_pair(mixsolver_doublet.results[0], mixsolver_doublet.results[1])));
            }
        }
        
        string bcstr = bc2str(x->first);
        if (bcstr == "AAACCCAGTTCGGTCG"){
            fprintf(stderr, "here\n");
            //exit(0);
        }
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0){
            assn_llr.emplace(x->first, llr);
            pair<int, int> key;
            if (lls[0].second != -1){
                if (top_idx < lls[0].second){
                    key = make_pair(top_idx, lls[0].second);
                }
                else{
                    key = make_pair(lls[0].second, top_idx);
                }
                if (labels[top_idx] < labels[lls[0].second]){
                    assn1.emplace(x->first, labels[top_idx]);
                    assn2.emplace(x->first, labels[lls[0].second]);
                }
                else{
                    assn1.emplace(x->first, labels[lls[0].second]);
                    assn2.emplace(x->first, labels[top_idx]);
                }
            }
            else{
                key = make_pair(top_idx, -1);
                assn1.emplace(x->first, labels[top_idx]);
            }
            if (allfrac.count(key) == 0){
                allfrac.insert(make_pair(key, make_pair(0.0, 0.0)));
                allcount.insert(make_pair(key, 0.0));
            }
            allcount[key]++;
            allfrac[key].first += idx_frac[lls[0].second].first;
            allfrac[key].second += idx_frac[lls[0].second].second;
        }
    }
    
    for (map<pair<int, int>, double>::iterator ac = allcount.begin(); ac != allcount.end(); ++ac){
        double f1 = allfrac[ac->first].first / ac->second;
        double f2 = allfrac[ac->first].second / ac->second;
        double fbg = 1.0 - f1 - f2;
        string name;
        if (ac->first.second == -1){
            name = labels[ac->first.first];
        }
        else{
            if (labels[ac->first.first] < labels[ac->first.second]){
                name = labels[ac->first.first] + "+" + labels[ac->first.second];
            }
            else{
                name = labels[ac->first.second] + "+" + labels[ac->first.first];
            }
        }
        fprintf(stdout, "%s\t%f\t%f\t%f\n", name.c_str(), f1, f2, fbg);
    }
    exit(0);

    return;


    for (int i = 0; i < labels.size(); ++i){
        fprintf(stderr, "SIZE %s) %f\n", labels[i].c_str(), (double)id_tots[i] / (double)id_counts[i]); 
    }

    vector<mixtureDist> dists;
    vector<string> names1;
    vector<string> names2;
    vector<bool> dists_doublet;

    for (int i = 0; i < labels.size(); ++i){
        vector<double> params;
        double frac_bg_i = solver.results[i];
        if (id_counts.count(i) == 0 || id_counts[i] == 0){
            continue;
        }
        for (int x = 0; x < labels.size(); ++x){
            if (x == i){
                params.push_back(frac_bg_i * solver.results[x+labels.size()] + 
                    (1.0 - frac_bg_i)*1.0);
            }
            else{
                params.push_back(frac_bg_i * solver.results[x+labels.size()]);
            }
        }
        
        double meansize = (double)id_tots[i] / (double)id_counts[i];
        mixtureDist dist(vector<string>{ "multinomial", "poisson" }, vector<vector<double> >{ params, vector<double>{meansize} });
        dist.set_num_inputs(0, labels.size());
        dist.name = labels[i];
        dists.push_back(dist);
        names1.push_back(labels[i]);
        names2.push_back("");
        dists_doublet.push_back(false);
        
        /* 
        // Self-self doublets
         mixtureDist dist(vector<string>{ "multinomial", "poisson" }, vector<vector<double> >{ params, vector<double>{2*meansize} });
        dist.set_num_inputs(0, labels.size());
        dist.name = labels[i];
        dists.push_back(dist);
        names1.push_back(labels[i]);
        names2.push_back("");
        dists_doublet.push_back(false);
        */

        for (int j = i + 1; j < labels.size(); ++j){
            if (id_counts.count(j) == 0 || id_counts[j] == 0){
                continue;
            }
            double meansize2 = (double)id_tots[j] / (double)id_counts[j];
            double prop1 = meansize/(meansize+meansize2);

            double frac_bg_j = solver.results[j];
            params.clear();
            for (int x = 0; x < labels.size(); ++x){

                double frac1;
                double frac2;
                
                if (x == i){
                    frac1 = frac_bg_i * solver.results[labels.size()+x] + (1.0-frac_bg_i)*1.0;
                }
                else{
                    frac1 = frac_bg_i * solver.results[labels.size()+x];
                }
                if (x == j){
                    frac2 = frac_bg_j * solver.results[labels.size()+x] + (1.0-frac_bg_j)*1.0;
                }
                else{
                    frac2 = frac_bg_j * solver.results[labels.size()+x];
                }
                
                params.push_back(prop1*frac1 + (1.0-prop1)*frac2);

            }
            mixtureDist distD(vector<string>{ "multinomial", "poisson"}, vector<vector<double> >{ params, vector<double>{ meansize+meansize2 }});
            distD.set_num_inputs(0, labels.size());
            string name;
            string name1;
            string name2;
            if (labels[i] < labels[j]){
                name = labels[i] + "+" + labels[j];
                name1 = labels[i];
                name2 = labels[j];
            }
            else{
                name = labels[j] + "+" + labels[i];
                name1 = labels[j];
                name2 = labels[i];
            }
            distD.name = name;
            dists.push_back(distD);
            names1.push_back(name1);
            names2.push_back(name2);
            dists_doublet.push_back(true);
        }
    }
    
    // Make assignments
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        vector<double> countvec;
        for (int i = 0; i < labels.size(); ++i){
            countvec.push_back(0.0);
        }
        double count_tot = 0.0;
        for (vector<pair<int, string> >::iterator y = x->second.begin(); y != x->second.end();
            ++y){
            if (y->second != "" && label2idx.count(y->second) > 0){
                countvec[label2idx[y->second]] = -(double)y->first;
                count_tot += -(double)y->first;
            }
        }
        // Add data for Poisson
        countvec.push_back(count_tot);

        vector<pair<double, int> > lls;
        for (int i = 0; i < dists.size(); ++i){
            double ll = dists[i].loglik(countvec);
            lls.push_back(make_pair(-ll, i));
        }
        sort(lls.begin(), lls.end());
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0.0){
            int idx = lls[0].second;
            assn_llr.emplace(x->first, llr);
            if (dists_doublet[idx]){
                assn1.emplace(x->first, names1[idx]);
                assn2.emplace(x->first, names2[idx]);
            }
            else{
                assn1.emplace(x->first, names1[idx]);
            }
        }
    }
    return;

    


    map<string, vector<double> > dists_fit;
    map<string, double> weight_tot;
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = 
        bc_ms_counts.begin(); bmc != bc_ms_counts.end(); ++bmc){
        string name;
        double llr;
        if (assn1.count(bmc->first) > 0){
            if (assn2.count(bmc->first) > 0){
                name = assn1[bmc->first] + "+" + assn2[bmc->first];
            }
            else{
                name = assn1[bmc->first];
            }
            llr = assn_llr[bmc->first];
        }
        else{
            name = "Unassigned";
            llr = 1.0;
        }
        if (dists_fit.count(name) == 0){
            vector<double> v;
            dists_fit.insert(make_pair(name, v));
            for (int i = 0; i < label2idx.size(); ++i){
                dists_fit[name].push_back(0.0);
            }
            weight_tot.insert(make_pair(name, 0));
        }
        
        for (vector<pair<int, string> >::iterator c = bmc->second.begin(); c != bmc->second.end(); ++c){
            // skip "pseudocounts" added earlier
            if (label2idx.count(c->second) == 0){
                continue;
            }
            int idx = label2idx[c->second];
            dists_fit[name][idx] += -(double)c->first * llr;
            weight_tot[name] += llr;
        }
    }
    
    vector<mixtureDist> dists_new;

    for (map<string, vector<double> >::iterator df = dists_fit.begin(); df != dists_fit.end(); ++df){
         
        double tot = 0;
        for (int i = 0; i < df->second.size(); ++i){
            if (df->second[i] == 0){
                df->second[i] = 1.0;
            }
            tot += df->second[i];        
        }
        for (int i = 0; i < df->second.size(); ++i){
            df->second[i] /= tot;
        }
        
        mixtureDist dist("multinomial", df->second);
        dist.set_num_inputs(label2idx.size());
        dist.name = df->first;
        dists_new.push_back(dist);
        dists_new[dists_new.size()-1].print(0);
    }
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = bc_ms_counts.begin();
        bmc != bc_ms_counts.end(); ++bmc){
        
        vector<double> obs;
        for (int i = 0; i < label2idx.size(); ++i){
            obs.push_back(0.0);
        }
        for (vector<pair<int, string> >::iterator bmc2 = bmc->second.begin(); bmc2 != 
            bmc->second.end(); ++bmc2){
            if (label2idx.count(bmc2->second) > 0){
                obs[label2idx[bmc2->second]] = -(double)bmc2->first;
            }
        }
        
        vector<pair<double, string> > lls;
        for (int i = 0; i < dists_new.size(); ++i){
            double ll = dists_new[i].loglik(obs);
            lls.push_back(make_pair(-ll, dists_new[i].name));
        }
        sort(lls.begin(), lls.end());
        string assn = lls[0].second;
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0){
            string bc_str = bc2str(bmc->first);
            string name_old;
            if (assn1.count(bmc->first) > 0){
                if (assn2.count(bmc->first) > 0){
                    name_old = assn1[bmc->first] + "+" + assn2[bmc->first];
                }
                else{
                    name_old = assn1[bmc->first];
                }
            }
            else{
                name_old = "Unassigned";
            }
            fprintf(stdout, "%s\t%s\t%f\t%s\t%f\n", bc_str.c_str(), name_old.c_str(),
                assn_llr[bmc->first], assn.c_str(), llr);
        }
    }

    exit(0);


    // Get distributions for off-target
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn1.begin(); a != assn1.end(); ++a){
        int tot = 0;
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                tot += -c->first;
            }
        }
        
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                fprintf(stdout, "%s\t%f\n", c->second.c_str(), (double)-c->first / (double)tot);        
            }
        }
    }
    exit(0);
    /* 
    return; 
    sum1 /= tot1;
    sum2 /= tot2;
    sum3 /= tot3;

    fprintf(stderr, "1 %f 2 %f 3 %f\n", sum1, sum2, sum3);

    double med_noise = median(noise);
    double l_noise = (med_noise - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "noise %f\n", l_noise);
    double med_valid = median(valid);
    double l_valid = (med_valid - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "valid %f\n", l_valid); 
    double med_doub = median(doub);
    double l_doub = (med_doub - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "doub %f\n", l_doub); 
    exit(0);
    // Optional second pass through data: find mean successful count and mean unsuccessful
    // count; throw out identifications where successful count is more likely under the
    // data set-wide unsuccessful count distribution.

    // This seems unnecessary - now, LLRs reflect what we can say about the observed data; 
    // we don't need to care what other cells look like.
    
    bool second_pass = false;

    if (second_pass){
        double mean_valid = 0;
        double mean_noise = 0;

        double valid_weight = 1.0/(double)valid.size();
        for (int i = 0; i < valid.size(); ++i){
            mean_valid += valid_weight*(double)valid[i];    
        }
        double noise_weight = 1.0/(double)noise.size();
        for (int i = 0; i < noise.size(); ++i){
            mean_noise += noise_weight*(double)noise[i];
        }
        
        for (robin_hood::unordered_map<unsigned long, double>::iterator a = assn_llr.begin();
            a != assn_llr.end();){
            bool rm = false;
            double llr2 = 0.0;
            if (assn2.count(a->first) > 0){
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) -
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
                llr2 += dpois(-bc_ms_counts[a->first][1].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][1].first, mean_noise);
            }
            else{
                int count1 = -bc_ms_counts[a->first][0].first;
                int count2 = -bc_ms_counts[a->first][1].first;
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
            }
            
            double llr_new = a->second + llr2;
            if (llr_new < 0){
                // Eliminate
                assn1.erase(a->first);
                if (assn2.count(a->first) > 0){
                    assn2.erase(a->first);
                }
                rm = true;
            }
            else{
                a->second = llr_new;
            }

            if (rm){
                assn_llr.erase(a++);
            }
            else{
                a++;
            }
        }
    }
    */
}
void assign_ids_new(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned,
    vector<string>& labels){

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
         
        // If we don't have a positive count for the second BC, add a pseudocount
        bool pseudo2 = false;
        // If we don't have a positive count for the third BC, add a pseudocount
        bool pseudo3 = false;

        sort(x->second.begin(), x->second.end());
        if (x->second.size() == 0){
            // No counts
            continue;
        }
        else if (x->second.size() == 1){
            pseudo2 = true;
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
            x->second.push_back(make_pair(-1, ""));
        }
        else if (x->second.size() == 2){
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
        }
        
        int tot = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
        }
        double ll1 = 0;
        double ll2 = 0;
        double ll3 = 0;

        for (int ngroups = 1; ngroups <= 3; ++ngroups){
            double mean_true = 0;
            double tot_true = 0;
            for (int i = 0; i < ngroups; ++i){
                mean_true += -(double)x->second[i].first;
                tot_true++;
            }
            if (tot_true == 0){
                tot_true = 1;
            }
            mean_true /= tot_true;

            double mean_err = 0;
            double tot_err = 0;
            for (int i = ngroups; i < x->second.size(); ++i){
                mean_err += -(double)x->second[i].first;
                tot_err++;
            }
            if (tot_err == 0){
                tot_err = 1;
            }
            mean_err /= tot_err;
            if (mean_err == 0){
                // Prevent issues using Poisson dist
                mean_err = 1;
            }
            double ll = 0;
            for (int i = 0; i < ngroups; ++i){
                ll += dpois(-x->second[i].first, mean_true);
            }
            for (int i = ngroups; i < x->second.size(); ++i){
                ll += dpois(-x->second[i].first, mean_err);
            }
            
            if (ngroups == 1){
                ll1 = ll;
            }
            else if (ngroups == 2){
                ll2 = ll;
            }
            else if (ngroups == 3){
                ll3 = ll;
            }
        }
        
        if (pseudo2){
            tot--;
        }
        if (pseudo3){
            tot--;
        }
        
        double llr;
        if (ll3 > ll1 && ll3 > ll2){
            // 3 or more == noise. Don't make an assignment.        
        }
        else if (ll2 > ll1 && ll2 > ll3){
            if (!pseudo2){
                // Doublet.
                if (ll1 > ll3){
                    llr = ll2 - ll1;
                }
                else{
                    llr = ll2 - ll3;
                }
                if (x->second[0].second < x->second[1].second){
                    assn1.emplace(x->first, x->second[0].second);
                    assn2.emplace(x->first, x->second[1].second);
                }
                else{
                    assn1.emplace(x->first, x->second[1].second);
                    assn2.emplace(x->first, x->second[0].second);
                }
                assn_llr.emplace(x->first, llr);
            }
            else{
                // Noise
            }
        }
        else if (ll1 > ll2 && ll1 > ll3){
            llr = ll1 - ll2;
            assn1.emplace(x->first, x->second[0].second);
            assn_llr.emplace(x->first, llr);
        }
        else{
        }
    }

    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }

    map<string, vector<double> > dists_fit;
    map<string, double> weight_tot;
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = 
        bc_ms_counts.begin(); bmc != bc_ms_counts.end(); ++bmc){
        string name;
        double llr;
        if (assn1.count(bmc->first) > 0){
            if (assn2.count(bmc->first) > 0){
                name = assn1[bmc->first] + "+" + assn2[bmc->first];
            }
            else{
                name = assn1[bmc->first];
            }
            llr = assn_llr[bmc->first];
        }
        else{
            name = "Unassigned";
            llr = 1.0;
        }
        if (dists_fit.count(name) == 0){
            vector<double> v;
            dists_fit.insert(make_pair(name, v));
            for (int i = 0; i < label2idx.size(); ++i){
                dists_fit[name].push_back(0.0);
            }
            weight_tot.insert(make_pair(name, 0));
        }
        
        for (vector<pair<int, string> >::iterator c = bmc->second.begin(); c != bmc->second.end(); ++c){
            // skip "pseudocounts" added earlier
            if (label2idx.count(c->second) == 0){
                continue;
            }
            int idx = label2idx[c->second];
            dists_fit[name][idx] += -(double)c->first * llr;
            weight_tot[name] += llr;
        }
    }
    
    vector<mixtureDist> dists_new;

    for (map<string, vector<double> >::iterator df = dists_fit.begin(); df != dists_fit.end(); ++df){
         
        double tot = 0;
        for (int i = 0; i < df->second.size(); ++i){
            if (df->second[i] == 0){
                df->second[i] = 1.0;
            }
            tot += df->second[i];        
        }
        for (int i = 0; i < df->second.size(); ++i){
            df->second[i] /= tot;
        }
        
        mixtureDist dist("multinomial", df->second);
        dist.set_num_inputs(label2idx.size());
        dist.name = df->first;
        dists_new.push_back(dist);
        dists_new[dists_new.size()-1].print(0);
    }
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = bc_ms_counts.begin();
        bmc != bc_ms_counts.end(); ++bmc){
        
        vector<double> obs;
        for (int i = 0; i < label2idx.size(); ++i){
            obs.push_back(0.0);
        }
        for (vector<pair<int, string> >::iterator bmc2 = bmc->second.begin(); bmc2 != 
            bmc->second.end(); ++bmc2){
            if (label2idx.count(bmc2->second) > 0){
                obs[label2idx[bmc2->second]] = -(double)bmc2->first;
            }
        }
        
        vector<pair<double, string> > lls;
        for (int i = 0; i < dists_new.size(); ++i){
            double ll = dists_new[i].loglik(obs);
            lls.push_back(make_pair(-ll, dists_new[i].name));
        }
        sort(lls.begin(), lls.end());
        string assn = lls[0].second;
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0){
            string bc_str = bc2str(bmc->first);
            string name_old;
            if (assn1.count(bmc->first) > 0){
                if (assn2.count(bmc->first) > 0){
                    name_old = assn1[bmc->first] + "+" + assn2[bmc->first];
                }
                else{
                    name_old = assn1[bmc->first];
                }
            }
            else{
                name_old = "Unassigned";
            }
            fprintf(stdout, "%s\t%s\t%f\t%s\t%f\n", bc_str.c_str(), name_old.c_str(),
                assn_llr[bmc->first], assn.c_str(), llr);
        }
    }

    exit(0);


    // Get distributions for off-target
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn1.begin(); a != assn1.end(); ++a){
        int tot = 0;
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                tot += -c->first;
            }
        }
        
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                fprintf(stdout, "%s\t%f\n", c->second.c_str(), (double)-c->first / (double)tot);        
            }
        }
    }
    exit(0);
    /* 
    return; 
    sum1 /= tot1;
    sum2 /= tot2;
    sum3 /= tot3;

    fprintf(stderr, "1 %f 2 %f 3 %f\n", sum1, sum2, sum3);

    double med_noise = median(noise);
    double l_noise = (med_noise - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "noise %f\n", l_noise);
    double med_valid = median(valid);
    double l_valid = (med_valid - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "valid %f\n", l_valid); 
    double med_doub = median(doub);
    double l_doub = (med_doub - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "doub %f\n", l_doub); 
    exit(0);
    // Optional second pass through data: find mean successful count and mean unsuccessful
    // count; throw out identifications where successful count is more likely under the
    // data set-wide unsuccessful count distribution.

    // This seems unnecessary - now, LLRs reflect what we can say about the observed data; 
    // we don't need to care what other cells look like.
    
    bool second_pass = false;

    if (second_pass){
        double mean_valid = 0;
        double mean_noise = 0;

        double valid_weight = 1.0/(double)valid.size();
        for (int i = 0; i < valid.size(); ++i){
            mean_valid += valid_weight*(double)valid[i];    
        }
        double noise_weight = 1.0/(double)noise.size();
        for (int i = 0; i < noise.size(); ++i){
            mean_noise += noise_weight*(double)noise[i];
        }
        
        for (robin_hood::unordered_map<unsigned long, double>::iterator a = assn_llr.begin();
            a != assn_llr.end();){
            bool rm = false;
            double llr2 = 0.0;
            if (assn2.count(a->first) > 0){
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) -
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
                llr2 += dpois(-bc_ms_counts[a->first][1].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][1].first, mean_noise);
            }
            else{
                int count1 = -bc_ms_counts[a->first][0].first;
                int count2 = -bc_ms_counts[a->first][1].first;
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
            }
            
            double llr_new = a->second + llr2;
            if (llr_new < 0){
                // Eliminate
                assn1.erase(a->first);
                if (assn2.count(a->first) > 0){
                    assn2.erase(a->first);
                }
                rm = true;
            }
            else{
                a->second = llr_new;
            }

            if (rm){
                assn_llr.erase(a++);
            }
            else{
                a++;
            }
        }
    }
    */
}

void assign_ids2(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned,
    vector<string>& labels){

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
         
        // If we don't have a positive count for the second BC, add a pseudocount
        bool pseudo2 = false;
        // If we don't have a positive count for the third BC, add a pseudocount
        bool pseudo3 = false;

        sort(x->second.begin(), x->second.end());
        if (x->second.size() == 0){
            // No counts
            continue;
        }
        else if (x->second.size() == 1){
            pseudo2 = true;
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
            x->second.push_back(make_pair(-1, ""));
        }
        else if (x->second.size() == 2){
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
        }
        
        int tot = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
        }
        double ll1 = 0;
        double ll2 = 0;
        double ll3 = 0;

        for (int ngroups = 1; ngroups <= 3; ++ngroups){
            double mean_true = 0;
            double tot_true = 0;
            for (int i = 0; i < ngroups; ++i){
                mean_true += -(double)x->second[i].first;
                tot_true++;
            }
            if (tot_true == 0){
                tot_true = 1;
            }
            mean_true /= tot_true;

            double mean_err = 0;
            double tot_err = 0;
            for (int i = ngroups; i < x->second.size(); ++i){
                mean_err += -(double)x->second[i].first;
                tot_err++;
            }
            if (tot_err == 0){
                tot_err = 1;
            }
            mean_err /= tot_err;
            if (mean_err == 0){
                // Prevent issues using Poisson dist
                mean_err = 1;
            }
            double ll = 0;
            for (int i = 0; i < ngroups; ++i){
                ll += dpois(-x->second[i].first, mean_true);
            }
            for (int i = ngroups; i < x->second.size(); ++i){
                ll += dpois(-x->second[i].first, mean_err);
            }
            
            if (ngroups == 1){
                ll1 = ll;
            }
            else if (ngroups == 2){
                ll2 = ll;
            }
            else if (ngroups == 3){
                ll3 = ll;
            }
        }
        
        if (pseudo2){
            tot--;
        }
        if (pseudo3){
            tot--;
        }
        
        double llr;
        if (ll3 > ll1 && ll3 > ll2){
            // 3 or more == noise. Don't make an assignment.        
        }
        else if (ll2 > ll1 && ll2 > ll3){
            if (!pseudo2){
                // Doublet.
                if (ll1 > ll3){
                    llr = ll2 - ll1;
                }
                else{
                    llr = ll2 - ll3;
                }
                if (x->second[0].second < x->second[1].second){
                    assn1.emplace(x->first, x->second[0].second);
                    assn2.emplace(x->first, x->second[1].second);
                }
                else{
                    assn1.emplace(x->first, x->second[1].second);
                    assn2.emplace(x->first, x->second[0].second);
                }
                assn_llr.emplace(x->first, llr);
            }
            else{
                // Noise
            }
        }
        else if (ll1 > ll2 && ll1 > ll3){
            llr = ll1 - ll2;
            assn1.emplace(x->first, x->second[0].second);
            assn_llr.emplace(x->first, llr);
        }
        else{
        }
    }

    map<string, int> label2idx;
    for (int i = 0; i < labels.size(); ++i){
        label2idx.insert(make_pair(labels[i], i));
    }

    map<string, vector<double> > dists_fit;
    map<string, double> weight_tot;
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = 
        bc_ms_counts.begin(); bmc != bc_ms_counts.end(); ++bmc){
        string name;
        double llr;
        if (assn1.count(bmc->first) > 0){
            if (assn2.count(bmc->first) > 0){
                name = assn1[bmc->first] + "+" + assn2[bmc->first];
            }
            else{
                name = assn1[bmc->first];
            }
            llr = assn_llr[bmc->first];
        }
        else{
            name = "Unassigned";
            llr = 1.0;
        }
        if (dists_fit.count(name) == 0){
            vector<double> v;
            dists_fit.insert(make_pair(name, v));
            for (int i = 0; i < label2idx.size(); ++i){
                dists_fit[name].push_back(0.0);
            }
            weight_tot.insert(make_pair(name, 0));
        }
        
        for (vector<pair<int, string> >::iterator c = bmc->second.begin(); c != bmc->second.end(); ++c){
            // skip "pseudocounts" added earlier
            if (label2idx.count(c->second) == 0){
                continue;
            }
            int idx = label2idx[c->second];
            dists_fit[name][idx] += -(double)c->first * llr;
            weight_tot[name] += llr;
        }
    }
    
    vector<mixtureDist> dists_new;

    for (map<string, vector<double> >::iterator df = dists_fit.begin(); df != dists_fit.end(); ++df){
         
        double tot = 0;
        for (int i = 0; i < df->second.size(); ++i){
            if (df->second[i] == 0){
                df->second[i] = 1.0;
            }
            tot += df->second[i];        
        }
        for (int i = 0; i < df->second.size(); ++i){
            df->second[i] /= tot;
        }
        
        mixtureDist dist("multinomial", df->second);
        dist.set_num_inputs(label2idx.size());
        dist.name = df->first;
        dists_new.push_back(dist);
        dists_new[dists_new.size()-1].print(0);
    }
    
    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator bmc = bc_ms_counts.begin();
        bmc != bc_ms_counts.end(); ++bmc){
        
        vector<double> obs;
        for (int i = 0; i < label2idx.size(); ++i){
            obs.push_back(0.0);
        }
        for (vector<pair<int, string> >::iterator bmc2 = bmc->second.begin(); bmc2 != 
            bmc->second.end(); ++bmc2){
            if (label2idx.count(bmc2->second) > 0){
                obs[label2idx[bmc2->second]] = -(double)bmc2->first;
            }
        }
        
        vector<pair<double, string> > lls;
        for (int i = 0; i < dists_new.size(); ++i){
            double ll = dists_new[i].loglik(obs);
            lls.push_back(make_pair(-ll, dists_new[i].name));
        }
        sort(lls.begin(), lls.end());
        string assn = lls[0].second;
        double llr = -lls[0].first - -lls[1].first;
        if (llr > 0){
            string bc_str = bc2str(bmc->first);
            string name_old;
            if (assn1.count(bmc->first) > 0){
                if (assn2.count(bmc->first) > 0){
                    name_old = assn1[bmc->first] + "+" + assn2[bmc->first];
                }
                else{
                    name_old = assn1[bmc->first];
                }
            }
            else{
                name_old = "Unassigned";
            }
            fprintf(stdout, "%s\t%s\t%f\t%s\t%f\n", bc_str.c_str(), name_old.c_str(),
                assn_llr[bmc->first], assn.c_str(), llr);
        }
    }

    exit(0);


    // Get distributions for off-target
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn1.begin(); a != assn1.end(); ++a){
        int tot = 0;
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                tot += -c->first;
            }
        }
        
        for (vector<pair<int, string> >::iterator c = bc_ms_counts[a->first].begin(); c != 
            bc_ms_counts[a->first].end(); ++c){
            if (c->second != a->second){
                fprintf(stdout, "%s\t%f\n", c->second.c_str(), (double)-c->first / (double)tot);        
            }
        }
    }
    exit(0);
    /* 
    return; 
    sum1 /= tot1;
    sum2 /= tot2;
    sum3 /= tot3;

    fprintf(stderr, "1 %f 2 %f 3 %f\n", sum1, sum2, sum3);

    double med_noise = median(noise);
    double l_noise = (med_noise - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "noise %f\n", l_noise);
    double med_valid = median(valid);
    double l_valid = (med_valid - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "valid %f\n", l_valid); 
    double med_doub = median(doub);
    double l_doub = (med_doub - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "doub %f\n", l_doub); 
    exit(0);
    // Optional second pass through data: find mean successful count and mean unsuccessful
    // count; throw out identifications where successful count is more likely under the
    // data set-wide unsuccessful count distribution.

    // This seems unnecessary - now, LLRs reflect what we can say about the observed data; 
    // we don't need to care what other cells look like.
    
    bool second_pass = false;

    if (second_pass){
        double mean_valid = 0;
        double mean_noise = 0;

        double valid_weight = 1.0/(double)valid.size();
        for (int i = 0; i < valid.size(); ++i){
            mean_valid += valid_weight*(double)valid[i];    
        }
        double noise_weight = 1.0/(double)noise.size();
        for (int i = 0; i < noise.size(); ++i){
            mean_noise += noise_weight*(double)noise[i];
        }
        
        for (robin_hood::unordered_map<unsigned long, double>::iterator a = assn_llr.begin();
            a != assn_llr.end();){
            bool rm = false;
            double llr2 = 0.0;
            if (assn2.count(a->first) > 0){
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) -
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
                llr2 += dpois(-bc_ms_counts[a->first][1].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][1].first, mean_noise);
            }
            else{
                int count1 = -bc_ms_counts[a->first][0].first;
                int count2 = -bc_ms_counts[a->first][1].first;
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
            }
            
            double llr_new = a->second + llr2;
            if (llr_new < 0){
                // Eliminate
                assn1.erase(a->first);
                if (assn2.count(a->first) > 0){
                    assn2.erase(a->first);
                }
                rm = true;
            }
            else{
                a->second = llr_new;
            }

            if (rm){
                assn_llr.erase(a++);
            }
            else{
                a++;
            }
        }
    }
    */
}

void write_assignments(string filename, 
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    FILE* f_out = fopen(filename.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn1.begin(); 
        a != assn1.end(); ++a){
        char type = 'S';
        string name = a->second;
        if (assn2.count(a->first) > 0){
            type = 'D';
            name = a->second + "+" + assn2[a->first];
        }
        string bc_str = bc2str(a->first);
        fprintf(f_out, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type, assn_llr[a->first]); 
    }
    fclose(f_out);
}

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

double ll_multinom_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double xsum = 1;
    double term2 = 0;
    double term3 = 0;
    
    vector<string> keys = { "n1", "n2", "n3", "n4", "n5", "n6" };
    for (int i = 0; i < keys.size(); ++i){
        double x = (double)data_i.at(keys[i]);
        xsum += x;
        term2 += lgammaf(x + 1);
        term3 += x * log(params[i]);
    }
    double term1 = lgammaf(xsum);
    
    /*    
    if (term1 - term2 + term3 > 0){
        fprintf(stderr, "term1 %f term2 %f term3 %f\n", term1, term2, term3);
        fprintf(stderr, "%f\n", term1 - term2 + term3);
        for (int i = 0; i < keys.size(); ++i){
            fprintf(stderr, "%d %f\n", data_i.at(keys[i]), params[i]);
        }
        exit(1);
    }
    */
    return (term1 - term2 + term3); 
    vector<double> n;
    n.push_back(data_i.at("n1"));
    n.push_back(data_i.at("n2"));
    n.push_back(data_i.at("n3"));
    n.push_back(data_i.at("n4"));
    n.push_back(data_i.at("n5"));
    n.push_back(data_i.at("n6"));
    vector<double> p = params;
    return dmultinom(n, p)/log2(exp(1));
}

void dll_multinom_test(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    results[0] = (double)data_i.at("n1")/params[0];
    results[1] = (double)data_i.at("n2")/params[1];
    results[2] = (double)data_i.at("n3")/params[2];
    results[3] = (double)data_i.at("n4")/params[3];
    results[4] = (double)data_i.at("n5")/params[4];
    results[5] = (double)data_i.at("n6")/params[5];
    double tot = (double)(data_i.at("n1") + data_i.at("n2") + data_i.at("n3") + data_i.at("n4") + data_i.at("n5") + data_i.at("n6"));
    results[0] -= tot;
    results[1] -= tot;
    results[2] -= tot;
    results[3] -= tot;
    results[4] -= tot;
    results[5] -= tot;
}

int main(int argc, char *argv[]) {    
    /* 
    vector<double> paramgrp{ 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };
    //optimML::multivar_ml_solver solver( vector<double>{ }, ll_multinom_test, dll_multinom_test);
    optimML::multivar_ml_solver solver( vector<double>{ 0.01 }, ll_new, dll_new2); 
    solver.constrain_01(0);
    solver.add_param_grp(paramgrp);
    //solver.add_one_param(500.0);

    ifstream inf("multinom_test2.txt");
    int num1;
    int num2;
    int num3;
    int num4;
    int num5;
    int num6;
    
    vector<int> nums1;
    vector<int> nums2;
    vector<int> nums3;
    vector<int> nums4;
    vector<int> nums5;
    vector<int> nums6;
    
    vector<int> trueidx;

    while (inf >> num1 >> num2 >> num3 >> num4 >> num5 >> num6){
        nums1.push_back(num1);
        nums2.push_back(num2);
        nums3.push_back(num3);
        nums4.push_back(num4);
        nums5.push_back(num5);
        nums6.push_back(num6);
        
        int max = num1;
        if (num2 > max){
            max = num2;
        }
        if (num3 > max){
            max = num3;
        }
        if (num4 > max){
            max = num4;
        }
        if (num5 > max){
            max = num5;
        }
        if (num6 > max){
            max = num6;
        }    

        if (max == num1){
            trueidx.push_back(0);
        }
        else if (max == num2){
            trueidx.push_back(1);
        }
        else if (max == num3){
            trueidx.push_back(2);
        }
        else if (max == num4){
            trueidx.push_back(3);
        }
        else if (max == num5){
            trueidx.push_back(4);
        }
        else if (max == num6){
            trueidx.push_back(5);
        }

    } 
    
    solver.add_data_fixed("ndim", 6);
    solver.add_data("true_idx", trueidx);
    solver.add_data("n_0", nums1);
    solver.add_data("n_1", nums2);
    solver.add_data("n_2", nums3);
    solver.add_data("n_3", nums4);
    solver.add_data("n_4", nums5);
    solver.add_data("n_5", nums6);
    fprintf(stderr, "solving\n");
    solver.solve();
    fprintf(stderr, "RESULTS %f %f %f %f %f %f %f\n", solver.results[0], solver.results[1], 
        solver.results[2], solver.results[3], solver.results[4], solver.results[5], 
        solver.results[6]);
    //fprintf(stderr, "lambda %f\n", solver.results[7]);
    exit(1);
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
    int mismatches = -1;
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
    robin_hood::unordered_map<unsigned long, string> assn1;
    robin_hood::unordered_map<unsigned long, string> assn2;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    //assign_ids_new2(bc_ms_counts, labels, assn1, assn2, assn_llr, output_unassigned);
    //assign_ids_new_2step(bc_ms_counts, assn1, assn2, assn_llr, output_unassigned, labels);  
    
    assign_ids_brandnew(bc_ms_counts, assn1, assn2, assn_llr, output_unassigned, labels);
    //assign_ids_yetagain(bc_ms_counts, assn1, assn2, assn_llr, output_unassigned, labels);


    string assnfilename = output_prefix + ".assignments";
    write_assignments(assnfilename, assn1, assn2, assn_llr);

    return 0;  
}
