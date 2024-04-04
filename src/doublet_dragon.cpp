#include <getopt.h>
#include <argp.h>
#include <string>
#include <algorithm>
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
#include <optimML/multivar_ml.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "doublet_dragon [FILES]\n");
    fprintf(stderr, "Given multiple .assignments files from the same library, estimates\n");
    fprintf(stderr, "the likeliest underlying doublet rate (via maximum likelihood).\n");
    fprintf(stderr, "Also infers the frequency of each individual type in the data set.\n");
    fprintf(stderr, "\nEXAMPLE:\n");
    fprintf(stderr, "If you pooled data from three species, separated these data using\n");
    fprintf(stderr, "demux_species, inferred mitochondrial haplotypes via demux_mt,\n");
    fprintf(stderr, "demultiplexed individuals from a VCF using demux_vcf, and processed\n");
    fprintf(stderr, "MULTIseq data on the same library using demux_multiseq, then you can\n");
    fprintf(stderr, "run doublet_dragon on all these output files to estimate the underlying\n");
    fprintf(stderr, "global doublet rate.\n");
    exit(code);
}

/**
 * Returns total number of cells
 */
int parse_assignments(string filename, 
    map<string, int>& counts_singlet, 
    map<string, int>& counts_doublet,
    map<string, double>& weights){
    
    ifstream inf(filename.c_str());
    string bc_str;
    string id;
    char type;
    double llr;
    
    map<string, int> id_count;
    
    int ntot = 0;

    while (inf >> bc_str >> id >> type >> llr){
        ntot++;
        if (id_count.count(id) == 0){
            id_count.insert(make_pair(id, 0));
            weights.insert(make_pair(id, 0.0));
        }    
        
        id_count[id]++;
        weights[id] += llr;

        if (type == 'S'){
            if (counts_singlet.count(id) == 0){
                counts_singlet.insert(make_pair(id, 0));
            }
            counts_singlet[id]++;
        }
        else{
            if (counts_doublet.count(id) == 0){
                counts_doublet.insert(make_pair(id, 0));
            }
            counts_doublet[id]++;
        }
    }

    double wsum = 0.0;
    for (map<string, double>::iterator w = weights.begin(); w != weights.end(); ++w){
        w->second = w->second / (double)id_count[w->first];
        wsum += w->second;
    }
    for (map<string, double>::iterator w = weights.begin(); w != weights.end(); ++w){
        w->second /= wsum;
    }

    return ntot;
}

double loglik(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double d = params[0];
    double n = (double)data_i.at("n");
    double k = (double)data_i.at("k");
    int p_idx = data_i.at("p_idx");
    int p_idx2 = data_i.at("p_idx2");
    int bias_idx = data_i.at("bias_idx");
    double dx = params[bias_idx];

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d*dx)*params[p_idx] + d*dx*params[p_idx] * params[p_idx];
    }
    else{
        // Doublet
        p = d*dx*params[p_idx]*params[p_idx2];
    }
    return logbinom(n, k, p);
}

void dloglik_dx(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double d = params[0];
    double n = (double)data_i.at("n");
    double k = (double)data_i.at("k");
    int p_idx = data_i.at("p_idx");
    int p_idx2 = data_i.at("p_idx2");
    int bias_idx = data_i.at("bias_idx");
    double dx = params[bias_idx];

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d*dx)*params[p_idx] + d*dx*params[p_idx] * params[p_idx];
    }
    else{
        // Doublet
        p = d*dx*params[p_idx]*params[p_idx2];
    }

    double dy_dp = (k - n*p)/(p - p*p);

    if (p_idx2 == -1){
        double a = params[p_idx];
        double a2 = a*a;
        results[0] += dy_dp * (-dx*a + dx*a2);
        //results[bias_idx] += dy_dp * (-d*a + d*a2);
        results[bias_idx] = 0;
        results[p_idx] += dy_dp * (1.0 - d*dx + 2*d*dx*a);
    }
    else{
        double a = params[p_idx];
        double b = params[p_idx2];
        results[0] += dy_dp * (dx*a*b);
        //results[bias_idx] += dy_dp * (d*a*b);
        results[bias_idx] = 0;
        results[p_idx] += dy_dp * (d*dx*b);
        results[p_idx2] += dy_dp * (d*dx*a);
    }
}

int main(int argc, char *argv[]) {    
    
    if (argc < 2){
        help(1);
    } 
    
    vector<int> n;
    vector<int> k;
    vector<int> bias_idx;
    vector<int> p_idx;
    vector<int> p_idx2;
    vector<double> weights;

    vector<double> params;
    params.push_back(0.01); // doublet rate
    
    vector<int> constr_01 = {0};
    vector<int> constr_pos;
     
    for (int i = 1; i < argc; ++i){
        map<string, int> counts_singlet;
        map<string, int> counts_doublet;
        map<string, double> weights_this;
        int totcells = parse_assignments(argv[i], counts_singlet, counts_doublet, weights_this);
        if (counts_doublet.size() == 0){
            fprintf(stderr, "ERROR: no doublet identifications in %s\n", argv[i]);
            exit(1);
        }
        int bias_idx_this = params.size();
        constr_pos.push_back(bias_idx_this);    
        // doublet rate bias (multiplier)
        params.push_back(1.0);
        
        int tot_singlet = 0;
        for (map<string, int>::iterator cs = counts_singlet.begin();
            cs != counts_singlet.end(); ++cs){
            tot_singlet += cs->second;
        }

        map<string, int> p2idx;
        for (map<string, int>::iterator cs = counts_singlet.begin(); 
            cs != counts_singlet.end(); ++cs){
            int p_idx_this = params.size();
            constr_01.push_back(p_idx_this);
            params.push_back((double)cs->second / (double)tot_singlet);
            p2idx.insert(make_pair(cs->first, p_idx_this));
            n.push_back(totcells);
            k.push_back(cs->second);
            bias_idx.push_back(bias_idx_this);
            p_idx.push_back(p_idx_this);
            p_idx2.push_back(-1);
            weights.push_back(weights_this[cs->first]);
        }
        for (map<string, int>::iterator cd = counts_doublet.begin();
            cd != counts_doublet.end(); ++cd){
            size_t pos = cd->first.find("+");
            if (pos == string::npos){
                pos = cd->first.find("x");
            }
            if (pos == string::npos){
                pos = cd->first.find("X");
            }
            if (pos == string::npos){
                fprintf(stderr, "ERROR: unable to identify components of doublet ID %s\n", cd->first.c_str());
                exit(1);
            }
            string id1 = cd->first.substr(0, pos);
            string id2 = cd->first.substr(pos+1, cd->first.length()-pos-1);
            if (p2idx.count(id1) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s in %s\n", id1.c_str(), argv[i]);
                exit(1);
            }
            if (p2idx.count(id2) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s in %s\n", id2.c_str(), argv[i]);
                exit(1);
            }
            n.push_back(totcells);
            k.push_back(cd->second);
            bias_idx.push_back(bias_idx_this);
            p_idx.push_back(p2idx[id1]);
            p_idx2.push_back(p2idx[id2]);
            weights.push_back(weights_this[cd->first]);
        }
    }
    
    optimML::multivar_ml_solver solver(params, loglik, dloglik_dx);
    for (vector<int>::iterator cp = constr_pos.begin(); cp != constr_pos.end(); ++cp){
        solver.constrain_pos(*cp);
    }
    for (vector<int>::iterator c01 = constr_01.begin(); c01 != constr_01.end(); ++c01){
        solver.constrain_01(*c01);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("bias_idx", bias_idx);
    solver.add_data("p_idx", p_idx);
    solver.add_data("p_idx2", p_idx2);
    solver.add_weights(weights);
    solver.solve();

    fprintf(stderr, "here\n");
    for (int i = 0; i < solver.results.size(); ++i){
        fprintf(stderr, "%d) %f\n", i, solver.results[i]);
    }
    return 0;
}
