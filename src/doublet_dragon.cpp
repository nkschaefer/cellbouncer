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
#include <optimML/functions.h>
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
    fprintf(stderr, "Given one or more .assignments files from the same library, estimates\n");
    fprintf(stderr, "the likeliest underlying doublet rate (via maximum likelihood).\n");
    fprintf(stderr, "Also infers the frequency of each individual type in the data set.\n");
    fprintf(stderr, "If the same names are encountered in multiple assignments files, it is\n");
    fprintf(stderr, "assumed they correspond to the same individual/ID. If this is not the case\n");
    fprintf(stderr, "(for example, if you have clustered mitochondria from multiple species in the\n");
    fprintf(stderr, "same library, and they both have numeric IDs), be sure to change the IDs in one\n");
    fprintf(stderr, "of the files (in this example, create .ids files for both).\n");
    fprintf(stderr, "\nEXAMPLE 1:\n");
    fprintf(stderr, "If you pooled data from three species, separated these data using\n");
    fprintf(stderr, "demux_species, inferred mitochondrial haplotypes via demux_mt,\n");
    fprintf(stderr, "demultiplexed individuals from a VCF using demux_vcf, and processed\n");
    fprintf(stderr, "MULTIseq data on the same library using demux_multiseq, then you can\n");
    fprintf(stderr, "run doublet_dragon on all these output files to estimate the underlying\n");
    fprintf(stderr, "global doublet rate.\n");
    fprintf(stderr, "\nEXAMPLE 2:\n");
    fprintf(stderr, "If you identified individuals in the same data set multiple ways (i.e. via\n");
    fprintf(stderr, "cell-hashing, SNP-based demultiplexing, and mitochondrial clustering), then\n");
    fprintf(stderr, "make sure all of these .assignments files have the same names for the individuals\n");
    fprintf(stderr, "and run doublet_dragon to get a global doublet rate and find over/underrepresentation\n");
    fprintf(stderr, "of individual IDs in all of the .assignments files.\n");
    fprintf(stderr, "\nEXAMPLE 3:\n");
    fprintf(stderr, "You have just run demux_vcf on a single data set and want to find the most likely\n");
    fprintf(stderr, "doublet rate, including self-self doublets.\n");
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
    int grp_start = data_i.at("grp_start");
    int grp_end = data_i.at("grp_end");
    
    double grptot = 0;
    for (int i = grp_start; i <= grp_end; ++i){
        grptot += expit(params[i]);
    }
    double p1 = expit(params[p_idx]) / grptot;

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d)*p1 + d*p1 * p1;
    }
    else{
        // Doublet
        double p2 = expit(params[p_idx2]) / grptot;
        p = d*p1*p2;
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
    int grp_start = data_i.at("grp_start");
    int grp_end = data_i.at("grp_end");
    
    double grptot = 0;
    for (int i = grp_start; i <= grp_end; ++i){
        grptot += expit(params[i]);
    }
    double p1 = expit(params[p_idx]) / grptot;

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d)*p1 + d*p1 * p1;
    }
    else{
        // Doublet
        double p2 = expit(params[p_idx2])/grptot;
        p = d*p1*p2;
    }

    double dy_dp = (k - n*p)/(p - p*p);
    
    double e_negx1 = exp(-params[p_idx]);
    double e_negx1_p1_2 = pow(e_negx1 + 1, 2);
    double e_negx1_p1_3 = e_negx1_p1_2 * (e_negx1 + 1);
    double der_comp_1 = e_negx1 / (e_negx1_p1_2 * grptot) - e_negx1 / (e_negx1_p1_3 *grptot *grptot);

    if (p_idx2 == -1){
        double a = expit(params[p_idx])/grptot;
        double a2 = a*a;
        results[0] += dy_dp * (-a + a2);

        results[p_idx] += dy_dp * (1.0 - d + 2*d*a) * der_comp_1;
    }
    else{
        double a = expit(params[p_idx])/grptot;
        double b = expit(params[p_idx2])/grptot;
        results[0] += dy_dp * (a*b);
        results[p_idx] += dy_dp * (d*b) * der_comp_1;

        double e_negx2 = exp(-params[p_idx2]);
        double e_negx2_p1_2 = pow(e_negx2 + 1, 2);
        double e_negx2_p1_3 = e_negx2_p1_2 * (e_negx2 + 1);
        double der_comp_2 = e_negx2 / (e_negx2_p1_2 * grptot) - e_negx2/ ( e_negx2_p1_3 * grptot * grptot);

        results[p_idx2] += dy_dp * (d*a) * der_comp_2;
    }
}

int main(int argc, char *argv[]) {    
    
    if (argc < 2){
        help(1);
    } 
    
    vector<int> n;
    vector<int> k;
    vector<int> p_idx;
    vector<int> p_idx2;
    vector<int> grp_start;
    vector<int> grp_end;

    vector<double> weights;
    
    vector<double> params;
    params.push_back(0.01); // doublet rate
    
    vector<int> constr_01 = {0};
    
    map<int, string> idx2name;
    map<string, int> name2idx;
    
    map<int, int> idx2group;
    vector<set<string> > groups;

    map<int, int> params_div;

    map<string, map<string, int> > file2singlets;
    map<string, map<string, int> > file2doublets;
    map<string, int> file_tot;

    for (int i = 1; i < argc; ++i){
        map<string, int> counts_singlet;
        map<string, int> counts_doublet;
        map<string, double> weights_this;
        int totcells = parse_assignments(argv[i], counts_singlet, counts_doublet, weights_this);
        if (counts_doublet.size() == 0){
            fprintf(stderr, "ERROR: no doublet identifications in %s\n", argv[i]);
            exit(1);
        }
        
        file2singlets.insert(make_pair(argv[i], counts_singlet));
        file2doublets.insert(make_pair(argv[i], counts_doublet));
        file_tot.insert(make_pair(argv[i], totcells));
            
        int tot_singlet = 0;
        for (map<string, int>::iterator cs = counts_singlet.begin();
            cs != counts_singlet.end(); ++cs){
            tot_singlet += cs->second;
        }
        
        int grp_idx_this = -1;

        for (map<string, int>::iterator cs = counts_singlet.begin(); 
            cs != counts_singlet.end(); ++cs){
            
            int p_idx_this;
            if (name2idx.count(cs->first) > 0){
                // If name matches an already-existing name, treat them as the
                // same variable.
                p_idx_this = name2idx[cs->first];

                // Make the initial value the mean initial value across all data sets
                if (params_div.count(p_idx_this) == 0){
                    params_div.insert(make_pair(p_idx_this, 1));
                }
                params_div[p_idx_this]++;
                params[p_idx_this] += ((double)cs->second / (double)tot_singlet);
                
                grp_idx_this = idx2group[p_idx_this];
            }
            else{
                p_idx_this = params.size();
                params.push_back((double)cs->second / (double)tot_singlet);
                constr_01.push_back(p_idx_this);
                name2idx.insert(make_pair(cs->first, p_idx_this));
                idx2name.insert(make_pair(p_idx_this, cs->first));
                
                if (grp_idx_this == -1){
                    grp_idx_this = groups.size();
                    set<string> s;
                    groups.push_back(s);
                    groups[groups.size()-1].insert(cs->first);
                }
                else{
                    groups[grp_idx_this].insert(cs->first);
                }
            }
                    
            n.push_back(totcells);
            k.push_back(cs->second);
            p_idx.push_back(p_idx_this);
            p_idx2.push_back(-1);
            weights.push_back(weights_this[cs->first]);
        }
        for (map<string, int>::iterator cd = counts_doublet.begin();
            cd != counts_doublet.end(); ++cd){
            size_t pos = cd->first.find("+");
            if (pos == string::npos){
                fprintf(stderr, "ERROR: unable to identify components of doublet ID %s\n", cd->first.c_str());
                exit(1);
            }
            string id1 = cd->first.substr(0, pos);
            string id2 = cd->first.substr(pos+1, cd->first.length()-pos-1);
            if (name2idx.count(id1) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s (%s)\n", id1.c_str(), argv[i]);
                exit(1);
            }
            if (name2idx.count(id2) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s (%s)\n", id2.c_str(), argv[i]);
                exit(1);
            }
            n.push_back(totcells);
            k.push_back(cd->second);
            p_idx.push_back(name2idx[id1]);
            p_idx2.push_back(name2idx[id2]);
            weights.push_back(weights_this[cd->first]);
        }
    }
    // Use means of initial guesses for each singlet proportion that was 
    // present in multiple data sets.
    for (map<int, int>::iterator pd = params_div.begin(); pd != params_div.end(); ++pd){
        params[pd->first] /= (double)pd->second;
    }
    
    for (int i = 0; i < p_idx.size(); ++i){
        int gi = idx2group[p_idx[i]];
        int gimax = -1;
        int gimin = -1;
        for (set<string>::iterator g = groups[gi].begin(); g != groups[gi].end(); ++g){
            int ni = name2idx[*g];
            if (gimax == -1 || ni > gimax){
                gimax = ni;
            }
            if (gimin == -1 || ni < gimin){
                gimin = ni;
            }
        }
        grp_start.push_back(gimin);
        grp_end.push_back(gimax);
    }

    // Make sure groups sum to 1 (and transform them)
    for (vector<set<string> >::iterator grp = groups.begin(); grp != groups.end(); ++grp){
        double tot = 0.0;
        for (set<string>::iterator member = grp->begin(); member != grp->end(); ++member){
            int idx = name2idx[*member];
            tot += params[idx];
        }
        for (set<string>::iterator member = grp->begin(); member != grp->end(); ++member){
            int idx = name2idx[*member];
            params[idx] /= tot;
            params[idx] = logit(params[idx]);
        }
    }
    
    optimML::multivar_ml_solver solver(params, loglik, dloglik_dx);
    
    // Constrain doublet rate; handle other constraints in the LL/gradient functions
    solver.constrain_01(0);
    
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_idx", p_idx);
    solver.add_data("p_idx2", p_idx2);
    solver.add_data("grp_start", grp_start);
    solver.add_data("grp_end", grp_end);

    solver.add_weights(weights);
    solver.solve();
    
    double doublet_rate = solver.results[0];
    map<string, double> indv_freq;
    for (map<string, int>::iterator ni = name2idx.begin(); ni != name2idx.end(); ++ni){
        indv_freq.insert(make_pair(ni->first, expit(solver.results[ni->second])));
        fprintf(stderr, "%s\t%f\n", ni->first.c_str(), indv_freq[ni->first]);
    }
    fprintf(stderr, "Global doublet rate: %f\n", doublet_rate);
    
    for (map<string, int>::iterator ft = file_tot.begin(); ft != file_tot.end(); ++ft){
        for (map<string, int>::iterator d = file2doublets[ft->first].begin(); d != 
            file2doublets[ft->first].end(); ++d){
            size_t pos = d->first.find("+");
            string id1 = d->first.substr(0, pos);
            string id2 = d->first.substr(pos+1, d->first.length()-pos-1);
            double p1 = indv_freq[id1];
            double p2 = indv_freq[id2];
            double expected = doublet_rate*p1*p2*(double)ft->second;
            fprintf(stderr, "%s\t%s\t%d\t%f\n", ft->first.c_str(), d->first.c_str(),
                d->second, expected);
        }
    }

    return 0;
}
