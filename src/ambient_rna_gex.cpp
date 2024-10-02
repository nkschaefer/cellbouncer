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
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <cfloat>
#include <random>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "ambient_rna_gex.h"

using std::cout;
using std::endl;
using namespace std;

contam_profiler_gex::contam_profiler_gex(robin_hood::unordered_map<unsigned long, double>& contam_rate,
    map<int, double>& contam_prof, 
    robin_hood::unordered_map<unsigned long, int>& assn,
    int n_samples, 
    bool doublets_as_indvs){
    
    this->contam_rate = &contam_rate;
    this->contam_prof = &contam_prof;
    this->assn = &assn;
    this->n_samples = n_samples;

    this->mtx = NULL;
    this->clusts = NULL;
    this->ase_mtx1 = NULL;
    this->ase_mtx2 = NULL;
    this->nthreads = 0;
    this->doublets_as_indvs = doublets_as_indvs;
    this->profile_set = false;
    this->has_ase = false;
    this->round_counts_mtx = false;
}

void contam_profiler_gex::set_clusts(robin_hood::unordered_map<unsigned long, int>& clusts,
    int n_clusts){
    
    this->clusts = &clusts;
    this->n_clusts = n_clusts;    

}

void contam_profiler_gex::set_mtx(robin_hood::unordered_map<unsigned long, map<int, long int> >& mtx,
    int n_features){
    
    this->mtx = &mtx;
    this->n_features = n_features;
}

void contam_profiler_gex::set_threads(int nt){
    if (nt < 2){
        nt = 0;
    }
    this->nthreads = nt;
}

void contam_profiler_gex::set_ase(robin_hood::unordered_map<unsigned long, map<int, long int> >& ase1,
    robin_hood::unordered_map<unsigned long, map<int, long int > >& ase2){
    
    this->ase_mtx1 = &ase1;
    this->ase_mtx2 = &ase2;
    has_ase = true;
}

void contam_profiler_gex::round_counts(){
    this->round_counts_mtx = true;
    //srand(time(NULL));
    rand_gen = mt19937(rand_dev());
    rand_dist = uniform_real_distribution<>(0.0, 1.0);
}

/**
 * optimML-format log likelihood function for contamination profile inference
 */
double ll_multinom(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int num = data_i.at("num");
    double c = data_d.at("c");
    int gi = data_i.at("grp_idx");

    vector<double> p;
    vector<double> n;
    
    char buf[30];
    string bufstr;
    
    double xsum = 1;
    double term2 = 0;
    double term3 = 0;
    double psum = 0.0;
    
    int signp;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        int idx = data_i.at(bufstr);
        if (idx == -1){
            break;
        }
        double p_i = c*params[idx] + (1.0-c)*params[gi*num+idx];
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        double n_i = data_d.at(bufstr);
        xsum += n_i;
        term2 += lgammaf_r(n_i + 1, &signp);
        term3 += n_i*log(p_i);
    }
    double term1 = lgammaf_r(xsum, &signp);
    return (term1 - term2 + term3);
}

/**
 * optimML-format derivative log likelihood function for contamination profile
 * inference
 */
void dll_multinom(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
     
    int num = data_i.at("num");
    double c = data_d.at("c");
    double tot = data_d.at("tot");
    int gi = data_i.at("grp_idx");

    char buf[30];
    string bufstr;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        int idx = data_i.at(bufstr);
        if (idx == -1){
            break;
            //results[idx] -= tot;
            //results[num+idx] -= tot;
        }
        else{
            double p_i = c*params[idx] + (1.0-c)*params[gi*num+idx];
            sprintf(&buf[0], "n_%d", i);
            bufstr = buf;
            double n_i = data_d.at(bufstr);
            double dLL_dpi = (n_i/p_i);
            results[idx] += (dLL_dpi)*c;
            results[gi*num+idx] += (dLL_dpi)*(1.0-c);
        }
    }
}

/**
 * Infers the parameters of a multinomial distribution representing
 * ambient RNA contamination, as well as the average expression of
 * each gene in each cluster.
 */
bool contam_profiler_gex::get_profile(){
    if (clusts == NULL || mtx == NULL){
        fprintf(stderr, "ERROR: clusters and matrix are required.\n");
        return false;
    }

    vector<double> grp1;
    vector<vector<double> > grps2;
    vector<vector<double> > ns;
    vector<vector<int> > ns_idx;
    vector<double> tots;
    vector<double> cs;
    vector<int> grp_idx;
    double sum_totxc = 0.0;
    double sum_totx1mc = 0.0;
    double ncsum = 0.0;
    
    vector<double> grps2tot;
    int n_grp2 = n_clusts;
    if (n_clusts == 0){
        n_grp2 = 1;
    }
    for (int i = 0; i < n_grp2; ++i){
        vector<double> v;
        grps2.push_back(v);
        grps2tot.push_back((double)n_features);
    }

    for (int i = 0; i < n_features; ++i){
        grp1.push_back(1.0);
        for (int j = 0; j < n_grp2; ++j){
            grps2[j].push_back(1.0);
        }
        vector<double> v;
        ns.push_back(v);
        vector<int> v2;
        ns_idx.push_back(v2);
    }
    double grp1tot = (double)n_features;
    int num_added = 0;
    double grp2tot = (double)n_features;
    
    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator m = mtx->begin();
        m != mtx->end(); ++m){

        int clust = -1;
        if (clusts->count(m->first) > 0){
            clust = (*clusts)[m->first];
        }
        int cell_assn = -1;
        if (assn->count(m->first) > 0){
            cell_assn = (*assn)[m->first];
        }

        double celltot = 0;
        
        int i = 0;
        int cprev = -1;
        for (map<int, long int>::iterator c = m->second.begin(); c != m->second.end(); ++c){
            double count = (double)c->second;
            
            celltot += count;
            if (cell_assn >= 0){
                if (cell_assn >= n_samples){
                    pair<int, int> comb = idx_to_hap_comb(cell_assn, n_samples);
                    grp1[comb.first] += 0.5*(*contam_prof)[comb.first]*count;
                    grp1[comb.second] += 0.5*(*contam_prof)[comb.second]*count;
                    grp1tot += 0.5*(*contam_prof)[comb.first]*count + 0.5*(*contam_prof)[comb.second]*count;
                }
                else{ 
                    grp1[c->first] += (*contam_prof)[cell_assn] * count;
                    grp1tot += (*contam_prof)[cell_assn] * count;
                }
            }
            if (clust >= 0){
                grps2[clust][c->first] += count;
                grps2tot[clust] += count;

                ns_idx[i].push_back(c->first);
                for (int j = cprev+1; j < c->first; ++j){
                    ns[j].push_back(0);
                }
                ns[c->first].push_back(count);
                i++;
                cprev = c->first;
            }
        }
        
        cells_tot.emplace(m->first, celltot);

        if (clust >= 0){
            for (int j = i; j < n_features; ++j){
                ns_idx[j].push_back(-1);
            }
            for (int j = cprev + 1; j < n_features; ++j){
                ns[j].push_back(0.0);
            }

            tots.push_back(celltot);
            cs.push_back((*contam_rate)[m->first]);
            grp_idx.push_back(clust + 1);
            ++num_added;
        }
    }
    
    for (int j = 0; j < grps2.size(); ++j){ 
        for (int i = 0; i < n_features; ++i){
            grps2[j][i] /= grp2tot;
            if (j == 0){
                grp1[i] /= grp1tot;
            }
        }
    }
    vector<double> params = {};
    optimML::multivar_ml_solver mnsolver(params, ll_multinom, dll_multinom);
    if (nthreads > 1){
        mnsolver.set_threads(nthreads);
        mnsolver.set_bfgs_threads(nthreads);
    }
    mnsolver.add_param_grp(grp1);
    for (int i = 0; i < grps2.size(); ++i){
        mnsolver.add_param_grp(grps2[i]);
    }
    mnsolver.add_data_fixed("num", n_features);
    mnsolver.add_data("tot", tots);
    mnsolver.add_data("c", cs);
    mnsolver.add_data("grp_idx", grp_idx);
    char buf[30];
    for (int i = 0; i < n_features; ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        mnsolver.add_data(bufstr, ns[i]);
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        mnsolver.add_data(bufstr, ns_idx[i]);
    }
    //mnsolver.set_delta(1);
    fprintf(stderr, "Inferring ambient RNA expression profile...\n");
    
    prof_ambient.clear();
    prof_clusts.clear();
    profile_set = false;

    try{
        mnsolver.solve();
        
        // Extract profiles and save them.
        for (int i = 0; i < n_grp2; ++i){
            vector<double> v;
            prof_clusts.push_back(v);
        }
        for (int i = 0; i < n_features; ++i){
            prof_ambient.push_back(mnsolver.results[i]);
            for (int j = 0; j < n_grp2; ++j){
                prof_clusts[j].push_back(mnsolver.results[n_features*(j+1) + i]);
            }
        }
        profile_set = true;
    } 
    catch(int errcode){
        // Something went wrong with solution.
        return false;
    }
    
    if (has_ase){
        // Infer ASE-ness of each gene
        infer_ase();
    }

    return true;
}

bool contam_profiler_gex::infer_ase(){
    if (!profile_set){
        fprintf(stderr, "ERROR: contam profile not yet inferred\n");
        return false;
    }
    if (!has_ase){
        fprintf(stderr, "ERROR: no ASE data provided\n");
        return false;
    }
    
    fprintf(stderr, "Inferring allelic bias of contamination...\n");
    ase_fracs.clear();
    for (int i = 0; i < n_features; ++i){
        // Initialize to 100% non-allelic (no SNPs)
        vector<double> row{ 0.0, 0.0, 1.0 };
        
        ase_fracs.push_back(row);
    }
    return true;    
}

/**
 * Modifies the gene expression matrix to remove contamination from it.
 */
bool contam_profiler_gex::decontam(){
    if (!profile_set){
        return false;
    }
    mtx_decontam.clear();
    fprintf(stderr, "Creating adjusted gene expression matrix...\n");
    
    char buf[30];
    string bufstr;

    for (robin_hood::unordered_map<unsigned long, map<int, long int> >::iterator m = mtx->begin();
        m != mtx->end(); ++m){
        map<int, double> mtmp;
        mtx_decontam.emplace(m->first, mtmp);
        
        double celltot = cells_tot[m->first];
        double c = (*contam_rate)[m->first];
        
        // How many total counts need to be removed from this cell?
        double to_rm = celltot*c;
        
        // How much of the ambient RNA gene expression vector is accounted for
        // with nonzero entries in the expression matrix?
        double ambtot = 0.0;
        
        for (map<int, long int>::iterator m2 = m->second.begin(); m2 != m->second.end(); ++m2){
            ambtot += prof_ambient[m2->first];
        }

        deque<pair<double, int> > rms;
        for (map<int, long int>::iterator m2 = m->second.begin(); m2 != m->second.end(); ++m2){
            double rm_this = to_rm * (prof_ambient[m2->first] / ambtot);
            double after_rm = (double)m2->second - rm_this;
            rms.push_back(make_pair(after_rm, m2->first)); 
        }
        sort(rms.begin(), rms.end());
        
        double removed = 0.0;

        bool elim_neg = rms[0].first <= 0;
        while (elim_neg){
            while (rms[0].first <= 0){
                // Anything negative is eliminated.
                //fprintf(stderr, "%d) %f -> 0\n", rms[0].second, (double)m->second[rms[0].second]);
                to_rm -= (double)m->second[rms[0].second];
                removed += (double)m->second[rms[0].second];
                rms.pop_front();
            }
            for (int i = 0; i < rms.size(); ++i){
                double rm_this = to_rm * prof_ambient[rms[i].second];
                rms[i].first = (double)m->second[rms[i].second] - rm_this;
            }
            sort(rms.begin(), rms.end());
            elim_neg = rms[0].first <= 0;
        }
        while (rms.size() > 0){
            
            if (round_counts_mtx){
                double remainder = rms[0].first - floor(rms[0].first);
                
                //double r = ((double) rand() / (RAND_MAX));
                double r = rand_dist(rand_gen);
                if (r < remainder){
                    rms[0].first = ceil(rms[0].first);
                }
                else{
                    rms[0].first = floor(rms[0].first);
                }
            }
            
            removed += ((double)m->second[rms[0].second] - rms[0].first);
            mtx_decontam[m->first].insert(make_pair(rms[0].second, rms[0].first)); 
            rms.pop_front();
        }
    }
    
    return true;
}


