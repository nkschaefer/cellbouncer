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
#include <mixtureDist/functions.h>
#include "common.h"
#include "demux_vcf_llr.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * ===== Contains functions relating to identifying individuals of  ===== 
 * ===== origin, used by demux_vcf.                                 =====
 */


llr_table::llr_table(int x){
    n_indvs = 0;
    int n_elt = x + (int)round(pow(2, binom_coef_log(x, 2)));
    included.reserve(n_elt);
    maxllr.reserve(n_elt);
};
llr_table::~llr_table(){
    lookup_llr.clear();
    included.clear();
    maxllr.clear();
    minllr.clear();
}

void llr_table::print(string& bc_str, vector<string>& samples){
    
    for (map<double, vector<pair<short, short> > >::iterator x = lookup_llr.begin();
        x != lookup_llr.end(); ++x){
        for (int i = 0; i < x->second.size(); ++i){
            if (included[x->second[i].first] && included[x->second[i].second]){
                string n1 = idx2name(x->second[i].first, samples);
                string n2 = idx2name(x->second[i].second, samples);
                fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                    n1.c_str(), n2.c_str(), -x->first);
                fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                    n2.c_str(), n1.c_str(), x->first);
            }
        }
    }
}

void llr_table::print_ranges(string& barcode, vector<string>& samples){
    int n_samples = samples.size();
    for (int i = 0; i < n_samples; ++i){
        if (included[i]){
            fprintf(stdout, "%s\t%s\t%f\t%f\n", barcode.c_str(), idx2name(i, samples).c_str(), 
                minllr[i], maxllr[i]);
        }
        for (int j = i + 1; j < n_samples; ++j){
            int k = hap_comb_to_idx(i, j, n_samples);
            fprintf(stdout, "%s\t%s\t%f\t%f\n", barcode.c_str(), idx2name(k, samples).c_str(),
                minllr[k], maxllr[k]);
        }
    }
}

void llr_table::insert(short i1, short i2, double llr){
      
    while(included.size() < i1+1){
        included.push_back(false);
        maxllr.push_back(0.0);
        minllr.push_back(0.0);
    }
    while(included.size() < i2+1){
        included.push_back(false);
        maxllr.push_back(0.0);
        minllr.push_back(0.0);
    }
    if (!included[i1]){
        ++n_indvs;
    }
    if (!included[i2]){
        ++n_indvs;
    }
    if (!included[i1] || maxllr[i1] < llr){
        maxllr[i1] = llr;
    }
    if (!included[i1] || minllr[i1] > llr){
        minllr[i1] = llr;
    }
    if (!included[i2] || maxllr[i2] < -llr){
        maxllr[i2] = -llr;
    }
    if (!included[i2] || minllr[i2] > -llr){
        minllr[i2] = -llr;
    }
    included[i1] = true;
    included[i2] = true;
    
    short low;
    short high;
    if (llr > 0){
        low = i2;
        high = i1;
        llr = -llr;
    }
    else{
        low = i1;
        high = i2;
    }
    if (lookup_llr.count(llr) == 0){
        vector<pair<short, short> > v;
        lookup_llr.emplace(llr, v);
    }
    lookup_llr[llr].push_back(make_pair(low, high));
};
        
void llr_table::disallow(short i){
    if (i < included.size()){
        if (included[i]){
            maxllr[i] = 0.0;
            minllr[i] = 0.0;
            --n_indvs;
        }
        included[i] = false;
    }
}

bool llr_table::del(int n_keep){

    if (n_indvs < n_keep){
        return false;
    }
    
    it = lookup_llr.begin();
    while (n_indvs > n_keep && it != lookup_llr.end()){
        bool del_all = it->second.size() <= (n_indvs - n_keep);
        map<double, int> maxllr_counts;
        int indv_tot = 0;
        for (vector<pair<short, short> >::iterator x = it->second.begin();
            x != it->second.end();){
            if (!included[x->first] || !included[x->second]){
                it->second.erase(x);
            }
            else if (del_all){
                if (included[x->first]){
                    --n_indvs;
                    included[x->first] = false;
                    minllr[x->first] = 0.0;
                    maxllr[x->first] = 0.0;
                }
                it->second.erase(x);
            }
            else{
                double mllr = maxllr[x->first];
                if (maxllr_counts.count(mllr) == 0){
                    maxllr_counts.insert(make_pair(mllr, 1));
                }
                else{
                    maxllr_counts[mllr]++;
                }
                indv_tot++;
                ++x;
            }
        }
        
        if (indv_tot <= n_indvs - n_keep){
            del_all = true;
        }

        if (!del_all){
            double cutoff = 0.0;
            int runtot = 0;
            for (map<double, int>::iterator m = maxllr_counts.begin(); m != maxllr_counts.end();
                ++m){
                runtot += m->second;
                cutoff = m->first;
                if (runtot >= n_indvs - n_keep){
                    break;
                }
            }
            for (vector<pair<short, short> >::iterator x = it->second.begin(); x != it->second.end(); ){
                if (maxllr[x->first] <= cutoff){
                    if (included[x->first]){
                        included[x->first] = false;
                        maxllr[x->first] = 0.0;
                        minllr[x->first] = 0.0;
                        n_indvs--;
                    }
                    it->second.erase(x);
                }
                else{
                    ++x;
                }
            }
        }
        else{
            for (vector<pair<short, short> >::iterator x = it->second.begin(); x != 
                it->second.end(); ){
                if (included[x->first]){
                    included[x->first] = false;
                    n_indvs--;
                    maxllr[x->first] = 0.0;
                    minllr[x->first] = 0.0;
                }
                it->second.erase(x);
            }
            lookup_llr.erase(it++);
        }
    }

    if (n_indvs != n_keep){
        return false;
    }
    return true;
}
        
void llr_table::get_max(int& best_idx, double& best_llr){
    if (n_indvs < 2){
        best_idx = -1;
        best_llr = 0.0;
        return;
    }
    best_idx = -1;
    best_llr = 0.0;
    bool success = del(2);
    if (!success){
        return;
    }
    if (n_indvs == 2){
        while (best_idx == -1 && it != lookup_llr.end()){
            for (vector<pair<short, short> >::iterator x = it->second.begin(); 
                x != it->second.end(); ){
                if (included[x->first] && included[x->second]){
                    best_idx = x->second;
                    best_llr = -it->first;
                    break;
                }
                it->second.erase(x);
            } 
            if (best_idx == -1){
                ++it;
            }
        }
    }
    return;
}

/** 
 * Helper function used by compute_kcomps().
 * 
 * Computes comparisons between two doublet combinations conditional on
 * a shared individual. 
 *
 * In other words, given that individual A is involved, can compute the
 * LLR of (A+B) vs (A+C), by considering the LLR of (A+B)/A vs (A+C)/A.
 */
void get_kcomps_cond(map<int, map<int, double> >& kcomps,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments){

    for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end();
        ++llr){

        vector<int> js;
        vector<int> comps;
        for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != llr->second.end();
            ++llr2){
            if (llr2->first >= n_samples){
                //if (allowed_assignments.size() == 0 || allowed_assignments.find(llr2->first) != 
                //    allowed_assignments.end()){
                
                if (!tab.included[llr2->first]){
                //    continue;
                }

                if (true){
                    pair<int, int> comb = idx_to_hap_comb(llr2->first, n_samples);
                    int j = -1;
                    if (comb.first == llr->first){
                        j = comb.second;
                    }
                    else{
                        j = comb.first;
                    }
                    js.push_back(j);
                    comps.push_back(-llr2->second);
                }
            }
        }
        
        if (comps.size() > 1){
            for (int i = 0; i < comps.size()-1; ++i){
                int idx1 = js[i];
                for (int j = i +1; j < comps.size(); ++j){
                    int idx2 = js[j];
                    double llr = comps[i] - comps[j];
                    if (idx1 < idx2){
                        if (kcomps.count(idx1) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx1, m));
                        }
                        if (kcomps[idx1].count(idx2) == 0){
                            kcomps[idx1].insert(make_pair(idx2, 0.0));
                        }
                        kcomps[idx1][idx2] += llr;
                    }
                    else{
                        if (kcomps.count(idx2) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx2, m));
                        }
                        if (kcomps[idx2].count(idx1) == 0){
                            kcomps[idx2].insert(make_pair(idx1, 0.0));
                        }
                        kcomps[idx2][idx1] += llr;
                    }
                }
            }
        }
    }
}

/**
 * Used by populate_llr_table(). 
 *
 * Whereas that function can directly
 * compare the likelihood of a doublet combination with one of its
 * singlet components (i.e. individual A+B vs individual A and 
 * individual A+B vs individual B), this function computes all other
 * possible comparisons involving doublet combinations (i.e.
 * individual A+B vs individual C, individual A+B vs individual A+D,
 * and individual A+B vs individual C+D.
 * 
 * Puts these comparisons into the same log likelihood ratio table
 * so a single decision can be made across all comparisons.
 */
void compute_k_comps(map<int, map<int, double> >& llrs,
    llr_table& tab,
    vector<int>& ks,
    int n_samples,
    set<int>& allowed_assignments,
    double doublet_rate){
    
    // Get a table of LLRs between different combination members, conditional on a cell
    // being half another member. In other words, keys/values here are
    // ID1 -> ID2A -> ID2B -> LLR
    // LLR is the log likelihood ratio of ID2A vs ID2B being included in a double model
    // with ID1. 
    map<int, map<int, map<int, double> > > kcomp_cond;
    for (map<int, map<int, double> >::iterator samp = llrs.begin(); samp != llrs.end(); ++samp){
        if (!tab.included[samp->first]){
            continue;
        }
        map<int, map<int, double> > llrstmp;
        llrstmp.insert(make_pair(samp->first, samp->second));
        map<int, map<int, double> > kcomp; 
        get_kcomps_cond(kcomp, llrstmp, tab, n_samples, allowed_assignments);
        kcomp_cond.insert(make_pair(samp->first, kcomp));
    }
    // compute all single vs double comparisons
    for (int i = 0; i < n_samples; ++i){
        if (allowed_assignments.size() != 0 && 
            allowed_assignments.find(i) == allowed_assignments.end()){
            continue;
        }
        else if (!tab.included[i]){
            continue;
        }

        for (int ki = 0; ki < ks.size(); ++ki){
            // NOTE: ks is already a filtered list (allowable)
            int k = ks[ki];
            pair<int, int> comb = idx_to_hap_comb(k, n_samples);
            if (comb.first != i && comb.second != i){
                double ll1 = llrs[comb.first][k];
                double ll2 = llrs[comb.second][k];
                if (comb.first < i){
                    ll1 += -llrs[comb.first][i];
                }
                else{
                    ll1 += llrs[i][comb.first];
                }
                if (comb.second < i){
                    ll2 += -llrs[comb.second][i];
                }
                else{
                    ll2 += llrs[i][comb.second];
                }
                double llr = 0.5*ll1 + 0.5*ll2;
                if (doublet_rate != 0.5 && doublet_rate < 1.0){
                    // i is singlet, k is doublet
                    llr += log2(1.0 - doublet_rate) - log2(doublet_rate);
                }
                tab.insert(i, k, llr);
            }
            else{
                // It's already been computed.
                double llr = llrs[i][k];
                if (doublet_rate != 0.5 && doublet_rate < 1.0){
                    llr += log2(1.0 - doublet_rate) - log2(doublet_rate);
                }
                tab.insert(i, k, llr);
            }
        }
    }

    // compute all double vs double comparisons
    for (int ki = 0; ki < ks.size()-1; ++ki){
        int k1 = ks[ki];
        
        pair<int, int> comb1 = idx_to_hap_comb(k1, n_samples);
        for (int kj = ki+1; kj < ks.size(); ++kj){
            int k2 = ks[kj];
            pair<int, int> comb2 = idx_to_hap_comb(k2, n_samples);
            vector<double> llr_parts;
            int a = comb1.first;
            int b = comb1.second;
            int c = comb2.first;
            int d = comb2.second;

            if (a == c){
                if (b < d){
                    llr_parts.push_back(kcomp_cond[a][b][d]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][d][b]);
                }
            }
            else if (a == d){
                if (b < c){
                    llr_parts.push_back(kcomp_cond[a][b][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][c][b]);
                }
            }
            else if (b == c){
                if (a < d){
                    llr_parts.push_back(kcomp_cond[b][a][d]);
                }         
                else{
                    llr_parts.push_back(-kcomp_cond[b][d][a]);
                }
            }
            else if (b == d){
                if (a < c){
                    llr_parts.push_back(kcomp_cond[b][a][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[b][c][a]);
                }
            }
            else{
                double llr1 = 0;
                double llr2 = 0;
                double llr3 = 0;
                double llr4 = 0;
                double llr5 = 0;
                double llr6 = 0;
                double llr7 = 0;
                double llr8 = 0;
                // Calculate (A+B)/(C+D) 
                if (b < c){
                    llr1 += kcomp_cond[a][b][c];
                }
                else{
                    llr1 += -kcomp_cond[a][c][b];
                }
                if (a < d){
                    llr1 += kcomp_cond[c][a][d];
                }
                else{
                    llr1 += -kcomp_cond[c][d][a];
                }
                
                if (a < c){
                    llr2 += kcomp_cond[b][a][c];
                }
                else{
                    llr2 += -kcomp_cond[b][c][a];
                }
                if (b < d){
                    llr2 += kcomp_cond[c][b][d];
                }
                else{
                    llr2 += -kcomp_cond[c][d][b];
                }

                if (b < d){
                    llr3 += kcomp_cond[a][b][d];
                }
                else{
                    llr3 += -kcomp_cond[a][d][b];
                }
                if (a < c){
                    llr3 += kcomp_cond[d][a][c];
                }
                else{
                    llr3 += -kcomp_cond[d][c][a];
                }

                if (a < d){
                    llr4 += kcomp_cond[b][a][d];
                }
                else{
                    llr4 += -kcomp_cond[b][d][a];
                }
                if (b < c){
                    llr4 += kcomp_cond[d][b][c];
                }
                else{
                    llr4 += -kcomp_cond[d][c][b];
                }

                if (d < a){
                    llr5 += kcomp_cond[c][d][a];
                }
                else{
                    llr5 += -kcomp_cond[c][a][d];
                }
                if (c < b){
                    llr5 += kcomp_cond[a][c][b];
                }
                else{
                    llr5 += -kcomp_cond[a][b][c];
                }
                llr5 = -llr5;

                if (c < a){
                    llr6 += kcomp_cond[d][c][a];
                }
                else{
                    llr6 += -kcomp_cond[d][a][c];
                }
                if (d < b){
                    llr6 += kcomp_cond[a][d][b];
                }
                else{
                    llr6 += -kcomp_cond[a][b][d];
                }
                llr6 = -llr6;

                if (d < b){
                    llr7 += kcomp_cond[c][d][b];
                }
                else{
                    llr7 += -kcomp_cond[c][b][d];
                }
                if (c < a){
                    llr7 += kcomp_cond[b][c][a];
                }
                else{
                    llr7 += -kcomp_cond[b][a][c];
                }
                llr7 = -llr7;

                if (c < b){
                    llr8 += kcomp_cond[d][c][b];
                }
                else{
                    llr8 += -kcomp_cond[d][b][c];
                }
                if (d < a){
                    llr8 += kcomp_cond[b][d][a];
                }
                else{
                    llr8 += -kcomp_cond[b][a][d];
                }
                llr8 = -llr8;
                llr_parts.push_back(llr1);
                llr_parts.push_back(llr2);
                llr_parts.push_back(llr3);
                llr_parts.push_back(llr4);
                llr_parts.push_back(llr5);
                llr_parts.push_back(llr6);
                llr_parts.push_back(llr7);
                llr_parts.push_back(llr8);
            }
            
            float llrsum = 0.0;
            float llrcount = 0;
            for (int idx = 0; idx < llr_parts.size(); ++idx){
                llrsum += llr_parts[idx];
                llrcount++;
            }
            double llr = llrsum / llrcount;
            tab.insert(k1, k2, llr);
        }
    }
}

/**
 * Given a set of allele counts (at all possible SNP types) for a single cell,
 * populates a log likelihood ratio table for that cell, which gives the LLR
 * of every possible identity vs every other possible identity.
 *
 * allowed_assignments is a filtered list of possible identities cells can take on.
 *  It must include all possible singlet identities (if the user has restricted the
 *  possible doublet identities, we must still have all singlet identities included
 *  in those doublet identities to build the LLR table).
 *
 * allowed_assignments2 is allowed_assignments, but without singlets that are only
 *  there because they make up allowed doublet identities (this is the actual 
 *  filtered list of possible identities).
 */
bool populate_llr_table(map<pair<int, int>, 
    map<pair<int, int>, pair<float, float> > >& counts,
    map<int, map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool incl_contam,
    double contam_rate,
    map<pair<int, int>, map<pair<int, int>, double> >* amb_fracs){
    
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
        counts.begin(); y != counts.end(); ++y){
        
        if (allowed_assignments.size() > 0 && allowed_assignments.find(y->first.first) == 
            allowed_assignments.end()){
            continue;
        }
        
        // Set default expectation for indv1
        // 0 = homozygous ref (~0% alt allele)
        // 1 = heterozygous (~50% alt allele)
        // 2 = homozygous alt (~100% alt allele)
        double exp1 = adjust_p_err((double)y->first.second / 2.0, error_rate_ref, error_rate_alt);
        /*
        float exp1 = error_rate_ref;
        if (y->first.second == 1){
            //exp1 = 0.5;
            exp1 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
        }
        else if (y->first.second == 2){
            exp1 = 1.0-error_rate_alt;
        }
        */

        for (map<pair<int, int>, pair<float, float> >::iterator z = 
            y->second.begin(); z != y->second.end(); ++z){
            
            if (allowed_assignments.size() > 0 && allowed_assignments.find(z->first.first) ==
                allowed_assignments.end()){
                continue;
            } 

            // If same site type, we can't distinguish between the
            // two individuals from this piece of information
            if (z->first.first != -1 && y->first.second != z->first.second){
                
                if (incl_contam){
                    exp1 = (1.0-contam_rate)*((double)y->first.second/2.0) + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp1 = adjust_p_err(exp1, error_rate_ref, error_rate_alt);
                }

                // Set default expectation for indv2
                double exp2 = adjust_p_err((double)z->first.second/2.0, error_rate_ref, error_rate_alt);
                /*
                float exp2 = error_rate_ref;
                if (z->first.second == 1){
                    //exp2 = 0.5;
                    exp2 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
                }
                else if (z->first.second == 2){
                    exp2 = 1.0-error_rate_alt;
                }
                */
                if (incl_contam){
                    exp2 = (1.0-contam_rate)*((double)z->first.second/2.0) + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp2 = adjust_p_err(exp2, error_rate_ref, error_rate_alt);
                }
                
                double exp3 = adjust_p_err((double)(y->first.second + z->first.second)/4.0, 
                    error_rate_ref, error_rate_alt);
                /*
                float exp3;
                if (y->first.second == 0 && z->first.second == 0){
                    exp3 = error_rate_ref;
                }
                else if (y->first.second == 2 && z->first.second == 2){
                    exp3 = 1.0 - error_rate_alt;
                }
                else if ((y->first.second == 1 && z->first.second == 1) ||
                        (y->first.second == 0 && z->first.second == 2) ||
                        (y->first.second == 2 && z->first.second == 0)){
                    exp3 = 0.5*( 1.0 - error_rate_alt + error_rate_ref);
                }
                else if ((y->first.second == 0 && z->first.second == 1) ||
                    (y->first.second == 1 && z->first.second == 0)){
                    exp3 = 0.25*(1 - error_rate_alt + 3*error_rate_ref);
                }
                else if ((y->first.second == 1 && z->first.second == 2) ||
                    (y->first.second == 2 && z->first.second == 1)){
                    exp3 = 0.25*(3.0 - 3*error_rate_alt + error_rate_ref);
                }
                */
                if (incl_contam){
                    exp3 = (1.0-contam_rate)*(double)(y->first.second + z->first.second)/4.0 + 
                        contam_rate*((*amb_fracs)[y->first][z->first]);
                    exp3 = adjust_p_err(exp3, error_rate_ref, error_rate_alt);
                }

                int i = y->first.first;
                int j = z->first.first;
                int k = hap_comb_to_idx(i, j, n_samples);

                int ref = (int)round(z->second.first);
                int alt = (int)round(z->second.second);
                
                double ll1 = dbinom(ref+alt, alt, exp1);   
                double ll2 = dbinom(ref+alt, alt, exp2);
                double ll3 = dbinom(ref+alt, alt, exp3);
                
                map<int, double> m;
                if (llrs.count(i) == 0){
                    llrs.insert(make_pair(i, m));
                }
                if (llrs.count(j) == 0){
                    llrs.insert(make_pair(j, m));
                }
                if (llrs[i].count(j) == 0){
                    llrs[i].insert(make_pair(j, 0.0));
                }
                llrs[i][j] += (ll1-ll2);
                if (doublet_rate > 0.0){
                    // Store comparisons between i and (i,j) combo and 
                    // between j and (i,j) combo
                    if (llrs[i].count(k) == 0){
                        llrs[i].insert(make_pair(k, 0.0));
                    }
                    if (llrs[j].count(k) == 0){
                        llrs[j].insert(make_pair(k, 0.0));
                    }
                    
                    llrs[i][k] += (ll1-ll3);
                    llrs[j][k] += (ll2-ll3);
                }
            }
        }
    }
    // Populate LLR table with singlet/singlet comparisons
    for (map<int, map<int, double> >::iterator x = llrs.begin(); x != llrs.end(); ++x){
        for (map<int, double>::iterator y = x->second.begin(); y != x->second.end(); ++y){
            if (y->first < n_samples){
                if (allowed_assignments.size() == 0 || 
                    (allowed_assignments.find(x->first) != allowed_assignments.end() &&
                     allowed_assignments.find(y->first) != allowed_assignments.end())){
                    tab.insert(x->first, y->first, y->second);
                }
            }
        }
    }
    if (doublet_rate > 0.0){ 
        
        // Toss out unlikely individuals (based on losing end of largest LLR
        // in the table, iteratively), until we are left with 10 individuals
        int n_target = 10;
        if (tab.n_indvs > n_target){
            bool success = tab.del(n_target);
            if (!success){
                return false;
            }
        }
        
        // If we started with more than 3 individuals and threw out enough to get down
        // below 3, give up trying to make an assignment. This will only happen if there
        // are lots of ties, which would be the result of very sparse data.

        if (tab.n_indvs < 3 && tab.n_indvs < n_samples){
            return false;
        }
        
        // Get a list of all possible double identities to consider.
        vector<int> ks;
        for (int i = 0; i < n_samples-1; ++i){
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            else if (!tab.included[i]){
                continue;
            }

            for (int j = i + 1; j < n_samples; ++j){
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                else if (!tab.included[j]){
                    continue;
                }
                int k = hap_comb_to_idx(i, j, n_samples);
                if (allowed_assignments.size() == 0 || allowed_assignments.find(k) != 
                    allowed_assignments.end()){
                    
                    ks.push_back(k);
                }
            }
        }
        
        if (ks.size() > 0){
            // put k indices in increasing order so we know what order to store comparisons
            // in the data structure
            
            sort(ks.begin(), ks.end());
            // With all possible values of k to consider, we already have all member component vs 
            // double model comparisons computed. We now need to compute all other possible
            // single vs double model comparisons, as well as double model vs double model comparisons.
            
            compute_k_comps(llrs, tab, ks, n_samples, allowed_assignments, doublet_rate);
            
        }
    }
    
    // Disallow impossible combinations. We should already have excluded disallowed k combinations
    // so only need to exclude single individuals.
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (int i = 0; i < n_samples; ++i){
            if (doublet_rate == 1 || (allowed_assignments2.size() > 0 && 
                allowed_assignments2.find(i) == allowed_assignments2.end())){
                tab.disallow(i);
            }
        }
    }

    return true;
}

/**
 * Adjust an expected allele matching rate given the ref and alt allele
 * misreading errors.
 */
double adjust_p_err(double p, double e_r, double e_a){
    return p - p*e_a + (1-p)*e_r;
}
