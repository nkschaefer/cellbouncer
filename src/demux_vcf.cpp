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
#include <sys/stat.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_io.h"
#include "demux_vcf_hts.h"

using std::cout;
using std::endl;
using namespace std;

class llr_table{
    private:
        // Set of pairs rather than map allows duplicate LLRs
        map<double, vector<pair<short, short> > > lookup_llr;
        map<double, vector<pair<short, short> > >::iterator it;
        vector<double> maxllr;
        vector<double> minllr;

    public:
        vector<bool> included;
        int n_indvs;
        llr_table(int x){
            n_indvs = 0;
            int n_elt = x + (int)round(pow(2, binom_coef_log(x, 2)));
            included.reserve(n_elt);
            maxllr.reserve(n_elt);
        };
        ~llr_table(){
            lookup_llr.clear();
            included.clear();
            maxllr.clear();
            minllr.clear();
        }
        void print(string& bc_str, vector<string>& samples){
            
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

        void print_ranges(string& barcode, vector<string>& samples){
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

        void insert(short i1, short i2, double llr){
              
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
        
        void disallow(short i){
            if (i < included.size()){
                if (included[i]){
                    maxllr[i] = 0.0;
                    minllr[i] = 0.0;
                    --n_indvs;
                }
                included[i] = false;
            }
        }

        bool del(int n_keep){

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
        
        void get_max(int& best_idx, double& best_llr){
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
        };
};

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
    vector<string>& samples,
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
    vector<string>& samples,
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
        get_kcomps_cond(kcomp, llrstmp, tab, samples.size(), samples, allowed_assignments);
        kcomp_cond.insert(make_pair(samp->first, kcomp));
    }
    // compute all single vs double comparisons
    for (int i = 0; i < samples.size(); ++i){
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
            pair<int, int> comb = idx_to_hap_comb(k, samples.size());
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
        
        pair<int, int> comb1 = idx_to_hap_comb(k1, samples.size());
        for (int kj = ki+1; kj < ks.size(); ++kj){
            int k2 = ks[kj];
            pair<int, int> comb2 = idx_to_hap_comb(k2, samples.size());
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
    vector<string>& samples,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt){
    
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
        float exp1 = error_rate_ref;
        if (y->first.second == 1){
            //exp1 = 0.5;
            exp1 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
        }
        else if (y->first.second == 2){
            exp1 = 1.0-error_rate_alt;
        }

        for (map<pair<int, int>, pair<float, float> >::iterator z = 
            y->second.begin(); z != y->second.end(); ++z){
            
            if (allowed_assignments.size() > 0 && allowed_assignments.find(z->first.first) ==
                allowed_assignments.end()){
                continue;
            } 

            // If same site type, we can't distinguish between the
            // two individuals from this piece of information
            if (z->first.first != -1 && y->first.second != z->first.second){
                
                // Set default expectation for indv2
                float exp2 = error_rate_ref;
                if (z->first.second == 1){
                    //exp2 = 0.5;
                    exp2 = 0.5*(1.0 - error_rate_alt + error_rate_ref);
                }
                else if (z->first.second == 2){
                    exp2 = 1.0-error_rate_alt;
                }
                
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

                int i = y->first.first;
                int j = z->first.first;
                int k = hap_comb_to_idx(i, j, samples.size());

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
            if (y->first < samples.size()){
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

        if (tab.n_indvs < 3 && tab.n_indvs < samples.size()){
            return false;
        }
        
        // Get a list of all possible double identities to consider.
        vector<int> ks;
        for (int i = 0; i < samples.size()-1; ++i){
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            else if (!tab.included[i]){
                continue;
            }

            for (int j = i + 1; j < samples.size(); ++j){
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                else if (!tab.included[j]){
                    continue;
                }
                int k = hap_comb_to_idx(i, j, samples.size());
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
            
            compute_k_comps(llrs, tab, ks, samples, allowed_assignments, doublet_rate);
            
        }
    }
    
    // Disallow impossible combinations. We should already have excluded disallowed k combinations
    // so only need to exclude single individuals.
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (int i = 0; i < samples.size(); ++i){
            if (doublet_rate == 1 || (allowed_assignments2.size() > 0 && 
                allowed_assignments2.find(i) == allowed_assignments2.end())){
                tab.disallow(i);
            }
        }
    }

    return true;
}

/**
 * Given genome-wide counts of different types of alleles, determines the 
 * most likely identity of each cell and stores in the data structure.
 */
void assign_ids(robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    set<int>& allowed_assignments,
    set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt){
    
    assignments.clear();
    assignments_llr.clear();

    // Left in for debugging purposes only: if print_llrs is true,
    // the user will be prompted for a cell barcode of interest.
    // Once found, the log likelihood ratio table for the selected
    // cell barcode will be printed to stdout and the program will
    // exit.
    bool print_llrs = false;
    unsigned long searchbc;
    if (print_llrs){
        string searchbc_str;
        cerr << "Enter barcode: ";
        cin >> searchbc_str;
        if (searchbc_str.length() != 16){
            cerr << "ERROR: input barcode length is not 16";
            exit(1);
        }
        cerr << "Searching for " << searchbc_str;
        cerr << "\n"; 
        bc searchbc_bin;
        str2bc(searchbc_str.c_str(), searchbc_bin);
        searchbc = searchbc_bin.to_ulong();
    }
    
    for (robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >::iterator x = 
            indv_allelecounts.begin(); x != indv_allelecounts.end(); ++x){
        
        if (print_llrs && x->first != searchbc){
            continue;
        }
        
        // Get a table of log likelihood ratios between every possible
        // pair of identities
        map<int, map<int, double> > llrs;
        llr_table tab(samples.size());
        bool success = populate_llr_table(x->second, llrs, tab, samples, allowed_assignments, 
            allowed_assignments2,doublet_rate, error_rate_ref, error_rate_alt);
        
        // Debugging only: print this table and quit
        if (print_llrs){
            string bc_str = bc2str(x->first);
            tab.print(bc_str, samples);
            
            int x;
            double y;
            tab.get_max(x, y);
            fprintf(stderr, "%s %f\n", idx2name(x, samples).c_str(), y);
            exit(0);
        
        }    
        
        // Now that we have a table of all pairwise LLRs between assignments for this cell,
        // choose the best identity.
        
        // This is done by iteratively finding the largest LLR between identities and eliminating
        // the least likely one, until only two identities are left. The final LLR is the LLR between
        // the best and second best assignment.
        double llr_final = 0.0;
        int assn = -1;
        if (success){
            tab.get_max(assn, llr_final);
        }
        
        // Only store information if an assignment has been made (don't accept equal
        // likelihood of two choices)
        if (llr_final > 0.0){
            assignments.emplace(x->first, assn);
            assignments_llr.emplace(x->first, llr_final);
        }
    }            
}

/**
 * Log likelihood function for computing error rates, for use by
 * multivar_ml_solver.
 */
double ll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p0 = data_d.at("exp");
    double e_r = params[0];
    double e_a = params[1];
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    double ll = dbinom(n, k, p)/log2(exp(1));
    return ll;
}

/**
 * Derivative of log likelihood function (wrt error rates), for use by
 * multivar_ml_solver, for re-estimating reference and alt allele
 * misreading error rates.
 */
void dll_err(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p0 = data_d.at("exp");
    double e_r = params[0];
    double e_a = params[1];
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    double ll = dbinom(n, k, p)/log2(exp(1));
    
    double dy_dp = (k-n*p)/(p-p*p);
    double dp_de_a = -p0;
    double dp_de_r = 1.0 - p0;
    results[0] = dy_dp * dp_de_r;
    results[1] = dy_dp * dp_de_a;
}

/**
 * After assigning all cells to identities, use the assignments
 * (weighted by log likelihood ratio of best assignment)
 * to re-infer error rates (rate of reading ref alleles as
 * alt, and rate of reading alt alleles as ref).
 *
 * Then we will re-assign identities using the newly-calculated
 * error rates.
 */
pair<double, double> infer_error_rates(robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
    map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    int n_samples,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    double error_ref,
    double error_alt,
    double error_sigma,
    vector<string>& samples){
   
    vector<double> n;
    vector<double> k;
    vector<double> expected;
    vector<double> weights_llr;

    pair<int, int> nullkey = make_pair(-1, -1);
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end();
        ++a){
        
        double weight = assn_llr[a->first];
        bool is_combo = false;
        pair<int, int> combo;
        if (a->second >= n_samples){
            is_combo = true;
            combo = idx_to_hap_comb(a->second, n_samples);
        }
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator x = 
            indv_allelecounts[a->first].begin(); x != indv_allelecounts[a->first].end(); ++x){
            if ((!is_combo && x->first.first == a->second) || (is_combo && x->first.first == combo.first)){
                if (is_combo){
                    for (map<pair<int, int>, pair<float, float> >::iterator y = x->second.begin(); y !=
                        x->second.end(); ++y){
                        if (y->first.first == combo.second){
                            double this_expected = (double)(x->first.second + y->first.second)/4.0;
                            expected.push_back(this_expected);
                            n.push_back(y->second.first + y->second.second);
                            k.push_back(y->second.second);
                            weights_llr.push_back(weight);
                        }
                    }
                }
                else{
                    double this_expected = (double)x->first.second / 2.0;
                    expected.push_back(this_expected);
                    n.push_back(x->second[nullkey].first + x->second[nullkey].second);
                    k.push_back(x->second[nullkey].second);
                    weights_llr.push_back(weight);
                }
            }
        }
    }

    optimML::multivar_ml_solver solver({error_ref, error_alt}, ll_err, dll_err);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("exp", expected);
    solver.add_weights(weights_llr);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_normal_prior(0, error_ref, error_sigma, 0.0, 1.0);
    solver.add_normal_prior(1, error_alt, error_sigma, 0.0, 1.0);
    solver.solve();
    
    return make_pair(solver.results[0], solver.results[1]);
}

/**
 * Set id2 == -1 for second ID = all other IDs in the set
 */
double mannwhitney_llr(int id1, 
    int id2, 
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    double n1 = 0.0;
    double n2 = 0.0;

    vector<pair<double, int> > llrs;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (a->second == id1){
                llrs.push_back(make_pair(assn_llr[a->first], id1));
                n1++;
            }
            else if (id2 == -1 || a->second == id2){
                llrs.push_back(make_pair(assn_llr[a->first], id2));
                n2++;
            }
        }
    }
    
    if (n1 == 0){
        return 0.0;
    }
    else if (n2 == 0){
        return 1.0;
    }

    sort(llrs.begin(), llrs.end());
    
    // Assign ranks
    vector<double> ranks;
    int rank = 1;
    for (int i = 0; i < llrs.size(); ++i){
        ranks.push_back((double)rank);
        rank++;
    }

    // Deal with ties in ranks
    double prevllr = 0.0;
    double ranksum = 0.0;
    int ranktot = 0;
    bool ties = false;

    vector<double> nties;

    for (int i = 0; i < llrs.size(); ++i){
        if (llrs[i].first == prevllr){
            ranksum += ranks[i];
            ranktot++;
        }
        else{
            if (ranktot > 1){
                nties.push_back((double)ranktot);
                ties = true;
                double rankmean = ranksum / (double)ranktot;
                for (int j = i - 1; j >= i - 1 - (ranktot-1); --j){
                    ranks[j] = rankmean;
                }
                ranksum = 0.0;
                ranktot = 0;
            }
            else{
                if (i > 0){
                    nties.push_back(0);
                    ranks[i-1] = ranksum;
                }
                ranksum = 0;
                ranktot = 0;
            }
            ranksum += ranks[i];
            ranktot++;

        }
        prevllr = llrs[i].first;    
    }
    // Handle last one.
    if (ranktot > 1){
        ties = true;
        nties.push_back(ranktot);
        double rankmean = ranksum / (double)ranktot;
        for (int j = llrs.size()-1; j >= llrs.size()-1 - (ranktot-1); --j){
            ranks[j] = rankmean;
        }
    }
    else{
        nties.push_back(0);
        ranks[llrs.size()-1] = ranksum;
    }
    
    double sum_id1 = 0;
    double sum_id2 = 0;
    for (int i = 0; i < ranks.size(); ++i){
        if (llrs[i].second == id1){
            sum_id1 += ranks[i];
        }
        else{
            sum_id2 += ranks[i];
        }
    } 
    
    double U1 = sum_id1 - (n1*(n1+1))/2.0;
    double U2 = sum_id2 - (n2*(n2+1))/2.0;
    double m_u = (n1*n2)/2.0;
    double sigma_u = sqrt((n1*n2*(n1+n2 +1))/(12.0));
    if (ties){
        double tsum = 0;
        for (int i = 0; i < nties.size(); ++i){
            if (nties[i] > 0){
                tsum += (pow(nties[i], 3) - nties[i]);
            }
        }
        double term1 = (n1*n2*(n1+n2+1))/12;
        double term2 = (n1*n2*tsum)/(12*(n1+n2)*(n1+n2-1));
        sigma_u = sqrt(term1-term2);
    }
    
    // cnorm(U1, m_u, sigma_u) is the same as 
    // 1.0 - cnorm(U2, m_u, sigma_u)

    // For testing alt: U1 < U2, want cnorm(U1, m_u, sigma_u)
    // For testing alt: U1 > U1, want 1.0 - cnorm(U1, m_u, sigma_u)

    if (n1 < 3 || n2 < 3){
        // Can't get a reliable sigma here.
        if (U1 < m_u){
            return 0.0;
        }
        else{
            return 1.0;
        }
    }

    // Test whether id1 < id2
    double p = pnorm(U1, m_u, sigma_u);
    return p;
}

/**
 * Does some QC on assignments - for each ID in the assignment file,
 * checks for significantly lower LLR distribution than the rest of
 * the data set. Also checks for significantly lower numbers of cells.
 */
void id_qc(robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<int, double>& pois_p,
    map<int, double>& mannwhitney_p){
    
    // Get total num cells for each ID 
    map<int, int> idsizes;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (assn_llr[a->first] > 0){
            if (idsizes.count(a->second) == 0){
                idsizes.insert(make_pair(a->second, 0));
            }
            idsizes[a->second]++;
        }
    }

    for (map<int, int>::iterator ids = idsizes.begin(); ids != idsizes.end(); ++ids){
        double mean_othersize = 0.0;
        double mean_othersize_tot = 0.0;
        for (map<int, int>::iterator ids2 = idsizes.begin(); ids2 != idsizes.end(); ++ids2){
            if (ids2->first != ids->first){
                mean_othersize += (double)ids2->second;
                mean_othersize_tot++;
            }
            
        }
        double p1 = ppois(ids->second, mean_othersize_tot);
        pois_p.insert(make_pair(ids->first, p1));
        double p2 = mannwhitney_llr(ids->first, -1, assn, assn_llr);
        mannwhitney_p.insert(make_pair(ids->first, p2));
    } 
}

bool filter_identities(robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr, 
    int n_samples,
    set<int>& allowed_ids, 
    set<int>& allowed_ids2){
    
    bool removed_ids = false;
    
    // Get number of cells per ID
    map<int, int> cells_per_id;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (cells_per_id.count(a->second) == 0){
            cells_per_id.insert(make_pair(a->second, 0));
        }
        cells_per_id[a->second]++;
    }

    // Get a list of identities not in the original filtered list, but added because
    // they're singlet components of allowed doublets
    set<int> candidate_ids;
    for (set<int>::iterator a = allowed_ids.begin(); a != allowed_ids.end(); ++a){
        if (allowed_ids2.find(*a) == allowed_ids2.end()){
            
            // Check whether this individual has significantly fewer cells than 
            // all "parent" individuals (doublets including this individual)
            vector<double> p_ncell;
            
            // Check whether cells assigned to this individual have significantly
            // lower LLRs than cells assigned to "parent" individuals
            vector<double> p_llr;
            
            for (set<int>::iterator a2 = allowed_ids2.begin(); a2 != allowed_ids2.end(); ++a2){
                if (*a2 >= n_samples){
                    pair<int, int> combo = idx_to_hap_comb(*a2, n_samples);
                    if (combo.first == *a || combo.second == *a){
                        double p1 = ppois((double)cells_per_id[*a], (double)cells_per_id[*a2]);
                        double p2 = mannwhitney_llr(*a, *a2, assn, assn_llr);
                        p_ncell.push_back(p1);
                        p_llr.push_back(p2);
                    }
                }
            }
            bool all_p_signif = true;
            for (int i = 0; i < p_ncell.size(); ++i){
                if (p_ncell[i] > 0.05 || p_llr[i] > 0.05){
                    all_p_signif = true;
                }
            }
            if (all_p_signif){
                // Remove a.
                removed_ids = true;
            }
            else{
                // Safe to keep.
                allowed_ids2.insert(*a);
            }
        }
    }
    return removed_ids;
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_vcf [OPTIONS]\n");
    fprintf(stderr, "Attempts to demultiplex cells in a BAM file based on genotype data \n");
    fprintf(stderr, "for source cell lines provided via VCF. Can infer both single-\n");
    fprintf(stderr, "individual cells and cells better modeled as a combination of two \n");
    fprintf(stderr, "individuals.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants. Only biallelic SNPs \n");
    fprintf(stderr, "       will be considered, and phasing will be ignored.\n");
    fprintf(stderr, "    --output_prefix -o Base name for output files\n");
    fprintf(stderr, "===== RECOMMENDED =====\n");
    fprintf(stderr, "    --barcodes -B To consider only barcodes that have passed \n");
    fprintf(stderr, "       filtering and represent true cells, provide a filtered list \n");
    fprintf(stderr, "       of cell barcodes here. If you mapped the data using \n");
    fprintf(stderr, "       cellranger, this file will be located in \n");
    fprintf(stderr, "       outs/filtered_feature_bc_matrix/barcodes.tsv or \n");
    fprintf(stderr, "       outs/filtered_feature_bc_matrix/barcodes.tsv.gz.\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "----- Algorithm options -----\n");
    fprintf(stderr, "    --doublet_rate -D Prior probability of a cell being a mixture \n");
    fprintf(stderr, "       of two individuals rather than a single individual. \n");
    fprintf(stderr, "       default = 0.5 (set to 0 to disable checking for doublet \n");
    fprintf(stderr, "       combinations and to 1 to disable checking for single indvs\n");
    fprintf(stderr, "    --error_ref -e Initial/prior error rate for misreading ref alleles as alt\n");
    fprintf(stderr, "       (default = 0.005)\n");
    fprintf(stderr, "    --error_alt -E Initial/prior error rate for misreading alt alleles as ref\n");
    fprintf(stderr, "       (default = 0.005)\n");
    fprintf(stderr, "    --error_sigma -s After the first round of assignments, error rates will be\n");
    fprintf(stderr, "       re-estimated from the assignments and used to make a second round of\n");
    fprintf(stderr, "       assignments using these new (maximum a posteriori) error rates. This\n");
    fprintf(stderr, "       estimation will use the initial values as means of prior distributions,\n");
    fprintf(stderr, "       which will be truncated normal distributions on (0, 1) with standard\n");
    fprintf(stderr, "       deviation equal to this parameter. A larger value gives more weight to\n");
    fprintf(stderr, "       the data and less to the prior. Default = 0.1\n");
    fprintf(stderr, "    --qual -q Minimum variant quality to consider (default 50)\n");
    fprintf(stderr, "    --ids -i If the VCF file contains individuals that you do not\n");
    fprintf(stderr, "       expect to see in your sample, specify individuals to include here.\n");
    fprintf(stderr, "       These identities will be considered, as well as all possible doublet\n");
    fprintf(stderr, "       combinations of them. Expects a text file, with one individual ID per\n");
    fprintf(stderr, "       line, matching individual names in the VCF.\n");
    fprintf(stderr, "    --ids_doublet -I Similar to --ids/-i argument above, but allows control\n");
    fprintf(stderr, "       over which doublet identities are considered. Here, you can specify\n");
    fprintf(stderr, "       individuals and combinations of two individuals to allow. Doublet\n");
    fprintf(stderr, "       combinations not specified in this file will not be considered.\n");
    fprintf(stderr, "       Single individuals involved in doublet combinations specified here\n");
    fprintf(stderr ,"       but not explicitly listed in the file will still be considered.\n");
    fprintf(stderr, "       Names of individuals must match those in the VCF, and combinations\n");
    fprintf(stderr, "       of two individuals can be specified by giving both individual names\n");
    fprintf(stderr, "       separated by \"+\", with names in either order.\n");
    fprintf(stderr, "----- I/O options -----\n");
    print_libname_help();
    fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    fprintf(stderr, "    --disable_conditional -f By default, the program will compute\n");
    fprintf(stderr, "       expected alt allele matching probabilities for SNPs of each type\n");
    fprintf(stderr, "       (homozygous ref, het, or homozygous alt in each individual), conditional\n");
    fprintf(stderr, "       on the cell's identity being each individual. This will be spilled\n");
    fprintf(stderr, "       to disk and can be used by quant_contam. If you do not plan on inferring\n");
    fprintf(stderr, "       ambient RNA contamination, setting this option can slightly speed things\n");
    fprintf(stderr, "       up by omitting this step.\n");
    fprintf(stderr, "   --dump_conditional -F If you ran once with -f (to disable computing\n");
    fprintf(stderr, "       alt allele matching probabilities conditional on each identity), and you\n");
    fprintf(stderr, "       now wish to run quant_contam, re-run this program with this option enabled\n");
    fprintf(stderr, "       to load the VCF, compute these values, and write them to disk. Use the\n");
    fprintf(stderr, "       same output prefix you used on the previous run.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"index_jump", no_argument, 0, 'j'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"qual", required_argument, 0, 'q'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", no_argument, 0, 'E'},
       {"error_sigma", required_argument, 0, 's'},
       {"disable_conditional", no_argument, 0, 'f'},
       {"dump_conditional", no_argument, 0, 'F'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string vcf_file = "";
    bool cell_barcode = false;
    string cell_barcode_file = "";
    bool stream = true; // Opposite of index-jumping
    string output_prefix = "";
    int vq = 50;
    string idfile;
    string idfile_doublet;
    bool idfile_given = false;
    bool idfile_doublet_given = false;
    double doublet_rate = 0.5;
    string barcode_group = "";
    double error_ref = 0.005;
    double error_alt = 0.005;
    double error_sigma = 0.1;
    
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;

    bool disable_conditional = false;
    bool dump_conditional = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:B:i:I:q:D:n:e:E:p:s:fFCSUjh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'b':
                bamfile = optarg;
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'n':
                barcode_group = optarg;
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
            case 'B':
                cell_barcode = true;
                cell_barcode_file = optarg;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'i':
                idfile_given = true;
                idfile = optarg;
                break;
            case 'I':
                idfile_doublet_given = true;
                idfile_doublet = optarg;
                break;
            case 'q':
                vq = atoi(optarg);
                break;
            case 'j':
                stream = false;
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 's':
                error_sigma = atof(optarg);
                break;
            case 'f':
                disable_conditional = true;
                break;
            case 'F':
                dump_conditional = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (vq < 0){
        fprintf(stderr, "ERROR: variant quality must be a positive integer (or 0 for no filter)\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        exit(1);
    }
    if (doublet_rate < 0 || doublet_rate > 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1, inclusive.\n");
        exit(1);
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i/-I is allowed.\n");
        exit(1);
    }
    if (disable_conditional && dump_conditional){
        fprintf(stderr, "ERROR: only one of -f/-F is allowed.\n");
        exit(1);
    }
    
    // Init BAM reader
    bam_reader reader = bam_reader();

    // Decide whether we will be loading counts or computing them
    bool load_counts = false;
    string countsfilename = output_prefix + ".counts";
    if (!dump_conditional && file_exists(countsfilename)){
        load_counts = true;
    }
    else if (!dump_conditional){
        if (bamfile.length() == 0){
            fprintf(stderr, "ERROR: bam file (--bam) required\n");
            exit(1);
        }    
        // We will be actually using the BAM reader, so we need to initialize properly.
        reader.set_file(bamfile);
    }

    // Store the names of all individuals in the VCF
    vector<string> samples;
    
    // Store data about SNPs
    map<int, map<int, var> > snpdat;
    
    int nsnps;
    
    bool samples_from_vcf = false;

    if (load_counts){
        // Load variant data from VCF file or from samples file, if it was written
        string samplesfile = output_prefix + ".samples";
        if (file_exists(samplesfile)){
           load_samples(samplesfile, samples); 
        }
        else{
            // Get samples from VCF header as a fallback
            if (vcf_file == ""){
                fprintf(stderr, "ERROR: vcf file is required\n");
                exit(1);
            }
            fprintf(stderr, "Reading VCF/BCF header...\n");
            nsnps = read_vcf(vcf_file, reader, samples, snpdat, vq, false, true);
            samples_from_vcf = true;
        }
    }
    else{
        if (vcf_file == ""){
            fprintf(stderr, "ERROR: vcf file is required\n");
            exit(1);
        }
        fprintf(stderr, "Loading variants from VCF/BCF...\n");
        nsnps = read_vcf(vcf_file, reader, samples, snpdat, vq, dump_conditional, load_counts);
        // Compute conditional matching fractions for quant_contam (can disable)
        if (!disable_conditional){
            map<pair<int, int>, map<int, float> > conditional_match_fracs;
            fprintf(stderr, "Computing conditional alt allele probabilities...\n");
            get_conditional_match_fracs(snpdat, conditional_match_fracs, samples.size());
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            dump_exp_fracs(outf, conditional_match_fracs);
            fclose(outf); 
            if (dump_conditional){
                // Nothing left to do.
                return 0;
            }
        }
        samples_from_vcf = true;
    }

    set<int> allowed_ids;
    set<int> allowed_ids2;

    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile.c_str());
        }
    }
    if (idfile_doublet_given){
        parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile_doublet.c_str());
        }
    }
    
    if (samples_from_vcf){
        // Store these to disk in case we run ambient RNA contamination finding later
        string samplesfile = output_prefix + ".samples"; 
        write_samples(samplesfile, samples);
    }
    
    set<unsigned long> cell_barcodes;
    if (cell_barcode){
        parse_barcode_file(cell_barcode_file, cell_barcodes);
        if (cell_barcodes.size() == 0){
            // Did not read barcodes correctly
            fprintf(stderr, "ERROR reading cell barcode list\n");
            exit(1);
        }
    }

    // Data structure to store allele counts at SNPs of each possible type.
    // SNPs are defined by their allelic state in each pair of 2 individuals.
    
    // In other words, one category is all SNPs at which individual A is 
    // homozygous for the reference allele and individual B is heterozygous.
    // Another category is all SNPs at which individual C is homozygous for
    // the alt allele and individual D is homozygous for the reference allele.
    
    // For each cell, we sum reference and alt allele counts across all SNPs
    // falling into each category.
    
    // Here, cell barcodes are represented as unsigned long (see bc_hash.cpp in 
    // htswrapper library)
    
    // Both pair<int, int> are keys of (individual ID, nalt) where nalt is 
    // number of alt alleles for that individual in SNPs of that type
    
    // Each combination is stored only once - lower index individual always
    // comes first

    // The final two numbers are ref and alt allele counts -- where each count
    // is actually Probability(mapping of read is correct), determined by
    // map quality

    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;

    // Store counts for currently-tracked SNPs
    map<int, map<int, robin_hood::unordered_map<unsigned long, 
        pair<float, float> > > > varcounts_site;
    
    robin_hood::unordered_map<unsigned long, map<int, float> > fracs;
    
    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;

    // Load pre-computed allele counts from file, if it already exists
    if (load_counts){
        // Figure out the appropriate file name from the previous run and
        // load the counts, instead of reading through the BAM file.
        fprintf(stderr, "Loading counts...\n");
        load_counts_from_file(indv_allelecounts, samples, countsfilename, allowed_ids);
    }
    else{

        fprintf(stderr, "Counting alleles in BAM file...\n");

        // initialize the BAM reader
        reader.set_file(bamfile);
        
        // retrieve cell barcodes
        reader.set_cb();
        
        // Get a mapping of chromosome names to internal numeric IDs in the 
        // BAM file. This is necessary to reconcile how HTSLib might represent
        // the same chromosome name in the BAM vs the VCF file.
        map<string, int> chrom2tid = reader.get_seq2tid();
        map<int, string> tid2chrom;
        for (map<string, int>::iterator ct = chrom2tid.begin(); ct != chrom2tid.end(); ++ct){
            tid2chrom.insert(make_pair(ct->second, ct->first));
        } 
        
        int nsnp_processed = 0;
        if (stream){

            // Read through entire BAM file and look for informative SNPs along the way
            // (default setting, appropriate for large numbers of SNPs).
            int curtid = -1;
            
            map<int, var>::iterator cursnp;
            while (reader.next()){
                if (curtid != reader.tid()){
                    // Started a new chromosome
                    if (curtid != -1){
                        while (cursnp != snpdat[curtid].end()){
                            if (snpdat[curtid].count(cursnp->first) > 0){        
                                dump_vcs_counts(varcounts_site[curtid][cursnp->first], 
                                    indv_allelecounts,
                                    cursnp->second, samples.size());
                            }
                            ++nsnp_processed;
                            ++cursnp;
                        }
                    }
                    cursnp = snpdat[reader.tid()].begin();
                    curtid = reader.tid();
                }
                // Advance to position within cur read
                while (cursnp != snpdat[reader.tid()].end() && 
                    cursnp->first < reader.reference_start){
                    
                    if (varcounts_site[reader.tid()].count(cursnp->first) > 0){
                        dump_vcs_counts(varcounts_site[reader.tid()][cursnp->first], 
                            indv_allelecounts, 
                            cursnp->second, samples.size());
                    }
                    ++nsnp_processed;
                    ++cursnp;
                }
                // Create a second iterator to look ahead for any additional SNPs 
                // within the current read
                map<int, var>::iterator cursnp2 = cursnp;
                while (cursnp2 != snpdat[reader.tid()].end() && 
                    cursnp2->first >= reader.reference_start && 
                    cursnp2->first <= reader.reference_end){   
                    process_bam_record(reader, cursnp2->first, cursnp2->second,
                        varcounts_site[reader.tid()], cell_barcode, cell_barcodes);
                    ++cursnp2;
                }
                if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                    fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, nsnps); 
                    last_print = nsnp_processed;
                }
            }
            
        }
        else{
            // Visit each SNP and index-jump in the BAM to it.
            for (map<int, map<int, var> >::iterator curchrom = snpdat.begin();
                curchrom != snpdat.end(); ++curchrom){

                int tid = curchrom->first;
                for (map<int, var>::iterator cursnp = curchrom->second.begin();
                    cursnp != curchrom->second.end(); ++cursnp){
                    int pos = cursnp->first;
                    
                    reader.set_query_site(tid, pos); 
                    
                    while(reader.next()){
                        process_bam_record(reader, cursnp->first, cursnp->second, 
                            varcounts_site[tid], cell_barcode, cell_barcodes);
                    }
                    
                    dump_vcs_counts(varcounts_site[tid][cursnp->first], 
                        indv_allelecounts, cursnp->second, samples.size());
                    
                    ++nsnp_processed;        
                
                    if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                        fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, 
                            nsnps); 
                        last_print = nsnp_processed;
                    }
                }
            }
        }
        fprintf(stderr, "Processed %d of %d SNPs\n", nsnp_processed, nsnps);
        
        // Write the data just compiled to disk.
        string fname = output_prefix + ".counts";
        //FILE* outf = fopen(fname.c_str(), "w");   
        gzFile outf = gzopen(fname.c_str(), "w");
        fprintf(stderr, "Writing allele counts to disk...\n");
        dump_cellcounts(outf, indv_allelecounts, samples);
        fprintf(stderr, "Done\n");
        //fclose(outf);
        gzclose(outf);
    }
        
    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    fprintf(stderr, "Finding likeliest identities of cells...\n");
    
    // Get assignments of cell barcodes
    
    assign_ids(indv_allelecounts, samples, assn, assn_llr, 
        allowed_ids, allowed_ids, doublet_rate, error_ref, error_alt);

    // Now, from these assignments, compute the likeliest error rate, weighting by
    //  log likelihood ratio of assignment.
    
    fprintf(stderr, "Finding likeliest alt/ref switch error rates...\n"); 
    pair<double, double> err_new = infer_error_rates(indv_allelecounts, samples.size(),
        assn, assn_llr, error_ref, error_alt, error_sigma, samples);
   
    double error_ref_posterior = err_new.first;
    double error_alt_posterior = err_new.second; 
    
    fprintf(stderr, "Posterior error rates:\n");
    fprintf(stderr, "\tref mismatch: %f\n", error_ref_posterior);
    fprintf(stderr, "\talt mismatch: %f\n", error_alt_posterior);
    
    // Re-assign individuals using posterior error rates
    fprintf(stderr, "Re-inferring identities of cells...\n");
    assign_ids(indv_allelecounts, samples, assn, assn_llr,
        allowed_ids, allowed_ids, doublet_rate, error_ref_posterior, error_alt_posterior);
    
    if (idfile_doublet_given){
        // The user gave an allowed list of specific doublet combinations, and we included
        // all possible singlets from the allowable doublets in the first round. If any of 
        // these singlet identities turned out to be uncommon, we will assume the user was
        // right and that those singlet identities truly do not exist in the pool -- and
        // perform one more round of assignments without them.
        if (allowed_ids.size() > allowed_ids2.size()){
            bool altered = filter_identities(assn, assn_llr, samples.size(), allowed_ids, 
                allowed_ids2);

            if (altered){
                fprintf(stderr, "Re-inferring with unlikely singlet identities removed...\n");
                assign_ids(indv_allelecounts, samples, assn, assn_llr,
                    allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior,
                    error_alt_posterior); 
            }
        }
    }
    
    map<int, double> p_ncell;
    map<int, double> p_llr;
    id_qc(assn, assn_llr, p_ncell, p_llr);

    // Write these best assignments to disk
    {
        string fname = output_prefix + ".assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing cell-individual assignments to disk...\n");
        dump_assignments(outf, assn, assn_llr, samples, barcode_group, cellranger, seurat, underscore);
        fclose(outf);
    }
    
    // Create summary file
    {
        string fname = output_prefix + ".summary";
        FILE* outf = fopen(fname.c_str(), "w");
        write_summary(outf, output_prefix, assn, samples, error_ref,
            error_alt, error_sigma, error_ref_posterior,
            error_alt_posterior, vcf_file, vq, doublet_rate,
            p_ncell, p_llr);
        fclose(outf);
    }
    
    return 0;
}
