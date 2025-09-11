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
#include <optimML/brent.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_llr.h"
#include "ambient_rna.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Functions related to inferring ambient RNA profile within cells as part of
 * demux_vcf.
 *
 * These are all built on the model where every observed reference and alt 
 * alt allele count for reads covering SNPs of a given category (i.e. SNPs
 * that are known from the VCF to be heterozygous in individual A) are modeled
 * as draws from a binomial distribution with parameters n, k, and p as follows:
 *
 * Let  n = ref + alt count
 *      k = alt count
 *      p = (1-c)*p_e + c*p_c
 * Where
 *      c = contamination rate for the cell
 *      p_e = expected fraction of alt reads given the category of SNPs
 *          (i.e. for cells from individual A at SNPs heterozygous in individual A,
 *          p_e = 0.5)
 *      p_c expected fraction of alt reads in ambient RNA given the category
 *          of SNPs
 *
 * To avoid expectations of 0 and 1 (i.e. where an individual is homozygous for
 *  reference alleles or homozygous for alt alleles), we also model to error rates
 *      e_r = Expected rate at which reference alleles are misread as alt
 *      e_a = Expected rate at which alt alleles are misread as reference
 * Then we adjust p_e as follows:
 *      p_e_adjusted = p_e - p_e*e_a + (1-p_e)*e_r
 *
 * These two parameters can be set by the user, but should be kept at a low value
 *  to reflect mostly sequencing error. Alternatively, if the variant calls are known
 *  to be very noisy, these can be set to higher values.
 *  At the end of the run, the user will be shown the residual e_r and e_a after 
 *  modeling contamination -- these should be low numbers unless variant calls are
 *  unreliable.
 */

// Define a value to bump 0 and 1 binomial expectations away from these
// numbers
double binom_0 = 1e-6;

// ===== Utility functions ===== 

/**
 * Integrate binomial log likelihood wrt p_c (expected alt allele match rate from contamination)
 * 
 * * NOT USED *
 */
double integral_ll_pc(double n, double k, double c, double p_e, double p_c){
    double x = c*p_c*(binom_coef_log(n,k)/log2(exp(1)) - n);
    double y = (k-n)*((c-1)*p_e - c*p_c + 1);
    double z = log((c-1)*p_e - c*p_c + 1);
    double zz = k*(-c*p_e + c*p_c + p_e)*log(-c*p_e + c*p_c + p_e);
    return (1/c)*(x + y*z + zz);
}

/**
 * Derivative wrt contam rate (c) of the integral of binomial log likelihood wrt p_c
 * (expected alt allele match rate from contamination)
 *
 * This allows us to find the likeliest c, over the range of all possible p_c.
 * 
 * * NOT USED *
 */
double d_integral_ll_pc_dc(double n, double k, double c, double p_e, double p_c){
    double x = (p_e-1)*(k-n)*log((c-1)*p_e - c*p_c + 1);
    double y = k*p_e*log(-c*p_e + c*p_c + p_e);
    double z = c*n*(p_c-p_e);
    return (x - y + z)/(c*c);
}

/**
 * * NOT USED *
 */
double integral_ll_c(double n, double k, double c, double p_e, double p_c){
    double denom = p_e - p_c;
    double x = -c*(p_e - p_c)*(n - binom_coef_log(n,k)/log2(exp(1)));
    double y = (k-n)*((c-1)*p_e - c*p_c + 1)*log((c-1)*p_e - c*p_c + 1);
    double z = k*((c-1)*p_e - c*p_c)*log(-c*p_e + c*p_c + p_e);
    
    return (x - y + z)/denom;
}

/**
 * * NOT USED *
 */
double d_integral_ll_c_dpc(double n, double k, double c, double p_e, double p_c){
    double denom = pow(p_e-p_c, 2);
    double x = (k-n)*(c*(p_e-p_c) + (p_e-1)*log((c-1)*p_e - c*p_c + 1));
    double y = -(k*(c*(p_e - p_c) + p_e*log(-c*p_e + c*p_c + p_e)));
    
    return (x+y)/denom;
}

// ===== Functions to use for Brent's root finding method - for maximizing
// ===== univariate log likelihood functions

/**
 * Log likelihood when estimating c
 */
double ll_c(double c, const map<string, double >& data_d, 
    const map<string, int >& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    return logbinom(n, k, binom_p);
}

/**
 * Derivative of log likelihood wrt c
 */
double dll_dc(double c, const map<string, double >& data_d, 
    const map<string, int >& data_i){

    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p * binom_p);
    double dp_dc = p_c - p_e;
    return dy_dp * dp_dc;
}

/**
 * Second derivative of log likelihood wrt c
 */
double d2ll_dc2(double c, const map<string, double >& data_d, 
    const map<string, int>& data_i){

    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double p_c = data_d.at("p_c");
    double binom_p = (1.0 - c)*p_e + c*p_c;
    double dp_dc = p_c - p_e;
    double d2y_dp2 = (k*(2*binom_p - 1) - n*binom_p*binom_p)/(pow(binom_p-1,2)*binom_p*binom_p);
    return d2y_dp2 * dp_dc * dp_dc;
}

// ===== Functions to use for multivar_ml_solver -- for optimizing
// ===== multivariate log likelihood functions using BFGS

/**
 * Log likelihood for estimating best set of p_c parameters
 * in ambient RNA
 */
double ll_ambmu(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = data_d.at("c");
    double binom_p = (1.0-c) * p_e + c*p_c;
    return logbinom(n, k, binom_p); 
}

/**
 * First derivative of log likelihood for estimating best set
 * of p_c parameters in ambient RNA
 */
void dll_ambmu(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = data_d.at("c");
    double binom_p = (1.0-c) * p_e + c*p_c;

    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    results[amb_idx] = dy_dp * c;
}

/**
 * Log likelihood function for finding best e_r & e_a
 */
double ll_err_rates(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    double e_r = params[0];
    double e_a = params[1];
    double p_c = data_d.at("p_c");
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double p_e_a = adjust_p_err(p_e, e_r, e_a);
    double binom_p = (1-c)*p_e_a + c*p_c;
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0-DBL_MIN*1e6){
        binom_p = 1.0-DBL_MIN*1e6;
    }
    
    return logbinom(n, k, binom_p);
}

/**
 * First derivative log likelihood function for finding best e_r & e_a
 */
void dll_err_rates(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double e_r = params[0];
    double e_a = params[1];
    double p_c = data_d.at("p_c");
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double p_e_a = adjust_p_err(p_e, e_r, e_a);
    double binom_p = (1-c)*p_e_a + c*p_c;
    if (binom_p < DBL_MIN*1e6){
        binom_p = DBL_MIN*1e6;
    }
    else if (binom_p > 1.0-DBL_MIN*1e6){
        binom_p = 1.0-DBL_MIN*1e6;
    }
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    
    results[0] += dy_dp * (1 - p_e + c*(p_e - p_c));
    results[1] += dy_dp * (c*(p_e - p_c) - p_e);
    
}

/**
 * optimML-format log likelihood function for updating the 
 * ambient RNA profile mixture proportions.
 */
double ll_amb_prof_mixture(const vector<double>& params, 
    const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double c;
    if (data_d.count("c") > 0){
        c = data_d.at("c");
    }    
    else{
        c = params[0];
    }
    double p_c = params[params.size()-1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    double binom_p = (1.0 - c) * p_e + c * p_c;
    
    if (isnan(logbinom(n,k,binom_p)) || isinf(logbinom(n,k,binom_p))){
       fprintf(stderr, "oops\n");
       exit(1);

    }
    return logbinom(n, k, binom_p);
}

/**
 * Derivative of log likelihood function for updating
 * ambient RNA profile mixture proportions
 * (optimML-format)
 */
void dll_amb_prof_mixture(const vector<double>& params, 
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double c;
    bool c_in_data = false;
    if (data_d.count("c") > 0){
        c = data_d.at("c");
        c_in_data = true;
    }
    else{
        c = params[0];
    }
    double p_c = params[params.size()-1];
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    
    double binom_p = (1.0 - c) * p_e + c * p_c;
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    
    results[results.size()-1] += dy_dp * c;
    if (c_in_data){
        results[0] += dy_dp * (p_c - p_e);
    }
}

// ===== contamFinder functions =====

/**
 * Class constructor
 *
 * Stores copies of data structures/values that will be used by
 * many functions.
 *
 * Compiles data in the form that will be needed by later functions.
 */
contamFinder::contamFinder(robin_hood::unordered_map<unsigned long, 
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    map<pair<int, int>, map<int, float> >& exp_match_fracs,
    int n_samples,
    set<int>& allowed_ids,
    set<int>& allowed_ids2){
    
    this->contam_prof_initialized = false;
    this->c_init = -1;

    // Copy external data structures that we will need in the future 
    this->assn = assn;
    this->assn_llr = assn_llr;
    this->indv_allelecounts = indv_allelecounts;
    this->n_samples = n_samples;
    this->n_mixprop_trials = 10;
    this->expfracs = exp_match_fracs;
    this->allowed_ids = allowed_ids;
    this->allowed_ids2 = allowed_ids2;
    
    // Reassign cells by default
    this->skip_reassign = false;

    // Don't worry about doublet rate by default
    this->doublet_rate = -1;
    
    this->ef_all_avg = true;

    llrtot = 0.0;

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        allowed_ids.insert(a->second);
        // Make sure we also allow sub-IDs of combinations
        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            allowed_ids.insert(combo.first);
            allowed_ids.insert(combo.second);
        }
        if (id_llrsum.count(a->second) == 0){
            id_llrsum.insert(make_pair(a->second, 0.0));
            id_count.insert(make_pair(a->second, 0.0));
        }
        id_llrsum[a->second] += assn_llr[a->first];
        id_count[a->second]++;
        llrtot += assn_llr[a->first];
        if (a->second >= n_samples){
            pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
            if (id_llrsum2.count(comb.first) == 0){
                id_llrsum2.insert(make_pair(comb.first, 0.0));
            }
            if (id_llrsum2.count(comb.second) == 0){
                id_llrsum2.insert(make_pair(comb.second, 0.0));
            }
            id_llrsum2[comb.first] += 0.5*assn_llr[a->first];
            id_llrsum2[comb.second] += 0.5*assn_llr[a->first];
        }
        else{
            if (id_llrsum2.count(a->second) == 0){
                id_llrsum2.insert(make_pair(a->second, 0.0));
            }
            id_llrsum2[a->second] += assn_llr[a->first];
        }
    }

    // For the ambient RNA mixture, we only want to model single (not doublet) 
    // individual genotypes that were allowed in this run (if no filter was given,
    // allow all single individuals from the VCF).

    if (allowed_ids.size() == 0){
        // Default to all possible individuals
        for (int i = 0; i < n_samples; ++i){
            idx2samp.push_back(i);
        }
    }
    else{
        for (set<int>::iterator a = allowed_ids.begin(); a != allowed_ids.end(); ++a){
            if (*a < n_samples){
                idx2samp.push_back(*a);
            }
        }
    }

    // Set some default parameter values.
    this->e_r = 1e-3;
    this->e_a = 1e-3;
    
    this->inter_species = false;

    this->contam_cell_prior = -1;
    this->contam_cell_prior_var = -1;

    this->delta_thresh = 0.1;
    this->maxits = 100;
    
    this->weighted = false;

    // Compile data in the format needed by other functions
    this->compile_data(assn, indv_allelecounts);
}

void contamFinder::set_init_contam_prof(map<int, double>& cp){
    this->contam_prof = cp;
    contam_prof_initialized = true;
}

/**
 * Functions to set parameter values
 */

void contamFinder::set_error_rates(double error_ref, double error_alt){
    this->e_r = error_ref;
    this->e_a = error_alt;
}

void contamFinder::set_doublet_rate(double d){
    this->doublet_rate = d;
}

void contamFinder::model_other_species(){
    this->inter_species = true;
}

void contamFinder::model_single_species(){
    this->inter_species = false;
}

void contamFinder::set_mixprop_trials(int nt){
    this->n_mixprop_trials = nt;
}

void contamFinder::set_delta(double delta){
    this->delta_thresh = delta;
}

void contamFinder::set_maxiter(int its){
    this->maxits = its;
}

void contamFinder::use_weights(){
    this->weighted = true;
}

void contamFinder::no_weights(){
    this->weighted = false;
}

void contamFinder::no_reassign(){
    this->skip_reassign = true;
}

void contamFinder::set_init_c(double c){
    c_init = c;
}

void contamFinder::set_num_threads(int nt){
    num_threads = nt;
}

/**
 * Populates internal data structures with data and expected values.
 */
void contamFinder::compile_data(robin_hood::unordered_map<unsigned long, int>& assn,
     robin_hood::unordered_map<unsigned long, 
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts){

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        vector<double> n;
        vector<double> k;
        vector<double> p_e;
        vector<pair<int, int> > type1;
        vector<pair<int, int> > type2; 

        this->get_reads_expectations(a->second, indv_allelecounts[a->first],
            n, k, p_e, type1, type2);
        
        for (int i = 0; i < n.size(); ++i){
            int idx_all = n_all.size();

            n_all.push_back(n[i]);
            k_all.push_back(k[i]);
            p_e_all.push_back(p_e[i]);
            type1_all.push_back(type1[i]);
            type2_all.push_back(type2[i]);

            if (expfrac_to_idx.count(type1[i]) == 0){
                map<pair<int, int>, vector<int> > m;
                expfrac_to_idx.insert(make_pair(type1[i], m));
            }
            if (expfrac_to_idx[type1[i]].count(type2[i]) == 0){
                vector<int> v;
                expfrac_to_idx[type1[i]].insert(make_pair(type2[i], v));
            }
            expfrac_to_idx[type1[i]][type2[i]].push_back(idx_all);
            
            idx_to_cell.insert(make_pair(idx_all, a->first));
            
            if (cell_to_idx.count(a->first) == 0){
                vector<int> v;
                cell_to_idx.insert(make_pair(a->first, v));
            }
            cell_to_idx[a->first].push_back(idx_all);
            
            pair<int, int> sitekey;
            if (type2[i].second == -1){
                sitekey = make_pair(type1[i].second, -1);
            }
            else{
                int thistype1 = type1[i].second;
                int thistype2 = type2[i].second;
                if (thistype2 < thistype1){
                    int tmp = thistype2;
                    thistype2 = thistype1;
                    thistype1 = tmp;
                }
                sitekey = make_pair(thistype1, thistype2);
            }
            if (sitecomb_type_to_idx.count(sitekey) == 0){
                vector<int> v;
                sitecomb_type_to_idx.insert(make_pair(sitekey, v));
            }
            sitecomb_type_to_idx[sitekey].push_back(idx_all);
        }
    }
}

void contamFinder::clear_data(){
    n_all.clear();
    k_all.clear();
    p_e_all.clear();
    type1_all.clear();
    type2_all.clear();
    expfrac_to_idx.clear();
    idx_to_cell.clear();
    cell_to_idx.clear();
    sitecomb_type_to_idx.clear();
}

/**
 * Helper function for compile_data.
 *
 * Given an assigned identity for a cell, collects all reads that should contain
 * each different type of allelic state and organize data in a way that is amenable
 * to solving for the parameters.
 */
void contamFinder::get_reads_expectations(int ident,
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > >& allelecounts,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<pair<int, int> >& type1,
    vector<pair<int, int> >& type2){
    
    static pair<int, int> nullkey = make_pair(-1, -1);

    bool is_combo = false;
    pair<int, int> combo;
    if (ident >= n_samples){
        combo = idx_to_hap_comb(ident, n_samples);
        is_combo = true;
    }
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac = 
        allelecounts.begin(); ac != allelecounts.end(); ++ac){
        if (is_combo){
            if (ac->first.first == combo.first){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = ac->second.begin();
                    ac2 != ac->second.end(); ++ac2){
                    if (ac2->first.first == combo.second){
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        if (ref + alt > 0){

                            // Store expectations without error incorporated (yet)
                            double expected = (double)(ac->first.second + ac2->first.second) / 4.0;
                            
                            n.push_back(ref+alt);
                            k.push_back(alt);
                            p_e.push_back(expected);
                            type1.push_back(ac->first);
                            type2.push_back(ac2->first);
                        }
                    }
                    else if (ac2->first.first > combo.second){
                        break;
                    }
                    else{
                        continue;
                    }
                }
            }
            else if (ac->first.first > combo.first){
                break;
            }
            else{
                continue;
            }
        }
        else{
            if (ac->first.first == ident){
                 // Store expectations without error incorporated (for now)
                double expected = (double)ac->first.second / 2.0;
                double ref = ac->second[nullkey].first;
                double alt = ac->second[nullkey].second;
                if (ref+alt > 0){
                    n.push_back(ref+alt);
                    k.push_back(alt);
                    p_e.push_back(expected);
                    type1.push_back(ac->first);
                    type2.push_back(nullkey);
                }
            }
            else if (ac->first.first > ident){
                break;
            }
            else{
                continue;
            }
        }
    }
}

/**
 * When starting out, we don't know what the initial guess of contamination rate
 * should be. This is a hard problem, since we also don't know what the contamination
 * profile looks like yet.
 *
 * For each type of SNP, we can estimate the minimum bound on contamination rate
 *   in a given cell by assuming that ambient RNA has an alternate allele frequency
 *   as far as possible in the direction needed to "fix" the measured frequency.
 *
 *   In other words, if a cell belongs to individual 1 and we measure an alt allele
 *     frequency of 0.9 in sites where individual 1 is heterozygous alt, a minimum bound
 *     on the contamination rate could be found by assuming ambient RNA has 100% reference
 *     alleles at these sites (and c = 0.1). 
 *
 *   If a cell has alt allele frequency 0.67 at sites where the individual is heterozygous,
 *     then we can imagine ambient RNA is 100% alt alleles at these sites.
 *
 *     We then compute c = (measured-expected)/(ambient-expected), giving c = 0.34.
 *
 *     This gives us many measurements of c, and we return a weighted average.
 */
double contamFinder::est_min_c(){
    map<pair<int, int>, map<pair<int, int>, double> > mean_f;
    map<pair<int, int>, map<pair<int, int>, double> > counts;
    
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            for (int nalt1 = 0; nalt1 <= 2; ++nalt1){
                pair<int, int> key1 = make_pair(combo.first, nalt1);
                for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                    pair<int, int> key2 = make_pair(combo.second, nalt2);
                    double ref = indv_allelecounts[a->first][key1][key2].first;
                    double alt = indv_allelecounts[a->first][key1][key2].second;
                    if (ref + alt == 0){
                        continue;
                    }
                    double f = alt/(ref+alt);
                    double p0 = (double)(nalt1 + nalt2)/4.0;
                    if (mean_f.count(key1) == 0){
                        map<pair<int, int>, double> m;
                        mean_f.insert(make_pair(key1, m));
                        counts.insert(make_pair(key1, m));
                    }
                    if (mean_f[key1].count(key2) == 0){
                        mean_f[key1].insert(make_pair(key2, 0));
                        counts[key1].insert(make_pair(key2, 0));
                    }
                    if (!weighted){
                        mean_f[key1][key2] += f;    
                        counts[key1][key2]++;
                    }
                    else{
                        mean_f[key1][key2] += f * assn_llr[a->first];
                        counts[key1][key2] += assn_llr[a->first];
                    }
                }
            }
        }
        else{
            pair<int, int> nullkey = make_pair(-1, -1);
            for (int nalt = 0; nalt <= 2; ++nalt){
                pair<int, int> key = make_pair(a->second, nalt);
                double ref = indv_allelecounts[a->first][key][nullkey].first;
                double alt = indv_allelecounts[a->first][key][nullkey].second;
                if (ref + alt == 0){
                    continue;
                }
                double f = alt/(ref+alt);
                double p0 = (double)nalt/2.0;
                if (mean_f.count(key) == 0){
                    map<pair<int, int>, double> m;
                    mean_f.insert(make_pair(key, m));
                    counts.insert(make_pair(key, m));
                }
                if (mean_f[key].count(nullkey) == 0){
                    mean_f[key].insert(make_pair(nullkey, 0));
                    counts[key].insert(make_pair(nullkey, 0));
                }
                if (!weighted){
                    mean_f[key][nullkey] += f;
                    counts[key][nullkey]++;
                }
                else{
                    mean_f[key][nullkey] += f*assn_llr[a->first];
                    counts[key][nullkey] += assn_llr[a->first];
                }
            }
        }
    }

    double minc_mean = 0.0;
    double minc_weightsum = 0.0;
    
    double minc_min = -1;
    double minc_max = -1;
    
    map<int, double> minc_by_id;
    map<int, double> minc_by_id_count;

    for (map<pair<int, int>, map<pair<int, int>, double> >::iterator d = mean_f.begin();
        d != mean_f.end(); ++d){
        for (map<pair<int, int>, double>::iterator d2 = d->second.begin(); d2 != d->second.end();
            ++d2){
            double expec;
            if (d2->first.first == -1){
                expec = d->first.second / 2.0;
            }
            else{
                expec = (d->first.second + d2->first.second)/4.0;
            }
            expec = adjust_p_err(expec, e_r, e_a);
            if (counts[d->first][d2->first] > 0){
               
                double inf = d2->second / counts[d->first][d2->first];
                double minc;
                if (inf > expec){
                    // treat p_c = 1
                    minc = (inf - expec)/(1.0 - expec);
                }   
                else{
                    // treat p_c = 0
                    minc = (inf - expec)/(0.0 - expec);
                }
                if (minc_min == -1 || minc < minc_min){
                    minc_min = minc;
                }
                if (minc_max == -1 || minc > minc_max){
                    minc_max = minc;
                }
                minc_mean += minc * counts[d->first][d2->first];
                minc_weightsum += counts[d->first][d2->first];
                
                if (minc_by_id.count(d->first.first) == 0){
                    minc_by_id.insert(make_pair(d->first.first, 0));
                    minc_by_id_count.insert(make_pair(d->first.first, 0));
                }
                
                if (d2->first.first != -1){
                    double w = 0.5*(double)counts[d->first][d2->first];
                    minc_by_id[d->first.first] += w*minc;
                    minc_by_id_count[d->first.first] += w;
                    if (minc_by_id.count(d2->first.first) == 0){
                        minc_by_id.insert(make_pair(d2->first.first, 0));
                        minc_by_id_count.insert(make_pair(d2->first.first, 0));
                    }
                    minc_by_id[d2->first.first] += w*minc;
                    minc_by_id_count[d2->first.first] += w;
                }
                else{
                    double w = (double)counts[d->first][d2->first];
                    minc_by_id[d->first.first] += w*minc;
                    minc_by_id_count[d->first.first] += w;
                }
            }
        }
    }
    
    double vsum = 0.0;
    double vcount = 0.0;
    for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
        vsum += mi->second/minc_by_id_count[mi->first];
        vcount++;
    }
    double c_est = vsum/(vcount-1.0);
    
    contam_prof.clear();
    double minval = 0.01;
    double denom = 0.0;
    for (map<int, double>::iterator mi = minc_by_id.begin(); mi != minc_by_id.end(); ++mi){
        double val = mi->second/minc_by_id_count[mi->first];
        double frac = 1.0 - val/c_est;
        if (frac < minval){
            frac = minval;
        }
        denom += frac;
        contam_prof.insert(make_pair(mi->first, frac));
    }
    if (inter_species){
        contam_prof.insert(make_pair(-1, 1.0/((double)minc_by_id.size() + 1.0)));
        denom += contam_prof[-1];
    }
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        cp->second /= denom;
    }
    
    return c_est;

    /*
    minc_mean /= minc_weightsum;
    fprintf(stderr, "mean %f\n", minc_mean);
    exit(0);

    return minc_mean;
    */
}

/**
 * When running the first time, create a starting guess of contamination profile
 * and infer the likeliest mixture of individuals.
 */
double contamFinder::init_params(double& init_c){
    
    if (!contam_prof_initialized){
        // Init contam prof
        //contam_prof.clear();
        
        // This will compute an initial contamination profile based on the counts of different
        // types of cells.

        // The commented-out block will compute an initial contamination profile where
        // all proportions are equal.
        /*
        map<int, int> acounts;
        int atots = 0;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (a->second >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                if (acounts.count(comb.first) == 0){
                    acounts.insert(make_pair(comb.first, 0));
                }
                if (acounts.count(comb.second) == 0){
                    acounts.insert(make_pair(comb.second, 0));
                }
                acounts[comb.first]++;
                acounts[comb.second]++;
                atots += 2;
            }
            else{
                if (acounts.count(a->second) == 0){
                    acounts.insert(make_pair(a->second, 0));
                }
                acounts[a->second] += 2;
                atots += 2;
            }
        }
        double minfrac = 0.001;
        double proptot = 0.0;
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            double prop = (double)acounts[samp]/(double)atots;
            if (prop < minfrac){
                prop = minfrac;
            }
            contam_prof.insert(make_pair(samp, prop));
            proptot += prop;
        }
        if (inter_species){
            double prop = 1.0/(double)(idx2samp.size() + 1);
            contam_prof.insert(make_pair(-1, prop));
            proptot += prop;
        }
        for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end();
            ++cp){
            cp->second /= proptot;
        }
        */
        
        /*
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            if (inter_species){
                contam_prof.insert(make_pair(samp, 1.0/(double)(idx2samp.size() + 1)));
            }
            else{
                contam_prof.insert(make_pair(samp, 1.0/(double)idx2samp.size()));
            }
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, 1.0/(double)(idx2samp.size() + 1)));
        }
        */

    }
    return update_amb_prof_mixture(true, init_c, true);
}

/**
 * Produces a single, global contamination rate estimate.
 * Stores this as contam_cell_prior.
 */
void contamFinder::est_contam_cells_global(){
    
    contam_rate.clear();
    contam_rate_se.clear();
    
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    
    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin(); 
        ci != cell_to_idx.end(); ++ci){
        
        for (vector<int>::iterator i = ci->second.begin(); i != ci->second.end(); ++i){
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
            p_c.push_back(amb_mu[type1_all[*i]][type2_all[*i]]);
        }
    }
    optimML::brent_solver c_global(ll_c, dll_dc, d2ll_dc2);
    if (num_threads > 1){
        c_global.set_threads(num_threads);
    }
    c_global.add_data("n", n);
    c_global.add_data("k", k);
    c_global.add_data("p_e", p_e);
    c_global.add_data("p_c", p_c);
    c_global.constrain_01();
    c_global.set_maxiter(-1);
    contam_cell_prior = c_global.solve(0,1);
    contam_cell_prior_var = -1;
    fprintf(stderr, "MLE global contamination rate: %f\n", contam_cell_prior);
}

/**
 * Once an ambient RNA profile exists (we have estimates of p_c parameters for
 * every category of SNP, we can re-estimate the likeliest contamination rate per
 * cell given these parameters.
 *
 * Uses an Emprical Bayes method -- a prior distribution on contamination rate
 * per cell is based on the data set-wide distribution from the last round of 
 * estimate.
 *
 */
void contamFinder::est_contam_cells(){
    
    contam_rate.clear();
    contam_rate_se.clear();
    contam_rate_ll.clear();

    // Store all successful estimates of contamination rate for computing the
    // mean and variance across the data set 
    vector<double> cell_c_maps;
    // Store cell weights (log likelihood ratio of most likely ID assignment)
    // to use in calculating this mean & variance
    vector<double> cell_c_llr;

    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin(); 
        ci != cell_to_idx.end(); ++ci){
        
        // Compile data for this cell    
        vector<double> n;
        vector<double> k;
        vector<double> p_e;
        vector<double> p_c;

        for (vector<int>::iterator i = ci->second.begin(); i != ci->second.end(); ++i){
            n.push_back(n_all[*i]);
            k.push_back(k_all[*i]);
            p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
            p_c.push_back(amb_mu[type1_all[*i]][type2_all[*i]]);
        }
        optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
        if (num_threads > 1){
            //c_cell.set_threads(num_threads);
        }

        c_cell.add_data("n", n);
        c_cell.add_data("k", k);
        c_cell.add_data("p_e", p_e);
        c_cell.add_data("p_c", p_c);
        c_cell.constrain_01();
        if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
            // We already ran once to get all per-cell maximum likelihood estimates.
            // Now fit a Beta distribution to these values and use it as a prior to shrink
            // per-cell estimates (Empirical Bayes)
            pair<double, double> bm = beta_moments(contam_cell_prior, contam_cell_prior_var);
            c_cell.add_beta_prior(bm.first, bm.second);
        }
        c_cell.set_maxiter(-1);
        
        bool root_found = 0.0;
        double se = 0.0;
        double c_cell_map = 1.0;
        double ll = 0.0;
        try{
            c_cell_map = c_cell.solve(0,1);
            if (c_cell.root_found){
                ll = c_cell.log_likelihood;
                if (c_cell.se_found){
                    se = c_cell.se;
                }
            }
        }
        catch (int exc){
            // pass
        }
        // Set SE to 0 if we did not find a maximum-likelihood estimate in the range (0,1) 
        contam_rate.emplace(ci->first, c_cell_map);
        cell_c_maps.push_back(c_cell_map);
        if (weighted){
            double weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
            cell_c_llr.push_back(weight);
        }
        contam_rate_se.emplace(ci->first, se);
        contam_rate_ll.emplace(ci->first, ll);
    }
    
    // Re-compute data set-wide distribution
    pair<double, double> mu_var;
    if (weighted){
        mu_var = welford_weights(cell_c_maps, cell_c_llr, false);
    }
    else{
        mu_var = welford(cell_c_maps);
    }

    if (weighted && mu_var.second < 1e-3){
        // If the variance is too low, try to re-compute without using weights
        mu_var = welford(cell_c_maps);
    }
    if (mu_var.second > 1e-6){
        if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
            fprintf(stderr, "Shrunken per-cell contamination rates:\n");
        }
        else{
            fprintf(stderr, "Per-cell contamination rates:\n");
        }
        fprintf(stderr, "  Mean: %f Std dev: %f\n", mu_var.first, sqrt(mu_var.second));
        contam_cell_prior = mu_var.first;
        contam_cell_prior_var = mu_var.second;
    }

    map<int, double> idcsum;
    map<int, double> idccount;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (a->second >= n_samples){
            pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
            if (idcsum.count(comb.first) == 0){
                idcsum.insert(make_pair(comb.first, 0));
                idccount.insert(make_pair(comb.first,0));
            }
            if (idcsum.count(comb.second) == 0){
                idcsum.insert(make_pair(comb.second, 0));
                idccount.insert(make_pair(comb.second,0));
            }
            idcsum[comb.first] += 0.5*contam_rate[a->first];
            idcsum[comb.second] += 0.5*contam_rate[a->first];
            idccount[comb.first] += 0.5;
            idccount[comb.second] += 0.5;
        }
        else{
            if (idcsum.count(a->second) == 0){
                idcsum.insert(make_pair(a->second, 0));
                idccount.insert(make_pair(a->second, 0));
            }
            idcsum[a->second] += contam_rate[a->first];
            idccount[a->second] += 1.0;
        }
    }
    for (map<int, double>::iterator c = idcsum.begin(); c != idcsum.end(); ++c){
        //fprintf(stderr, "contam mean %d) %f\n", c->first, c->second/idccount[c->first]);
    }
}

/**
 * Updates ambient RNA profile based on current cell-contamination estimates,
 * without consideration of the fraction of individuals making up ambient
 * RNA.
 *
 * Currently not used.
 */
double contamFinder::update_ambient_profile(bool global_c){
    
    vector<pair<int, int> > idx2expfrac1;
    vector<pair<int, int> > idx2expfrac2;
    
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    vector<int> efp_idx;
    vector<double> c;
    vector<double> weights;

    map<int, int> idx2idx;
    vector<double> efparams;
    map<int, pair<int, int> > efp2ef1;
    map<int, pair<int, int> > efp2ef2;

    double floor = 1e-3;
    double ceil = 1-1e-3;

    for (map<pair<int, int>, map<pair<int, int>, vector<int> > >::iterator ei1 = 
        expfrac_to_idx.begin(); ei1 != expfrac_to_idx.end(); ++ei1){
        for (map<pair<int, int>, vector<int> >::iterator ei2 = ei1->second.begin(); ei2 != 
            ei1->second.end(); ++ei2){
            
            int ef_param_idx = efparams.size();
            double p_c = amb_mu[ei1->first][ei2->first];
            if (p_c < floor){
                p_c = floor;
            }
            else if (p_c > ceil){
                p_c = ceil;
            }
            efparams.push_back(p_c);
            efp2ef1.insert(make_pair(ef_param_idx, ei1->first));
            efp2ef2.insert(make_pair(ef_param_idx, ei2->first));
            
            for (vector<int>::iterator i = ei2->second.begin(); i != ei2->second.end(); ++i){
                if (global_c || (idx_to_cell.count(*i) > 0 && 
                    contam_rate.count(idx_to_cell[*i]) > 0)){ 
                    
                    idx2idx.insert(make_pair(*i, n.size()));    
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                    efp_idx.push_back(ef_param_idx);
                    if (global_c){
                        c.push_back(contam_cell_prior);
                    }
                    else{
                        c.push_back(contam_rate[idx_to_cell[*i]]);
                    }
                    double weight = 1.0;
                    if (weighted){
                        double llr = assn_llr[idx_to_cell[*i]];
                        double llrtot = id_llrsum[assn[idx_to_cell[*i]]];
                        weight = llr / llrtot;
                    }
                    weights.push_back(weight);
                }
            }        
        }
    }
        
    optimML::multivar_ml_solver solver(efparams, ll_ambmu, dll_ambmu);
    if (num_threads > 1){
        solver.set_threads(num_threads);
        // Many parameters - use threads
        solver.set_bfgs_threads(num_threads);
    }
    for (int i = 0; i < efparams.size(); ++i){
        solver.constrain_01(i);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_data("ef_idx", efp_idx);
    solver.add_data("c", c);
    solver.add_weights(weights);
    solver.solve();
    
    for (int i = 0; i < efparams.size(); ++i){
        double updated = solver.results[i];
        if (amb_mu.count(efp2ef1[i]) == 0){
            map<pair<int, int>, double> m;
            amb_mu.insert(make_pair(efp2ef1[i], m));
        }
        if (amb_mu[efp2ef1[i]].count(efp2ef2[i]) == 0){
            amb_mu[efp2ef1[i]].insert(make_pair(efp2ef2[i], 0.0));
        }
        amb_mu[efp2ef1[i]][efp2ef2[i]] = updated;
    }
    return solver.log_likelihood;
}

/**
 * Compiles data for update_amb_prof_mixture(). Puts necessary values
 * in vectors, which can then be bootstrapped.
 */
void contamFinder::compile_amb_prof_dat(bool solve_for_c, 
    bool use_global_c,
    vector<vector<double> >& mixfracs,
    vector<double>& weights,
    vector<double>& n,
    vector<double>& k,
    vector<double>& p_e,
    vector<double>& c){
    
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        
        bool is_comb = false;
        pair<int, int> comb;
        if (a->second >= n_samples){
            is_comb = true;
            comb = idx_to_hap_comb(a->second, n_samples);
        }
        
        // get weight
        double weight = 1.0;
        if (weighted){
            weight = assn_llr[a->first] / id_llrsum[a->second];
        }

        double this_c;
        if (!solve_for_c){
            if (use_global_c){
                this_c = contam_cell_prior;
            }
            else{
                this_c = contam_rate[a->first];
            }
        }
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac1 = 
            indv_allelecounts[a->first].begin(); ac1 != indv_allelecounts[a->first].end(); 
            ++ac1){
            
            if ((!is_comb && ac1->first.first == a->second) || 
                (is_comb && ac1->first.first == comb.first)){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = 
                    ac1->second.begin(); ac2 != ac1->second.end(); ++ac2){
                    
                    if ((!is_comb && ac2->first.first == -1) || 
                        (is_comb && ac2->first.first == comb.second)){
                        
                        double expected;
                        if (!is_comb){
                            expected = adjust_p_err((double)ac1->first.second / 2.0, 
                                e_r, e_a);
                        }
                        else{
                            expected = adjust_p_err((double)(ac1->first.second + 
                                ac2->first.second)/4.0, e_r, e_a);
                        }
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        
                        n.push_back(ref+alt);
                        k.push_back(alt);
                        weights.push_back(weight);
                        if (!solve_for_c){
                            c.push_back(this_c);
                        }
                        p_e.push_back(expected);
                        vector<double> mixfrac_row;
                        for (int i = 0; i < idx2samp.size(); ++i){
                            int samp = idx2samp[i];
                            if (ac2->first.first == -1){
                                mixfrac_row.push_back(expfracs[ac1->first][samp]);
                            }
                            else{
                                if (ef_all_avg && ac1->first.first == samp){
                                    mixfrac_row.push_back(adjust_p_err(
                                        ac1->first.second / 2.0, e_r, e_a));
                                }
                                else if (ef_all_avg && ac2->first.first == samp){
                                    mixfrac_row.push_back(adjust_p_err(
                                        ac2->first.second / 2.0, e_r, e_a));
                                }
                                else{
                                    mixfrac_row.push_back(0.5 * expfracs[ac1->first][samp] + 
                                        0.5 * expfracs[ac2->first][samp]);
                                }
                            }
                        }
                        if (inter_species){
                            // Reference alleles
                            mixfrac_row.push_back(adjust_p_err(0.0, e_r, e_a));
                        }
                        mixfracs.push_back(mixfrac_row);
                    }
                }
            }
        }
    }
}

/**
 * Models ambient RNA as a mixture of individuals. Updates the 
 * ambient RNA alt allele fractions so est_contam_cells() can find
 * contam levels in cells.
 *
 * Optionally can solve for global contam rate, and can choose whether to
 * use global contam rate (stored) or per-cell estimates in calculations.
 *
 * Returns log likelihood.
 * Edits init_c to updated value (if solve_for_c == true)
 */
double contamFinder::update_amb_prof_mixture(bool solve_for_c, double& init_c, bool use_global_c){
    
    vector<double> params;
    if (solve_for_c){
        params.push_back(init_c);
    }
    optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture, dll_amb_prof_mixture);
    if (num_threads > 1){
        // Avoid multi-threading for evaluation for mixture proportion problems
        //solver.set_threads(num_threads);
        solver.set_bfgs_threads(num_threads);
    }

    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    // Get data
    compile_amb_prof_dat(solve_for_c, use_global_c, mixfracs, weights,
        n, k, p_e, c); 
    
    // Starting proportions should be those previously set
    vector<double> startprops;
    for (map<int, double>::iterator cp = contam_prof.begin(); 
        cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            startprops.push_back(cp->second);
        }
    }
    if (contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
    }
    
    // Set up ML solver 
    solver.add_mixcomp(mixfracs);
    solver.add_mixcomp_fracs(startprops);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_weights(weights);
    if (solve_for_c){
        solver.constrain_01(0);
    }
    else{
        solver.add_data("c", c);
    }
    //solver.set_delta(1e-6);

    if (solve_for_c){
        // First time. Try a few starting conditions and get the maximum LL.
        solver.solve();
        
        vector<double> lls;
        vector<vector<double> > mcs;
        vector<double> cs;
        int maxidx = 0;
        double maxll = solver.log_likelihood;
        lls.push_back(solver.log_likelihood);
        mcs.push_back(solver.results_mixcomp);
        cs.push_back(solver.results[0]);
        
        vector<double> mptest;
        for (int i = 0; i < startprops.size(); ++i){
            mptest.push_back(1.0/startprops.size());
        }
        solver.add_mixcomp_fracs(mptest);
        solver.solve();
        if (solver.log_likelihood > maxll){
            maxll = solver.log_likelihood;
        }
        lls.push_back(solver.log_likelihood);
        mcs.push_back(solver.results_mixcomp);
        cs.push_back(solver.results[0]);

        int ntrials = contam_prof.size() * n_mixprop_trials;
        //vector<double> trialprops;
        for (int n = 0; n < ntrials; ++n){
            solver.randomize_mixcomps();
             
            solver.solve();
            if (solver.log_likelihood > maxll){
                maxll = solver.log_likelihood;
                maxidx = lls.size();
                lls.push_back(solver.log_likelihood);
                mcs.push_back(solver.results_mixcomp);
                cs.push_back(solver.results[0]);
            }
        }
        solver.log_likelihood = maxll;
        solver.results[0] = cs[maxidx];
        solver.results_mixcomp = mcs[maxidx];
        // Update stored contamination profile (mixture of individuals)
        contam_prof.clear();
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            contam_prof.insert(make_pair(samp, solver.results_mixcomp[i]));
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, 
                solver.results_mixcomp[solver.results_mixcomp.size()-1]));
        }
    }
    else{
        vector<double> trialprops;
        solver.solve();
        vector<double> maxres = solver.results_mixcomp;
        double maxll = solver.log_likelihood;
        /* 
        for (int i = 0; i < contam_prof.size() * n_mixprop_trials; ++i){
            //rdirichlet(startprops, trialprops);
            //solver.add_mixcomp_fracs(trialprops);
            solver.randomize_mixcomps();
            solver.solve();
            if (maxll == 0 || solver.log_likelihood > maxll){
                maxres = solver.results_mixcomp;
                maxll = solver.log_likelihood;
            }
        }
        */
        // Update stored contamination profile (mixture of individuals)
        contam_prof.clear();
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            contam_prof.insert(make_pair(samp, maxres[i]));
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, 
                maxres[maxres.size()-1]));
        }
    }
    
    // Update stored ambient RNA profile (allele matching fractions)
    for (int x = 0; x < idx2samp.size(); ++x){
        int i = idx2samp[x];
        for (int nalt = 0; nalt <= 2; ++nalt){
            pair<int, int> key = make_pair(i, nalt);
            if (amb_mu.count(key) == 0){
                map<pair<int, int>, double> m;
                amb_mu.insert(make_pair(key, m));
            }
            pair<int, int> nullkey = make_pair(-1,-1);
            if (amb_mu[key].count(nullkey) == 0){
                amb_mu[key].insert(make_pair(nullkey, 0));
            }
            double val = 0.0;
            for (map<int, double>::iterator cp = contam_prof.begin(); 
                cp != contam_prof.end(); ++cp){
                
                if (cp->first == -1){
                    val += cp->second * 0.0;
                }
                else if (ef_all_avg && cp->first == key.first){
                    val += cp->second * ((double)key.second/2.0);
                }
                else{
                    val += cp->second * expfracs[key][cp->first];
                }
            }
            amb_mu[key][nullkey] = val;
            for (int y = x + 1; y < idx2samp.size(); ++y){
                int j = idx2samp[y];
                for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                    pair<int, int> key2 = make_pair(j, nalt2);
                    if (amb_mu[key].count(key2) == 0){
                        amb_mu[key].insert(make_pair(key2, 0));
                    }
                    double val = 0.0;
                    for (map<int, double>::iterator cp = contam_prof.begin(); 
                        cp != contam_prof.end(); ++cp){
                        
                        if (cp->first == -1){
                            val += cp->second * adjust_p_err(0.0, e_r, e_a);
                        }
                        else{
                            if (ef_all_avg && cp->first == key.first){
                                val += cp->second * adjust_p_err(
                                    (double)key.second / 2.0, e_r, e_a);
                            }
                            else if (ef_all_avg && cp->first == key2.first){
                                val += cp->second * adjust_p_err(
                                    (double)key2.second / 2.0, e_r, e_a);
                            }
                            else{
                                val += cp->second *
                                    (0.5*expfracs[key][cp->first] + 0.5*expfracs[key2][cp->first]);
                            }
                        }
                    }
                    amb_mu[key][key2] = val;
                }
            }
        }
    }
    if (solve_for_c){
        init_c = solver.results[0];
    }
    return solver.log_likelihood;
}

/**
 * Gets variances on contamination profile proportions by bootstrapping.
 * Fits a Dirichlet distribution to all bootstrap samples; Dirichlet
 * concentration parameters will then be reported in output files.
 */
void contamFinder::bootstrap_amb_prof(int n_boots, map<int, double>& dirichlet_params){
    // Assumes we have already solved everything and that this is at the end.

    // Compile everything we need
    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    compile_amb_prof_dat(false, false, mixfracs, weights, n, k, p_e, c);
   
    // Store MLEs from bootstrap samples, which will serve as samples from
    // a Dirichlet distribution 
    vector<vector<double> > dirprops;
    
    // Store MLE solution from contam_prof
    vector<double> mle_fracs;

    // Starting proportions should be those previously set
    vector<double> startprops;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            startprops.push_back(cp->second);
            vector<double> v;
            dirprops.push_back(v);
            mle_fracs.push_back(cp->second);
        }
    }
    if (contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
        vector<double> v;
        dirprops.push_back(v);
        mle_fracs.push_back(contam_prof[-1]);
    }
    
    // Initialize random stuff
    static random_device dev;
    static mt19937 rand_gen = mt19937(dev());
    static uniform_int_distribution<int> uni_dist(0, n.size()-1);
    
    for (int b = 0; b < n_boots; ++b){
        fprintf(stderr, "Bootstrap sample %d...\r", b+1);
        
        // Re-sample
        vector<vector<double> > mixfracs_boot;
        vector<double> weights_boot;
        vector<double> n_boot;
        vector<double> k_boot;
        vector<double> p_e_boot;
        vector<double> c_boot;
        
        for (int x = 0; x < n.size(); ++x){
            int r = uni_dist(rand_gen);
            mixfracs_boot.push_back(mixfracs[r]);
            weights_boot.push_back(weights[r]);
            n_boot.push_back(n[r]);
            k_boot.push_back(k[r]);
            p_e_boot.push_back(p_e[r]);
            c_boot.push_back(c[r]);
        }    
        
        // Solve
        vector<double> params;
        optimML::multivar_ml_solver solver(params, ll_amb_prof_mixture, dll_amb_prof_mixture);
        if (num_threads > 1){
            // Avoid multithreading for evaluation for mixture proportion problems
            //solver.set_threads(num_threads);
            solver.set_bfgs_threads(num_threads);
        }
        solver.add_mixcomp(mixfracs_boot);
        solver.add_mixcomp_fracs(startprops);
        solver.add_data("n", n_boot);
        solver.add_data("k", k_boot);
        solver.add_data("p_e", p_e_boot);
        solver.add_data("c", c_boot);
        solver.solve();
        
        for (int x = 0; x < solver.results_mixcomp.size(); ++x){
            dirprops[x].push_back(solver.results_mixcomp[x]);
        }
        
    }
    fprintf(stderr, "\n");
    
    // Now fit Dirichlet MLE
    vector<double> dirichlet_soln;
    fit_dirichlet(mle_fracs, dirprops, dirichlet_soln);
    
    dirichlet_params.clear();
    for (int i = 0; i < idx2samp.size(); ++i){
        int samp = idx2samp[i];
        dirichlet_params.insert(make_pair(samp, dirichlet_soln[i]));
    }
    if (inter_species){
        dirichlet_params.insert(make_pair(-1, 
           dirichlet_soln[dirichlet_soln.size()-1]));
    }
}

/**
 * Returns whether anything was changed.
 */
bool contamFinder::reclassify_cells(){
    // Assume uniform global contam rate
    double c = contam_cell_prior;
    bool changed = false;
    
    int n_reassigned = 0;
    
    bool reweight_doublets = (doublet_rate > 0 && doublet_rate < 1);
    
    pair<double, double> betaparams;
    if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
        betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
    }

    map<int, double> priorweights;
    map<int, double>* priorweights_ptr = NULL;
    if (reweight_doublets){
        
        // Notify the populate_llr_table function that we have LLR elements
        // to add into the table
        priorweights_ptr = &priorweights;
        
        // Compute total fraction of RNA per identity (as if bulk)
        // Store total # cells mapped to each identity, counting halves of doublets
        // once and singlets twice
        map<int, int> acounts;

        // Compute total number of cells with each identity
        // Store total # cells mapped to each identity, including doublets
        map<int, int> acounts_incld;
        int atot = 0;
        int atot_incld = 0;
        
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); 
            a != assn.end(); ++a){
            if (a->second >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                if (acounts.count(comb.first) == 0){
                    acounts.insert(make_pair(comb.first, 0));
                }
                if (acounts.count(comb.second) == 0){
                    acounts.insert(make_pair(comb.second, 0));
                }
                acounts[comb.first]++;
                acounts[comb.second]++;
                atot += 2;
                if (acounts_incld.count(a->second) == 0){
                    acounts_incld.insert(make_pair(a->second, 0));
                }
                acounts_incld[a->second]++;
                atot_incld++;
            }

            else{
                if (acounts.count(a->second) == 0){
                    acounts.insert(make_pair(a->second, 0));
                }
                acounts[a->second] += 2;
                atot += 2;   
                if (acounts_incld.count(a->second) == 0){
                   acounts_incld.insert(make_pair(a->second, 0));
                }
                acounts_incld[a->second]++;
                atot_incld++; 
            }
        }
        vector<int> samps;
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            samps.push_back(samp);
        }
        sort(samps.begin(), samps.end());
        
        for (int i = 0; i < samps.size(); ++i){
            int si = samps[i];
            
            // Compute expected fraction of this singlet identity
            // if d is doublet rate and x is the bulk fraction of this
            // identity, then p = (1-d)*s + d*s*s
            double si_p = (double)acounts[si]/(double)atot;
            double prob = (1.0 - doublet_rate)*si_p + doublet_rate*si_p*si_p;
            
            // Add in log likelihood of present number of counts under the model
            double ll = dbinom(atot_incld, acounts_incld[si], prob);
            priorweights.insert(make_pair(si, ll));
            
            for (int j = i + 1; j < samps.size(); ++j){
                int sj = samps[j];
                int sk = hap_comb_to_idx(si, sj, n_samples);

                // Compute expected fraction of this doublet identity
                // if d is doublet rate, x is bulk fraction of ID1 and y is
                // bulk fraction of ID2, then p = 2*d*x*y
                double sj_p = (double)acounts[sj]/(double)atot;
                double prob = 2.0*doublet_rate*si_p*sj_p;
                
                double ll = dbinom(atot_incld, acounts_incld[sk], prob);
                priorweights.insert(make_pair(sk, ll));
            }
            
        }
    }
    else if (false){
        priorweights_ptr = &priorweights;
        pair<double, double> betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
        map<int, double> crmean;
        map<int, int> crtot;
        map<int, vector<double> > crmed;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (contam_rate.count(a->first) > 0){
                if (crmean.count(a->second) > 0){
                    crmean.insert(make_pair(a->second, 0.0));
                    crtot.insert(make_pair(a->second, 0));
                    vector<double> v;
                    crmed.insert(make_pair(a->second, v));
                }
                crmean[a->second] += contam_rate[a->first];
                crtot[a->second]++;
                crmed[a->second].push_back(contam_rate[a->first]);
            }
        }
        for (map<int, double>::iterator crm = crmean.begin(); crm != crmean.end(); ++crm){
            crm->second /= (double)crtot[crm->first];
            /*
            double p = 1.0 - pbeta(crm->second, betaparams.first, betaparams.second);
            if (p == 0){
                p += 1e-6;
            }
            else if (p == 1){
                p -= 1e-6;
            }
            p = log2(p);
            */
            /*
            double med;
            if (crmed[crm->first].size() % 2 == 0){
                double m1 = crmed[crm->first][crmed[crm->first].size()/2];
                double m2 = crmed[crm->first][crmed[crm->first].size()/2-1];
                med = (m1+m2)/2.0;
            }
            else{
                med = crmed[crm->first][(crmed[crm->first].size()-1)/2];
            }
            */
            double p = dbeta(crm->second, betaparams.first, betaparams.second);
            priorweights.insert(make_pair(crm->first, p));
        }
    }
    
    // If the user has set doublet rate to 0 or 1, then instead of doing
    // complicated stuff above to try to enforce correct proportions of 
    // all individuals, we'll just pass that value to populate_llr_table,
    // which will exclude doublets or singlets accordingly.
    
    // Otherwise, no need to include doublet rate since we will already be
    // including it with the "prior_weights" parameter.
    double doub_rate_table = 0.5;
    if (doublet_rate == 0 || doublet_rate == 1){
        doub_rate_table = doublet_rate;
    }

    set<unsigned long> cell_rm;
   
    /*
    vector<pair<double, int> > cpsort;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            cpsort.push_back(make_pair(-cp->second, cp->first));
        }
    }
    sort(cpsort.begin(), cpsort.end());
    
    double newcsum = 0.0;
    double newccount = 0.0;

    int nproc = 0;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        
        vector<int> possibilities;
        vector<double> possibility_lls;
        vector<double> possibility_cs;
    
        string bcstr = bc2str(a->first);
        fprintf(stderr, "%s\n", bcstr.c_str());

        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            fprintf(stderr, "orig (%d %d) %f %f\n", combo.first, combo.second, contam_rate[a->first],
                contam_rate_ll[a->first]);
            if (allowed_ids2.size() == 0 || allowed_ids2.find(combo.first) != allowed_ids2.end()){
                possibilities.push_back(combo.first);
            }
            if (allowed_ids2.size() == 0 || allowed_ids2.find(combo.second) != allowed_ids2.end()){
                possibilities.push_back(combo.second);
            }
            for (int x = 0; x < cpsort.size(); ++x){
                int j = cpsort[x].second;
                if (j != combo.first && j != combo.second){
                    int k1;
                    int k2;
                    if (j < combo.first){
                        k1 = hap_comb_to_idx(j, combo.first, n_samples);
                    }
                    else{
                        k1 = hap_comb_to_idx(combo.first, j, n_samples);
                    }
                    if (j < combo.second){
                        k2 = hap_comb_to_idx(j, combo.second, n_samples);
                    }
                    else{
                        k2 = hap_comb_to_idx(combo.second, j, n_samples);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k1) != allowed_ids2.end()){
                        possibilities.push_back(k1);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k2) != allowed_ids2.end()){
                        possibilities.push_back(k2);
                    }
                }
            }
        }
        else{
            fprintf(stderr, "orig (%d) %f %f\n", a->second, contam_rate[a->first], contam_rate_ll[a->first]);
            //for (int x = 0; x < idx2samp.size(); ++x){
            for (int x = 0; x < cpsort.size(); ++x){
                int j = cpsort[x].second;    
            //int j = idx2samp[x];
                if (j != a->second){
                    int k;
                    if (j < a->second){
                        k = hap_comb_to_idx(j, a->second, n_samples);
                    }
                    else{
                        k = hap_comb_to_idx(a->second, k, n_samples);
                    }
                    if (allowed_ids2.size() == 0 || allowed_ids2.find(k) != allowed_ids2.end()){
                        possibilities.push_back(k);
                    }
                }
            }
        }
        int maxidx = -1;
        double maxll = contam_rate_ll[a->first];
        double maxc = contam_rate[a->first];
        double maxcse = contam_rate_se[a->first];
        for (int x = 0; x < possibilities.size(); ++x){
            int j = possibilities[x];
            bool is_comb = false;
            pair<int, int> comb;
            if (j >= n_samples){
                is_comb = true;
                comb = idx_to_hap_comb(j, n_samples);
            }
            vector<double> n;
            vector<double> k;
            vector<double> p_c;
            vector<double> p_e;
            for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac = 
                indv_allelecounts[a->first].begin(); ac != indv_allelecounts[a->first].end();
                ++ac){
                if ((!is_comb && ac->first.first == j) || (is_comb && ac->first.first == comb.first)){
                    for (map<pair<int, int>, pair<float, float> >::iterator ac2 = 
                        ac->second.begin(); ac2 != ac->second.end(); ++ac2){
                        bool breakout = false;
                        double p_e_this;
                        bool keep = false;
                        if ((!is_comb && ac2->first.first == -1)){
                            p_e_this = (double)ac->first.second/2.0;
                            keep = true;
                            breakout = true;
                        }
                        else if ((is_comb && ac2->first.first == comb.second)){
                            p_e_this = (double)(ac->first.second + ac2->first.second)/4.0;
                            keep = true;
                        }
                        if (keep){
                            double p_c_this = amb_mu[ac->first][ac2->first];
                            n.push_back(ac2->second.first + ac2->second.second);
                            k.push_back(ac2->second.second);
                            p_c.push_back(p_c_this);
                            p_e.push_back(p_e_this);
                            
                        }
                        if (breakout){
                            break;
                        }
                    }
                }
            }
            if (n.size() == 0){
                continue;
            }
            optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
            if (num_threads > 1){
                c_cell.set_threads(num_threads);
            }
            c_cell.constrain_01();
            c_cell.add_data("n", n);
            c_cell.add_data("k", k);
            c_cell.add_data("p_c", p_c);
            c_cell.add_data("p_e", p_e);
            if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                c_cell.add_beta_prior(betaparams.first, betaparams.second);
            }
            c_cell.set_maxiter(-1);
            double c_cell_map = c_cell.solve(contam_rate[a->first],1);
            //double c_cell_map = c_cell.solve(0,1);
            if (is_comb){
                fprintf(stderr, "  (%d %d) %f %f\n", comb.first, comb.second, c_cell_map,
                    c_cell.log_likelihood);
            }
            else{
                fprintf(stderr, "  (%d) %f %f\n", j, c_cell_map, c_cell.log_likelihood);
            }
            if (c_cell.root_found){
                if (maxll == 0 || c_cell.log_likelihood > maxll){
                    maxll = c_cell.log_likelihood;
                    maxidx = j;
                    maxc = c_cell_map;
                    maxcse = c_cell.se;
                }
                else if (a->second < n_samples || x > 1){
                    //break;
                }    
            }
        }
        if (maxll != 0 && maxidx != -1){
            ++n_reassigned;
            fprintf(stderr, "reassign %d / %d\n", n_reassigned, nproc);
            // Update
            contam_rate_ll[a->first] = maxll;
            contam_rate[a->first] = maxc;
            contam_rate_se[a->first] = maxcse;
            assn[a->first] = maxidx;
            
            newcsum += maxc;
            newccount++;
            
            fprintf(stderr, " c = %f\n", newcsum / newccount);
        }
        nproc++;
    }
    fprintf(stderr, " %d cells reassigned\n", n_reassigned);
    return n_reassigned > 0;
    */ 

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); 
        a != assn.end(); ++a){
        
        // Get a table of log likelihood ratios between every possible
        // pair of identities
        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        
        //c = contam_rate[a->first]; 
        bool success = populate_llr_table(indv_allelecounts[a->first], llrs, tab, n_samples, 
            allowed_ids, allowed_ids2, doub_rate_table, e_r, e_a, priorweights_ptr,
            true, contam_cell_prior, 0, &amb_mu);

        if (success){
            
            // If we reassigned an identity and it's different than the old one,
            // also infer a new contamination rate. Accept the change if the
            // overall log likelihood of the new identity + new contamination rate
            // beats the log likelihood of the old identity + old contamination rate.

            int a_new;
            double llr_new; 
            tab.get_max(a_new, llr_new);
            if (a_new != -1 && llr_new > 0){
                if (a_new != a->second){
                    
                    // Get new contam rate conditional on candidate new identity
                    vector<double> n;
                    vector<double> k;
                    vector<double> p_e;
                    vector<double> p_c;

                    if (a_new >= n_samples){
                        pair<int, int> comb = idx_to_hap_comb(a_new, n_samples);
                        for (int x = 0; x <= 2; ++x){
                            pair<int, int> k1 = make_pair(comb.first, x);
                            for (int y = 0; y <= 2; ++y){
                                pair<int, int> k2 = make_pair(comb.second, y);
                                if (true){
                                //if (x != y){
                                    double ref = indv_allelecounts[a->first][k1][k2].first;
                                    double alt = indv_allelecounts[a->first][k1][k2].second;
                                    n.push_back(ref+alt);
                                    k.push_back(alt);
                                    p_e.push_back(adjust_p_err((double)(x+y)/4.0, e_r, e_a));
                                    p_c.push_back(amb_mu[k1][k2]);
                                }
                            }
                        }
                    }
                    else{
                        pair<int, int> nullkey = make_pair(-1, -1);
                        for (int x = 0; x <= 2; ++x){
                            pair<int, int> key = make_pair(a_new, x);
                            double ref = indv_allelecounts[a->first][key][nullkey].first;
                            double alt = indv_allelecounts[a->first][key][nullkey].second;
                            n.push_back(ref+alt);
                            k.push_back(alt);
                            p_e.push_back(adjust_p_err((double)x/2.0, e_r, e_a));
                            p_c.push_back(amb_mu[key][nullkey]);
                        }
                    }
                    optimML::brent_solver c_cell(ll_c, dll_dc, d2ll_dc2);
                    
                    if (num_threads > 1){
                        c_cell.set_threads(num_threads);
                    }
                    c_cell.add_data("n", n);
                    c_cell.add_data("k", k);
                    c_cell.add_data("p_e", p_e);
                    c_cell.add_data("p_c", p_c);
                    c_cell.constrain_01();
                    if (contam_cell_prior > 0 && contam_cell_prior_var > 0){
                        c_cell.add_beta_prior(betaparams.first, betaparams.second);
                    }
                    c_cell.set_maxiter(-1);
                    double c_cell_map = c_cell.solve(0,1);
                    
                    // Only accept the change if the new assignment + contam rate inference has a higher
                    // log likelihood than the older assignment + contam rate inference
                    //if (c_cell.root_found && c_cell_map < contam_rate[a->first] && 
                    //    (contam_rate_ll[a->first] == 0 || 
                    //    c_cell.log_likelihood > contam_rate_ll[a->first])){
                    
                    //if (c_cell.root_found &&
                    //    (contam_rate_ll[a->first] == 0 || c_cell.log_likelihood < contam_rate_ll[a->first])){
  
                    if (c_cell.root_found &&
//                        c_cell_map >= contam_rate[a->first] &&  
                        (contam_rate_ll[a->first] == 0 || c_cell.log_likelihood >= contam_rate_ll[a->first])){
                        n_reassigned++;
                        
                        //fprintf(stdout, "%s\t%f\t%f\t%f\t%f\n", bc2str(a->first).c_str(),
                        //    contam_rate[a->first], c_cell_map, contam_rate_ll[a->first],
                        //    c_cell.log_likelihood);


                        changed = true;
                        a->second = a_new;
                        assn_llr[a->first] = llr_new;
                        contam_rate[a->first] = c_cell_map;
                        contam_rate_se[a->first] = c_cell.se;
                        contam_rate_ll[a->first] = c_cell.log_likelihood;
                    }
                }
            }
            else{
                cell_rm.insert(a->first);        
            }
        }
        else{
            // Keep original? Delete original?
            cell_rm.insert(a->first);
        }
    }
    
    fprintf(stderr, "  %d cells reassigned\n", n_reassigned);
    fprintf(stderr, "  %ld cells removed\n", cell_rm.size());

    if (cell_rm.size() > 0){
        changed = true;
    }
    for (set<unsigned long>::iterator rm = cell_rm.begin(); rm != cell_rm.end(); ++rm){
        assn.erase(*rm);
        assn_llr.erase(*rm);
        contam_rate.erase(*rm);
        contam_rate_se.erase(*rm);
    }
    return changed;
}

/**
 * Get best estimates of (residual) reference & alt allele misreading
 * error rates. This should be done after estimating all other parameters,
 * to see how close to zero these values are. High values indicate 
 * discordance between the data and variant calls not reflected in the
 * contamination model.
 */
pair<double, double> contamFinder::est_error_rates(bool init){
    
    vector<double> n;
    vector<double> k;
    vector<double> c;
    vector<double> p_e;
    vector<double> p_c; 
    vector<double> weights;

    for (map<unsigned long, vector<int> >::iterator ci = cell_to_idx.begin();
        ci != cell_to_idx.end(); ++ci){
        
        if (assn_llr.count(ci->first) > 0 && (init || contam_rate.count(ci->first) > 0)){
            for (vector<int>::iterator i = ci->second.begin(); i != 
                ci->second.end(); ++i){
                if (init){
                   c.push_back(0.0);
                }
                else{
                    c.push_back(contam_rate[ci->first]);
                }
                n.push_back(n_all[*i]);
                k.push_back(k_all[*i]);
                p_e.push_back(p_e_all[*i]);
                double weight = 1.0;
                if (weighted){
                    weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
                }
                weights.push_back(weight);
                if (init){
                    p_c.push_back(0.0);
                }
                else{
                    p_c.push_back(amb_mu[type1_all[*i]][type2_all[*i]]);
                }
            }
        }
    }
    
    optimML::multivar_ml_solver solver({e_r, e_a}, ll_err_rates, dll_err_rates);
    if (num_threads > 1){
        solver.set_threads(num_threads);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("c", c);
    solver.add_data("p_e", p_e);
    solver.add_data("p_c", p_c);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_weights(weights);

    double this_e_r = e_r;
    double this_e_a = e_a;

    try{
        solver.solve();
    
        this_e_r = solver.results[0];
        this_e_a = solver.results[1];
    
    }
    catch (int errcode){
        if (errcode != optimML::OPTIMML_MATH_ERR){
            fprintf(stderr, "unknown error\n");
            exit(1);
        }
        fprintf(stderr, "Error inferring error rates; keeping initial values.\n");
    }
    return make_pair(this_e_r, this_e_a);
}

/**
 * Compute log likelihood of data set with current parameters
 * Assumes pre-compiled data. This will need to be re-generated after
 * assignments are changed via reclassify_cells().
 */
double contamFinder::compute_ll(){
    // Compute log likelihood of data set
    double loglik = 0.0;
    pair<double, double> betaparams = beta_moments(contam_cell_prior, contam_cell_prior_var);
    for (int i = 0; i < n_all.size(); ++i){
        if (contam_rate.count(idx_to_cell[i]) > 0){
            double c = contam_rate[idx_to_cell[i]];
            double p_e = adjust_p_err(p_e_all[i], e_r, e_a);
            double p_c = amb_mu[type1_all[i]][type2_all[i]];
            double binom_p = (1.0-c)*p_e + c*p_c;
            loglik += logbinom(n_all[i], k_all[i], binom_p);
            
            // Add in prior prob of cell contam rates
            //double a = (c - contam_cell_prior) / contam_cell_prior_se;
            //loglik += (-0.5*pow(a, 2)) - log(contam_cell_prior_se) - log(sqrt(2*3.14159265358979));
            loglik += (dbeta(c, betaparams.first, betaparams.second)/log2(exp(1.0)));
        }
    }
    return loglik;
}

/**
 * Find maximum likelihood estimates of all parameters of interest.
 * Return overall log likelihood.
 */
void contamFinder::fit(){
    // Begin with minimum estimate of c (global contamination rate)
    // if not provided from a previous run
    if (c_init <= 0){
        // Get (estimated/averaged) min bound on c
        c_init = this->est_min_c();
    }
    double ll_init = init_params(c_init);
    fprintf(stderr, "Initial global contamination rate = %f\n", c_init);
    
    double dummy = -1.0;
    this->est_contam_cells_global();
    double loglik = this->update_amb_prof_mixture(false, dummy, true);

    // Once to set prior params
    this->est_contam_cells();
    // Again to use empirical Bayes shrinkage
    this->est_contam_cells();
    
    // Allow ambient prof to update without considering individuals of origin
    //this->update_ambient_profile(false);
    
    // Update mixture components using per-cell contam estimates
    //this->update_amb_prof_mixture(false, dummy, false);
    
    if (!skip_reassign){ 
        // Update cell identities, considering contamination profile 
        fprintf(stderr, "Reclassifying cells...\n");
        bool reclassified = this->reclassify_cells();
        
        if (reclassified){
            // Allows log likelihood computation
            clear_data();
            compile_data(assn, indv_allelecounts);
        }
    }
    /* 
    // Check how low the error rates have dropped after modeling contamination
    pair<double, double> err_final = this->est_error_rates(false);
    fprintf(stderr, "Residual error rates:\n");
    fprintf(stderr, "  Reference alleles: %f\n", err_final.first);
    fprintf(stderr, "  Alt alleles: %f\n", err_final.second);
    */
}

