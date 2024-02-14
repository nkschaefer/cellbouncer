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
#include <mixtureDist/functions.h>
#include <optimML/brent.h>
#include <optimML/multivar_ml.h>
#include "common.h"
#include "robin_hood.h"
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
 * Adjust an expected allele matching rate given the ref and alt allele
 * misreading errors.
 */
double adjust_p_err(double p, double e_r, double e_a){
    return p - p*e_a + (1-p)*e_r;
}

/**
 * Log PDF of binomial distribution wrt n, k, p
 */
double logbinom(double n, double k, double p){
    double ll = k * log(p) + (n-k)*log(1.0-p);
    // Compute log binomial coefficient
    if (k < n && k != 0){
        // Use Stirling's approximation
        double logn = log(n);
        double logk = log(k);
        double logn_k = log(n-k);
        ll += n*logn - k*logk - (n-k)*logn_k + 0.5*(logn - logk - logn_k - log(2*M_PI));
    }
    return ll;
}

/**
 * Integrate binomial log likelihood wrt p_c (expected alt allele match rate from contamination)
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
 */
double d_integral_ll_pc_dc(double n, double k, double c, double p_e, double p_c){
    double x = (p_e-1)*(k-n)*log((c-1)*p_e - c*p_c + 1);
    double y = k*p_e*log(-c*p_e + c*p_c + p_e);
    double z = c*n*(p_c-p_e);
    return (x - y + z)/(c*c);
}

double integral_ll_c(double n, double k, double c, double p_e, double p_c){
    double denom = p_e - p_c;
    double x = -c*(p_e - p_c)*(n - binom_coef_log(n,k)/log2(exp(1)));
    double y = (k-n)*((c-1)*p_e - c*p_c + 1)*log((c-1)*p_e - c*p_c + 1);
    double z = k*((c-1)*p_e - c*p_c)*log(-c*p_e + c*p_c + p_e);
    
    return (x - y + z)/denom;
}

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
double ll_ambmu(vector<double>& params, const map<string, double>& data_d, 
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

double ll_ambmu2(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = params[params.size()-1];
    double binom_p = (1.0-c) * p_e + c*p_c;
    return logbinom(n, k, binom_p); 
}

/**
 * First derivative of log likelihood for estimating best set
 * of p_c parameters in ambient RNA
 */
void dll_ambmu(vector<double>& params, const map<string, double>& data_d, 
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

void dll_ambmu2(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_e = data_d.at("p_e");
    int amb_idx = data_i.at("ef_idx");
    double p_c = params[amb_idx];
    double c = params[params.size()-1];
    double binom_p = (1.0-c) * p_e + c*p_c;

    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    results[amb_idx] = dy_dp * c;
    results[params.size()-1] = dy_dp * (p_c - p_e);
}

/**
 * Log likelihood function for finding best e_r & e_a
 */
double ll_err_rates(vector<double>& params, const map<string, double>& data_d, 
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
    return logbinom(n, k, binom_p);
}

/**
 * First derivative log likelihood function for finding best e_r & e_a
 */
void dll_err_rates(vector<double>& params, const map<string, double>& data_d, 
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
    
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    
    results[0] = dy_dp * (1 - p_e + c*(p_e - p_c));
    results[1] = dy_dp * (c*(p_e - p_c) - p_e);
}

/**
 * Log likelihood for computing expected p_c values in real cells given each
 * possible identity
 */
double ll_expfrac(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_c = data_d.at("p_c");
    double c = data_d.at("c");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");

    if (idx2 == -1){
        double binom_p = (1.0-c)*params[idx1] + p_c*c;
        return logbinom(n, k, binom_p);
    }
    else{
        double binom_p = (1.0 - c) * (0.5*params[idx1] + 0.5*params[idx2]) + p_c*c;
        return logbinom(n, k, binom_p);
    }
}

/**
 * Derivative of log likelihood for computing expected p_c values in real cells
 * given each possible identity
 */
void dll_expfrac(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    double p_c = data_d.at("p_c");
    double c = data_d.at("c");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    
    double dp_dparam;
    double binom_p;
    if (idx2 == -1){
        binom_p = (1.0-c)*params[idx1] + p_c*c;
        dp_dparam = 1.0 - c;
    }
    else{
        binom_p = (1.0 - c) * (0.5*params[idx1] + 0.5*params[idx2]) + p_c*c;
        dp_dparam = 0.5 - 0.5*c;
    }
    
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    results[idx1] = dy_dp * dp_dparam;
    if (idx2 != -1){
        results[idx2] = dy_dp * dp_dparam;
    }
}

/**
 * Log likelihood function for modeling ambient RNA as a mixture of 
 * individuals
 */
double ll_as_mixture(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double p_c = params[0];
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double binom_p = (1-c)*p_e + c*p_c;
    return logbinom(n, k, binom_p);
}

/**
 * First derivative log likelihood function for modeling ambient RNA
 * as a mixture of individuals
 */
void dll_as_mixture(vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double p_c = params[0];
    double p_e = data_d.at("p_e");
    double n = data_d.at("n");
    double k = data_d.at("k");
    double c = data_d.at("c");

    double binom_p = (1-c)*p_e + c*p_c;
    
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    // Only one variable to deal with
    double dp_dx = c;
    results[0] = dy_dp * dp_dx;
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
    int n_samples){
    
    // Copy external data structures that we will need in the future 
    this->assn = assn;
    this->assn_llr = assn_llr;
    this->indv_allelecounts = indv_allelecounts;
    this->n_samples = n_samples;
    this->n_mixprop_trials = 10;
     
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
        }
        id_llrsum[a->second] += assn_llr[a->first];
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
    
    this->max_cells_expfrac = 50;

    this->inter_species = false;

    this->contam_cell_prior = -1;
    this->contam_cell_prior_se = -1;

    this->delta_thresh = 0.1;
    this->maxits = 100;
    
    this->model_mixprop = true;
    
    this->weighted = false;

    // Compile data in the format needed by other functions
    this->compile_data(assn, indv_allelecounts);
}

/**
 * Functions to set parameter values
 */

void contamFinder::set_error_rates(double error_ref, double error_alt){
    this->e_r = error_ref;
    this->e_a = error_alt;
}

void contamFinder::model_other_species(){
    this->inter_species = true;
}

void contamFinder::model_single_species(){
    this->inter_species = false;
}

void contamFinder::model_mixture(){
    this->model_mixprop = true;
}

void contamFinder::skip_model_mixture(){
    this->model_mixprop = false;
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

void contamFinder::max_cells_for_expfracs(int mc){
    this->max_cells_expfrac = mc;
}

void contamFinder::use_weights(){
    this->weighted = true;
}

void contamFinder::no_weights(){
    this->weighted = false;
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
                double expected = ac->first.second / 2.0;
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
 * Get an initial estimate of the global contamination rate, as well as
 * all p_c parameters (expected match rate for alt alleles in ambient RNA)
 *
 * Solve jointly, using BFGS, with initial guesses of 0.001 for everything.
 */
void contamFinder::solve_params_init(){
    
    vector<pair<int, int> > idx2expfrac1;
    vector<pair<int, int> > idx2expfrac2;
    
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> p_c;
    vector<int> efp_idx;
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
                if (idx_to_cell.count(*i) > 0){ 
                    
                    idx2idx.insert(make_pair(*i, n.size()));    
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                    efp_idx.push_back(ef_param_idx);
                    
                    double weight = 1.0;
                    if (weighted){ 
                        // Weight by fraction of LLR sum of all cells assigned to this ID
                        // assigned to this cell
                        weight = assn_llr[idx_to_cell[*i]] / id_llrsum[assn[idx_to_cell[*i]]];
                    }
                    weights.push_back(weight);
                }
            }        
        }
    }
    efparams.push_back(0.01);
    optimML::multivar_ml_solver solver(efparams, ll_ambmu2, dll_ambmu2);
    for (int i = 0; i < efparams.size(); ++i){
        solver.constrain_01(i);
    }
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_data("ef_idx", efp_idx);
    solver.add_weights(weights);
    solver.solve();
    
    for (int i = 0; i < efparams.size()-1; ++i){
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
    fprintf(stderr, "Init global contamination estimate = %f\n", 
        solver.results[solver.results.size()-1]);
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
        c_cell.add_data("n", n);
        c_cell.add_data("k", k);
        c_cell.add_data("p_e", p_e);
        c_cell.add_data("p_c", p_c);
        c_cell.constrain_01();
        if (contam_cell_prior > 0 && contam_cell_prior_se > 0){
            // The last two parameters indicate that this is a truncated normal distribution
            // on the interval [0,1]
            c_cell.add_normal_prior(contam_cell_prior, contam_cell_prior_se, 0, 1);
        }
        double c_cell_map = c_cell.solve(0,1);
    
        // Only store estimates if a maximum-likelihood estimate was found in the range
        // (0,1)
        if (c_cell.root_found){
            contam_rate.emplace(ci->first, c_cell_map);
            cell_c_maps.push_back(c_cell_map);
            if (weighted){
                double weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
                cell_c_llr.push_back(weight);
            }
            if (c_cell.se_found){
                contam_rate_se.emplace(ci->first, c_cell.se);
            }
        }
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
        fprintf(stderr, "Contamination rate estimates:\n");
        fprintf(stderr, "  Mean: %f Std dev: %f\n", mu_var.first, sqrt(mu_var.second));
        contam_cell_prior = mu_var.first;
        contam_cell_prior_se = sqrt(mu_var.second);
    }
}

/**
 * After estimating contamination rate per cell, finds the likeliest ambient RNA
 * profile (set of p_c parameters for each type of SNP) based on the current
 * contamination rate estimates.
 */
double contamFinder::update_ambient_profile(){
    
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
                if (idx_to_cell.count(*i) > 0 && 
                    contam_rate.count(idx_to_cell[*i]) > 0){ 
                    
                    idx2idx.insert(make_pair(*i, n.size()));    
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                    efp_idx.push_back(ef_param_idx);
                    c.push_back(contam_rate[idx_to_cell[*i]]);
                    
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
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("c", c);
    solver.add_data("p_e", p_e);
    solver.add_data("p_c", p_c);
    solver.constrain_01(0);
    solver.constrain_01(1);
    solver.add_weights(weights);
    solver.solve();
    
    double this_e_r = solver.results[0];
    double this_e_a = solver.results[1];
    
    //fprintf(stderr, "Residual error rates:\n");
    //fprintf(stderr, "  Reference alleles: %f\n", this_e_r);
    //fprintf(stderr, "  Alternate alleles: %f\n", this_e_a);
    
    return make_pair(this_e_r, this_e_a);
}

/**
 * With the fit model, get expected alt allele matching rates of each type
 * of allele given each possible identity.
 */
void contamFinder::compute_expected_fracs_all_id(){
    expfracs.clear();
    
    // Figure out which cells to use -- take a maximum possible number per identity.
    // This will save on compute time (using the whole data set can take a very long
    // time, especially with large numbers of cells)

    set<unsigned long> cells_keep;
    if (max_cells_expfrac > 0){
        map<int, vector<pair<double, unsigned long> > > id2cells;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (id2cells.count(a->second) == 0){
                vector<pair<double, unsigned long> > v;
                id2cells.insert(make_pair(a->second, v));
            }
            id2cells[a->second].push_back(make_pair(-assn_llr[a->first], a->first));
        }
        for (map<int, vector<pair<double, unsigned long> > >::iterator i2c = id2cells.begin();
            i2c != id2cells.end(); ++i2c){
            // Sort by decreasing log likelihood ratio of assignment
            sort(i2c->second.begin(), i2c->second.end());
            for (int i = 0; i < max_cells_expfrac; ++i){
                cells_keep.insert(i2c->second[i].second);
            }
        }
    }

    // Set up parameters & compute initial values
    vector<double> params;
    vector<double> params_tot;

    vector<int> param2id;
    vector<pair<int, int> > param2type;
    map<pair<int, int>, map<int, int> > param_lookup;

    vector<double> n;
    vector<double> k;
    vector<double> c;
    vector<double> p_c;
    vector<int> pidx1;
    vector<int> pidx2;
    vector<double> pweight;
    
    pair<int, int> nullkey = make_pair(-1, -1);

    // Initialize stuff
    for (int i = 0; i < n_samples; ++i){
        if (allowed_ids.size() == 0 || allowed_ids.find(i) != allowed_ids.end()){
            for (int n = 0; n <= 2; ++n){
                pair<int, int> key = make_pair(i, n);
                if (param_lookup.count(key) == 0){
                    map<int, int> m;
                    param_lookup.insert(make_pair(key, m));
                }
                for (int j = 0; j < n_samples; ++j){
                    if (allowed_ids.size() == 0 || allowed_ids.find(j) != allowed_ids.end()){
                        if (param_lookup[key].count(j) == 0){
                            int idx = param2id.size();
                            param_lookup[key].insert(make_pair(j, idx));
                            param2id.push_back(j);
                            params.push_back(0.0);
                            params_tot.push_back(0.0);
                            param2type.push_back(key);
                        }
                    }
                }
            }
        }
    }
    
    // Get data
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); 
        a != assn.end(); ++a){
        
        // Skip cells that are not in the top however many we want
        if (cells_keep.size() > 0 && cells_keep.find(a->first) == cells_keep.end()){
            continue;
        }

        if (contam_rate.count(a->first) > 0){
            
            unsigned long cell = a->first;
            int cell_id = a->second;
            
            double weight = 1.0;
            if (weighted){
                //weight = assn_llr[a->first] / id_llrsum[a->second];
                // Don't normalize this because it's already by ID
                weight = assn_llr[a->first];
            }
            double contam_cell = contam_rate[a->first];

            bool is_combo = false;
            pair<int, int> combo;
            if (a->second >= n_samples){
                is_combo = true;
                combo = idx_to_hap_comb(a->second, n_samples);
                if (allowed_ids.size() > 0 && (allowed_ids.find(combo.first) == allowed_ids.end() ||
                    allowed_ids.find(combo.second) == allowed_ids.end())){
                    continue;
                }
            }
            else if (allowed_ids.size() > 0 && allowed_ids.find(a->second) == allowed_ids.end()){
                continue;
            }
            

            for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator x = 
                indv_allelecounts[a->first].begin(); x != indv_allelecounts[a->first].end(); ++x){
                
                double ref = x->second[nullkey].first;
                double alt = x->second[nullkey].second;
                double frac = alt/(ref+alt);
                    
                for (int i = 0; i < n_samples; ++i){
                    if (allowed_ids.size() == 0 || allowed_ids.find(i) != allowed_ids.end()){
                        
                        if (is_combo){
                            int idx1 = param_lookup[x->first][combo.first];
                            int idx2 = param_lookup[x->first][combo.second];
                            pidx1.push_back(idx1);
                            pidx2.push_back(idx2);
                            params[idx1] += frac * weight;
                            params_tot[idx1] += weight;
                            params[idx2] += frac * weight;
                            params_tot[idx2] += weight;
                            
                        }
                        else{
                            int idx = param_lookup[x->first][cell_id];
                            params[idx] += frac * weight;
                            params_tot[idx] += weight;
                            pidx1.push_back(idx);
                            pidx2.push_back(-1);
                            
                        }

                        n.push_back(ref+alt);
                        k.push_back(alt);
                        c.push_back(contam_cell);
                        pweight.push_back(weight);
                        p_c.push_back(amb_mu[x->first][nullkey]);
                    }
                }
            }
        }
    }

    // Get mean of each estimate
    for (int i = 0; i < params.size(); ++i){
        params[i] /= params_tot[i];
    }

    // Set up equation solver
    optimML::multivar_ml_solver solver(params, ll_expfrac, dll_expfrac);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("c", c);
    solver.add_data("p_c", p_c);
    solver.add_weights(pweight);
    solver.add_data("idx1", pidx1);
    solver.add_data("idx2", pidx2);
    for (int i = 0; i < params.size(); ++i){
        solver.constrain_01(i);
    }
    solver.solve();

    // Get results back
    for (int i = 0; i < params.size(); ++i){
        double frac = solver.results[i];
        pair<int, int> type = param2type[i];
        int id = param2id[i];
        if (expfracs.count(type) == 0){
            map<int, float> m;
            expfracs.insert(make_pair(type, m));
        }
        if (expfracs[type].count(id) == 0){
            expfracs[type].insert(make_pair(id, frac));
        }
        else{
            expfracs[type][id] = frac;
        }
    }
}

/**
 * Models the current ambient RNA profile as a mixture of individuals. Requires
 * that the amb_mu fractions have been estimated and that expected alt allele fractions
 * of each type under every possible identity have also been computed 
 * (compute_expected_fracs_all_id).
 *
 * Populates the contam_prof map, accessible outside the class.
 *
 */
double contamFinder::model_as_mixture(){

    vector<vector<double> > mixfrac_exp; 
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    vector<int> idx;

    for (map<pair<int, int>, map<pair<int, int>, vector<int> > >::iterator x = expfrac_to_idx.begin();
       x != expfrac_to_idx.end(); ++x){
        for (map<pair<int, int>, vector<int> >::iterator y = x->second.begin(); 
            y != x->second.end(); ++y){

            vector<double> mixfrac_row;
            for (int i = 0; i < idx2samp.size(); ++i){
                int samp = idx2samp[i];
                double expected;
                if (y->first.first == -1){
                    expected = expfracs[x->first][samp];
                }
                else{
                    double ef1 = expfracs[x->first][samp];
                    double ef2 = expfracs[y->first][samp];
                    expected = 0.5*ef1 + 0.5*ef2;
                }

                // These values are observed -- i.e. they already include errors.
                // Therefore, no need to adjust for errors.  
                
                mixfrac_row.push_back(expected);
            
            }
            if (inter_species){
                double expected = adjust_p_err(0.0, e_r, e_a);
                mixfrac_row.push_back(expected);
            }
            for (vector<int>::iterator i = y->second.begin(); i != y->second.end(); ++i){
                if (contam_rate.count(idx_to_cell[*i]) > 0){
                    idx.push_back(*i);
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                    c.push_back(contam_rate[idx_to_cell[*i]]);
                    mixfrac_exp.push_back(mixfrac_row);
                }
            }
        }
    }

    // Set up solver
    optimML::multivar_ml_solver mix({}, ll_as_mixture, dll_as_mixture);
    mix.add_data("n", n);
    mix.add_data("k", k);
    mix.add_data("c", c);
    mix.add_data("p_e", p_e);
    mix.add_mixcomp(mixfrac_exp);
    
    // Since MLE finding in this case can be very sensitive to initial
    // values, we will try once starting from an even mixture of cells.
    // We will then try again starting from individual proportions matching
    // their frequency in the group of cell->ID assignments. Finally, we
    // will randomly shuffle the starting proportions a set number of times.
    
    // We will end by taking the solution with the highest likelihood of
    // all those trials. 
    
    vector<double> solution_ll;
    vector<vector<double> > solution_cp;
    
    // Solve once using an even pool of each individual as initial guess
    mix.solve();

    solution_ll.push_back(mix.log_likelihood);
    solution_cp.push_back(mix.results_mixcomp);

    if (!inter_species){
        map<int, int> assncount;
        int assntot;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (assncount.count(a->second) == 0){
                assncount.insert(make_pair(a->second, 0));
            }
            assncount[a->second]++;
            assntot++;
        }
        vector<double> starting_props;
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            starting_props.push_back((double)assncount[samp]/(double)assntot);
        }

        // Solve another time using a pool with frequencies based on frequencies
        // of assignment
        mix.add_mixcomp_fracs(starting_props);
        mix.solve();
        solution_ll.push_back(mix.log_likelihood);
        solution_cp.push_back(mix.results_mixcomp);
    }

    // Finally, do a set number of random trials.
    for (int i = 0; i < n_mixprop_trials; ++i){
        mix.randomize_mixcomps();
        mix.solve();
        solution_ll.push_back(mix.log_likelihood);
        solution_cp.push_back(mix.results_mixcomp);
    }
    
    // Find maximum likelihood solution
    int max_sol_idx = -1;
    double max_ll;
    for (int i = 0; i < solution_ll.size(); ++i){
        if (max_sol_idx == -1 || solution_ll[i] > max_sol_idx){
            max_sol_idx = i;
            max_ll = solution_ll[i];
        }
    }

    // Update contamination profile
    contam_prof.clear();
    for (int i = 0; i < idx2samp.size(); ++i){
        int samp = idx2samp[i];
        contam_prof.insert(make_pair(samp, solution_cp[max_sol_idx][i]));
    }
    if (inter_species){
        contam_prof.insert(make_pair(-1, solution_cp[max_sol_idx][solution_cp[max_sol_idx].size()-1]));
    }

    return max_ll;
}

/**
 * Compute log likelihood of data set with current parameters
 */
double contamFinder::compute_ll(){
    // Compute log likelihood of data set
    double loglik = 0.0;
    for (int i = 0; i < n_all.size(); ++i){
        if (contam_rate.count(idx_to_cell[i]) > 0){
            double c = contam_rate[idx_to_cell[i]];
            double p_e = adjust_p_err(p_e_all[i], e_r, e_a);
            double p_c = amb_mu[type1_all[i]][type2_all[i]];
            double binom_p = (1.0-c)*p_e + c*p_c;
            loglik += logbinom(n_all[i], k_all[i], binom_p);
        }
    }
    return loglik;
}

/**
 * Find maximum likelihood estimates of all parameters of interest.
 * Return overall log likelihood.
 */
double contamFinder::fit(){
   
    fprintf(stderr, "Get initial parameter estimates...\n"); 
    // Get initial estimates for global contamination rate and all
    // p_c parameters
    solve_params_init();
    fprintf(stderr, "done\n");

    double delta = 999;
    nits = 0;
    double llprev = 0.0;
    
    while (delta > this->delta_thresh && nits < this->maxits){
        nits++;
        
        // Store backup of previous copies of parameter values in case
        // log likelihood decreases this iteration
        
        amb_mu_prev = amb_mu;
        contam_rate_prev = contam_rate;
        contam_rate_prev_se = contam_rate_se;
        
        // Get best estimates of contamination per cell, given current 
        // contamination profile
        this->est_contam_cells();
        
        // Get best estimate of contamination profile, given current
        // estimates of contamination per cell
        double loglik = this->update_ambient_profile();
        
        if (llprev != 0.0){
            delta = loglik - llprev;
            if (delta < 0){
                // Restore previous estimates
                contam_rate = contam_rate_prev;
                contam_rate_se = contam_rate_prev_se;
                amb_mu = amb_mu_prev;
                loglik = llprev;
                delta = 0.0;
            }
        }
        
        fprintf(stderr, " -- Log likelihood: %f\n", loglik);
        llprev = loglik;
    } 
    
    // Check how low the error rates have dropped after modeling contamination
    pair<double, double> err_final = this->est_error_rates(false);
    fprintf(stderr, "Residual error rates:\n");
    fprintf(stderr, "  Reference alleles: %f\n", err_final.first);
    fprintf(stderr, "  Alt alleles: %f\n", err_final.second);

    if (model_mixprop){
        this->compute_expected_fracs_all_id();
        double mixll = this->model_as_mixture();
    }

    return llprev;
}

