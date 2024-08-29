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

double ll_ambmu2(const vector<double>& params, const map<string, double>& data_d, 
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

void dll_ambmu2(const vector<double>& params, const map<string, double>& data_d, 
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
    
    results[0] = dy_dp * (1 - p_e + c*(p_e - p_c));
    results[1] = dy_dp * (c*(p_e - p_c) - p_e);
    
    //fprintf(stderr, "p_c %f p_e %f n %f k %f c %f | %f %f\n", p_c, p_e, n, k, c, params[0], params[1]); 
    //fprintf(stderr, " -- dy_dp %f [0] %f [1] %f\n", dy_dp, (1-p_e + c*(p_e-p_c)), (c*(p_e-p_c) - p_e));
}

/**
 * Log likelihood for computing expected p_c values in real cells given each
 * possible identity
 */
double ll_expfrac(const vector<double>& params, const map<string, double>& data_d, 
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
void dll_expfrac(const vector<double>& params, const map<string, double>& data_d, 
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
double ll_as_mixture(const vector<double>& params, const map<string, double>& data_d, 
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
void dll_as_mixture(const vector<double>& params, const map<string, double>& data_d, 
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
    
    this->inter_species = false;

    this->contam_cell_prior = -1;
    this->contam_cell_prior_se = -1;

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
                        
                        if (omit_hets && (ac->first.second == 1 || ac2->first.second == 1)){
                            continue;
                        }

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
            if (ac->first.first == ident && (!omit_hets || ac->first.second != 1)){
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
 * Get an initial estimate of the global contamination rate, as well as
 * all p_c parameters (expected match rate for alt alleles in ambient RNA)
 *
 * Solve jointly, using BFGS, with initial guesses of 0.001 for everything.
 */
double contamFinder::solve_params_init(){
    
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
            
            if (omit_hets && ei2->first.second == 1){
                continue;
            }

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

    // Add in global contam rate
    if (c_init_global > 0){
        efparams.push_back(c_init_global);
    }
    else{
        efparams.push_back(0.01);
    }
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
    
    return solver.log_likelihood;
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
        c_cell.set_maxiter(-1);
        double c_cell_map = c_cell.solve(0,1);
        
        // Set SE to 0 if we did not find a maximum-likelihood estimate in the range (0,1) 
        contam_rate.emplace(ci->first, c_cell_map);
        cell_c_maps.push_back(c_cell_map);
        if (weighted){
            double weight = assn_llr[ci->first] / id_llrsum[assn[ci->first]];
            cell_c_llr.push_back(weight);
        }
        if (c_cell.root_found && c_cell.se_found){
            contam_rate_se.emplace(ci->first, c_cell.se);
        }
        else{
            contam_rate_se.emplace(ci->first, 0.0);
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
        if (contam_cell_prior > 0 && contam_cell_prior_se > 0){
            fprintf(stderr, "Shrunken contamination rate estimates:\n");
        }
        else{
            fprintf(stderr, "Contamination rate estimates:\n");
        }
        fprintf(stderr, "  Mean: %f Std dev: %f\n", mu_var.first, sqrt(mu_var.second));
        contam_cell_prior = mu_var.first;
        contam_cell_prior_se = sqrt(mu_var.second);
    }
}

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
    c_global.add_data("n", n);
    c_global.add_data("k", k);
    c_global.add_data("p_e", p_e);
    c_global.add_data("p_c", p_c);
    c_global.constrain_01();
    c_global.set_maxiter(-1);
    contam_cell_prior = c_global.solve(0,1);
    //contam_cell_prior_se = c_global.se;
    contam_cell_prior_se = -1;
    fprintf(stderr, "Global contamination rate estimate: %f\n", contam_cell_prior);
}

/**
 * After estimating contamination rate per cell, finds the likeliest ambient RNA
 * profile (set of p_c parameters for each type of SNP) based on the current
 * contamination rate estimates.
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
            
            if (omit_hets && ei2->first.second == 1){
                continue;
            }

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
 *
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

        if (true){
        //if (contam_rate.count(a->first) > 0){
            
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
                
                if (ref + alt == 0){
                    continue;
                }

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
*/

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
                    expected = adjust_p_err(expfracs[x->first][samp], e_r, e_a);
                }
                else{
                    if (x->first.first == samp){
                        expected = adjust_p_err((double)x->first.second/2.0, e_r, e_a);
                    }
                    else if (y->first.first == samp){
                        expected = adjust_p_err((double)y->first.second/2.0, e_r, e_a);
                    }
                    else{
                        expected = adjust_p_err(0.5*expfracs[x->first][samp] + 
                            0.5*expfracs[y->first][samp], e_r, e_a);
                    }
                }

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
    if (contam_prof.size() > 0){
        vector<double> startprops;
        for (int i = 0; i < idx2samp.size(); ++i){
            startprops.push_back(contam_prof[i]);
        }
        if (inter_species){
            startprops.push_back(contam_prof[-1]);
        }
        mix.add_mixcomp_fracs(startprops);
    }
    mix.solve(); 
    //mix.explore_starting_mixcomps();
    
    contam_prof.clear();
    for (int i = 0; i < idx2samp.size(); ++i){
        int samp = idx2samp[i];
        contam_prof.insert(make_pair(samp, mix.results_mixcomp[i]));
    }
    if (inter_species){
        contam_prof.insert(make_pair(-1, mix.results_mixcomp[mix.results_mixcomp.size()-1]));
    }
    /*
    for (int i = 0; i < mix.results_mixcomp.size(); ++i){
        int key = i;
        if (i >= n_samples){
            key = -1;
        }
        contam_prof.insert(make_pair(i, mix.results_mixcomp[i]));
    }
    */
    return mix.log_likelihood;

    /*
    //mix.add_mixcomp(mixfrac_exp);
    
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
    
    // Do another trial where starting proportions are based on 1-(mean contamination rate) per ID
    // Intuition here: if some IDs tend to have much lower inferred contam rates than others, that
    // probably means those IDs more closely resemble the contamination itself.
    map<int, double> id_contam_mean;
    map<int, double> id_contam_count;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (a->second >= n_samples){
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            if (id_contam_count.count(combo.first) == 0){
                id_contam_mean.insert(make_pair(combo.first, 0.0));
                id_contam_count.insert(make_pair(combo.first, 0.0));
            }
            if (id_contam_count.count(combo.second) == 0){
                id_contam_mean.insert(make_pair(combo.second, 0.0));
                id_contam_count.insert(make_pair(combo.second, 0.0));
            }
            id_contam_mean[combo.first] += contam_rate[a->first];
            id_contam_count[combo.first]++;
            id_contam_mean[combo.second] += contam_rate[a->first];
            id_contam_count[combo.second]++;
        } 
        else{
            if (id_contam_count.count(a->second) == 0){
                id_contam_mean.insert(make_pair(a->second, 0.0));
                id_contam_count.insert(make_pair(a->second, 0.0));
            }
            id_contam_mean[a->second] += contam_rate[a->first];
            id_contam_count[a->second]++;
        }
    }
    vector<double> starting_props1;
    double starting_props1_tot = 0.0;
    double mean_max = 0.0;
    for (map<int, double>::iterator m = id_contam_mean.begin(); m != id_contam_mean.end(); ++m){
        m->second /= id_contam_count[m->first];
        if (m->second > mean_max){
            mean_max = m->second;
        }
    }
    double pseudocount = 0.01;

    for (map<int, double>::iterator m = id_contam_mean.begin(); m != id_contam_mean.end(); ++m){
        double val = (mean_max - m->second)/mean_max + pseudocount;
        starting_props1.push_back(val);
        starting_props1_tot += val;
    }
    for (int i = 0; i < starting_props1.size(); ++i){
        starting_props1[i] /= starting_props1_tot;
    }
    mix.add_mixcomp_fracs(starting_props1);
    mix.solve();
    solution_ll.push_back(mix.log_likelihood);
    solution_cp.push_back(mix.results_mixcomp);

    if (!inter_species){
        map<int, int> assncount;
        int assntot;
        for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
            a != assn.end(); ++a){
            if (a->second >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                if (assncount.count(comb.first) == 0){
                    assncount.insert(make_pair(comb.first, 0));
                }
                if (assncount.count(comb.second) == 0){
                    assncount.insert(make_pair(comb.second, 0));
                }
                assncount[comb.first]++;
                assncount[comb.second]++;
                assntot += 2;
            }
            else{
                if (assncount.count(a->second) == 0){
                    assncount.insert(make_pair(a->second, 0));
                }
                assncount[a->second] += 2;
                assntot += 2;
            }
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
    */
}

/**
 * Compute log likelihood of data set with current parameters
 */
double contamFinder::compute_ll(){
    // Compute log likelihood of data set
    double loglik = 0.0;
    
    /* 
    pair<int, int> nullkey = make_pair(-1, -1);
    
    double c = contam_cell_prior;

    // Can't use pre-compiled counts since we may have changed assignments
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end();
        ++a){

        bool is_combo = false;
        pair<int, int> combo;
        if (a->second >= n_samples){
            is_combo = true;
            combo = idx_to_hap_comb(a->second, n_samples);
        }

        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator x =
            indv_allelecounts[a->first].begin(); x != indv_allelecounts[a->first].end(); ++x){
            
            if ((!is_combo && x->first.first == a->second) || (is_combo && combo.first == x->first.first)){
                
                if (!is_combo){
                    double ref = x->second[nullkey].first;
                    double alt = x->second[nullkey].second;
                    double exp = adjust_p_err((double)x->first.second/2.0, e_r, e_a);
                    double p = (1.0-c)*exp + c*amb_mu[x->first][nullkey];
                    loglik += dbinom(ref+alt, alt, p);
                }   
                else{
                    for (int nalt2 = 0; nalt2 <= 2; ++nalt2){
                        pair<int, int> key2 = make_pair(combo.second, nalt2);
                        double ref = x->second[key2].first;
                        double alt = x->second[key2].second;
                        double exp = adjust_p_err((double)(x->first.second + nalt2)/4.0, e_r, e_a);
                        double p = (1.0-c)*exp + c*amb_mu[x->first][key2];
                        loglik += dbinom(ref+alt, alt, p);
                    }
                }
            }
        }
    }
    return loglik;
    */

    for (int i = 0; i < n_all.size(); ++i){
        if (contam_rate.count(idx_to_cell[i]) > 0){
            double c = contam_rate[idx_to_cell[i]];
            //c = contam_cell_prior;
            double p_e = adjust_p_err(p_e_all[i], e_r, e_a);
            //double p_e = p_e_all[i];
            double p_c = amb_mu[type1_all[i]][type2_all[i]];
            double binom_p = (1.0-c)*p_e + c*p_c;
            loglik += logbinom(n_all[i], k_all[i], binom_p);
        }
    }
    return loglik;
}

double ll_test_aux(double p_c, const map<string, double>& params_d, const map<string, int>& params_i){
    double n = params_d.at("n");
    double k = params_d.at("k");
    double c = params_d.at("c");
    double p_0 = params_d.at("p");
    double binom_p = (1.0 - c)*p_0 + c*p_c;
    return logbinom(n, k, binom_p);
}

double dll_test_aux(double p_c, const map<string, double>& params_d, const map<string, int>& params_i){
    double n = params_d.at("n");
    double k = params_d.at("k");
    double c = params_d.at("c");
    double p_0 = params_d.at("p");
    double binom_p = (1.0 - c)*p_0 + c*p_c;
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    return dy_dp * c;
}


double contamFinder::test_aux(double c){
 
    double llsum = 0.0;

    for (map<pair<int, int>, map<pair<int, int>, vector<int> > >::iterator ei1 = 
        expfrac_to_idx.begin(); ei1 != expfrac_to_idx.end(); ++ei1){
        for (map<pair<int, int>, vector<int> >::iterator ei2 = ei1->second.begin(); ei2 != 
            ei1->second.end(); ++ei2){
            
            // Solve for p_c given current c
            vector<double> n;
            vector<double> k;
            vector<double> p_e;
            vector<double> weights;
            
            for (vector<int>::iterator i = ei2->second.begin(); i != ei2->second.end(); ++i){
                if (idx_to_cell.count(*i) > 0){ 
                    
                    n.push_back(n_all[*i]);
                    k.push_back(k_all[*i]);
                    p_e.push_back(adjust_p_err(p_e_all[*i], e_r, e_a));
                            
                    double weight = 1.0;
                    if (weighted){ 
                        // Weight by fraction of LLR sum of all cells assigned to this ID
                        // assigned to this cell
                        weight = assn_llr[idx_to_cell[*i]] / id_llrsum[assn[idx_to_cell[*i]]];
                    }
                    weights.push_back(weight);
                }
            }        
            
            optimML::brent_solver solver(ll_test_aux, dll_test_aux);
            solver.add_weights(weights);
            solver.add_data_fixed("c", c);
            solver.add_data("n", n);
            solver.add_data("k", k);
            solver.add_data("p", p_e);
            solver.constrain_01();
            double res = solver.solve(0,1);
            
            amb_mu[ei1->first][ei2->first] = res;

            //if (ei2->first.first == -1){
            //    fprintf(stderr, "%d\t%d\t%d\t%d\t%f\n", ei1->first.first, ei1->first.second,
            //        ei2->first.first, ei2->first.second, res);
            //}
            //fprintf(stdout, "%f\t%f\n", c, solver.log_likelihood);
            
            //if (true){
            //    llsum += solver.log_likelihood; 
            //}

        }
    }
    
    return solve_params_init();

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end(); ++a){
        bool is_combo = false;
        pair<int, int> combo;
        if (a->second >= n_samples){
            is_combo = true;
            combo = idx_to_hap_comb(a->second, n_samples);
        }
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac1 = 
            indv_allelecounts[a->first].begin(); ac1 != indv_allelecounts[a->first].end(); ++ac1){
            for (map<pair<int, int>, pair<float, float> >::iterator ac2 = ac1->second.begin();
                ac2 != ac1->second.end(); ++ac2){
                
                if (!is_combo){
                    if (ac1->first.first == a->second && ac2->first.first == -1){
                        double expec = adjust_p_err((double)ac1->first.second / 2.0, e_r, e_a);
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        double p = (1-c)*expec + c*amb_mu[ac1->first][ac2->first];
                        llsum += logbinom(ref+alt, alt, p);
                    }
                }
                else{
                    if (ac1->first.first == combo.first && ac2->first.first == combo.second){
                        double expec = adjust_p_err((double)(ac1->first.second + ac2->first.second)/4.0,
                            e_r, e_a);
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        double p = (1-c)*expec + c*amb_mu[ac1->first][ac2->first];
                        llsum += logbinom(ref+alt, alt, p);
                    }
                }
            }
        }
    }

    return llsum;
}

void contamFinder::test(){
    double llprev = 0.0;
    for (double c = 0.01; c <= 0.99; c += 0.01){
        amb_mu.clear();
        double ll = test_aux(c);
        fprintf(stderr, "c = %f: %f\n", c, ll);
        if (llprev != 0 && (ll-llprev) < 1){
            //break;
        }
        llprev = ll;
    }
}

double ll_test_new(const vector<double>& params, 
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

void dll_test_new(const vector<double>& params, 
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
    double binom_p = (1.0 - c) * p_e + c * params[1];
    double dy_dp = (k-n*binom_p)/(binom_p - binom_p*binom_p);
    results[results.size()-1] += dy_dp * c;
    if (c_in_data){
        results[0] += dy_dp * (p_c - p_e);
    }
}

/**
 *
 * Returns log likelihood.
 * Edits init_c to updated value (if solve_for_c == true)
 */
double contamFinder::update_amb_prof_mixture(bool solve_for_c, double& init_c, bool use_global_c){
    vector<double> params;
    if (solve_for_c){
        params.push_back(init_c);
    }
    optimML::multivar_ml_solver solver(params, ll_test_new, dll_test_new);
    
    vector<vector<double> > mixfracs;
    vector<double> weights;
    vector<double> n;
    vector<double> k;
    vector<double> p_e;
    vector<double> c;
    
    // true means short-circuit expfrac calculation - if one of two individuals
    // for a SNP type matches conditional individual, treat as that individual
    bool ef_all_avg = true;

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
        /*
        vector<pair<float, float> > reads;
        vector<float> expectations;

        int nc = 3;
        int nt = 2;
        if (is_comb){
            nc = 5;
            nt = 4;
        }
        for (int i = 0; i < nc; ++i){
            reads.push_back(make_pair(0.0, 0.0));
            expectations.push_back(adjust_p_err((double)nc/(double)nt, e_r, e_a));
        }
        */
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator ac1 = 
            indv_allelecounts[a->first].begin(); ac1 != indv_allelecounts[a->first].end(); ++ac1){
            
            //if (true){
            //if (!is_comb || (ac1->first.first == comb.first || ac1->first.first == comb.second)){
            //if (!is_comb || ac1->first.first == comb.first){
            if ((!is_comb && ac1->first.first == a->second) || (is_comb && ac1->first.first == comb.first)){
                for (map<pair<int, int>, pair<float, float> >::iterator ac2 = ac1->second.begin();
                    ac2 != ac1->second.end(); ++ac2){
                    
                    //if (ac2->first.first == -1){
                    //if ((!is_comb && (ac1->first.first == a->second || ac2->first.first == a->second)) || 
                    //    (is_comb && ac1->first.first == comb.first && ac2->first.first == comb.second)){
                    if ((!is_comb && ac2->first.first == -1) || (is_comb && ac2->first.first == comb.second)){
                        
                        double expected;
                        if (!is_comb){
                            /*
                            if (ac1->first.first == a->second){
                                expected = adjust_p_err((double)ac1->first.second/2.0, e_r, e_a);
                            }
                            else{
                                expected = adjust_p_err((double)ac2->first.second/2.0, e_r, e_a);
                            }
                            */
                            expected = adjust_p_err((double)ac1->first.second / 2.0, e_r, e_a);
                        }
                        else{
                            expected = adjust_p_err((double)(ac1->first.second + ac2->first.second)/4.0, e_r, e_a);
                            /*
                            if (ac1->first.first == comb.first){
                                expected = 0.5*(double)ac1->first.second/2.0 + 0.5*expfracs[ac1->first][comb.second];
                            }
                            else{
                                expected = 0.5*(double)ac1->first.second/2.0 + 0.5*expfracs[ac1->first][comb.first];
                            }
                            expected = adjust_p_err(expected, e_r, e_a);
                            */
                        }
                        
                        double ref = ac2->second.first;
                        double alt = ac2->second.second;
                        
                        /*   
                        int n_alt_exp;
                        if (!is_comb){
                            n_alt_exp = ac1->first.second;
                        }
                        else{
                            n_alt_exp = ac1->first.second + ac2->first.second;
                        }
                        reads[n_alt_exp].first += ref;
                        reads[n_alt_exp].second += alt;
                        */

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
                                    mixfrac_row.push_back(adjust_p_err(ac1->first.second / 2.0, e_r, e_a));
                                }
                                else if (ef_all_avg && ac2->first.first == samp){
                                    mixfrac_row.push_back(adjust_p_err(ac2->first.second / 2.0, e_r, e_a));
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
    
    // Starting proportions should be those previously set
    vector<double> startprops;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        if (cp->first != -1){
            startprops.push_back(cp->second);
        }
    }
    if (contam_prof.count(-1) > 0){
        startprops.push_back(contam_prof[-1]);
    }

    solver.add_mixcomp(mixfracs);
    solver.add_mixcomp_fracs(startprops);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_e", p_e);
    solver.add_weights(weights);
    if (solve_for_c){
        solver.constrain_01(0);
    }
    
    vector<double> ll_trials;
    vector<vector<double> > cp_trials;
    vector<double> c_trials;

    //solver.solve();
    //solver.explore_starting_mixcomps();

    if (solve_for_c){
        // First time. Try a few starting conditions and get the maximum LL.
        /*
        double llmax = solver.log_likelihood;
        int idxmax = 0;
        cp_trials.push_back(solver.results_mixcomp);
        c_trials.push_back(solver.results[0]);

        // Try an even mixture
        vector<double> mceven;
        for (int i = 0; i < n_samples; ++i){
            mceven.push_back(1.0/(double)n_samples);
        }
        solver.add_mixcomp_fracs(mceven);
        solver.solve();
        if (solver.log_likelihood > llmax){
            llmax = solver.log_likelihood;
            idxmax = cp_trials.size();
        }
        cp_trials.push_back(solver.results_mixcomp);
        c_trials.push_back(solver.results[0]);
        */
        // Try some random trials
        /*
        for (int x = 0; x < n_mixprop_trials; ++x){
            solver.randomize_mixcomps();
            solver.solve();
            if (solver.log_likelihood > llmax){
                llmax = solver.log_likelihood;
                idxmax = cp_trials.size();
            }
            cp_trials.push_back(solver.results_mixcomp);
            c_trials.push_back(solver.results[0]);
        }
        */
        
        //solver.explore_starting_mixcomps();
        //vector<double> test{ 0.6382283502, 0.0661063368, 0.0009443762, 0.0009443762, 0.0259703466, 0.2668618378, 0.0009443762 };
        //solver.add_mixcomp_fracs(test);
        //solver.solve();
        
              
        solver.solve();
        
             
        vector<double> lls;
        vector<vector<double> > mcs;
        vector<double> cs;
        int maxidx = 0;
        double maxll = solver.log_likelihood;
        lls.push_back(solver.log_likelihood);
        mcs.push_back(solver.results_mixcomp);
        cs.push_back(solver.results[0]);
    
        int ntrials = contam_prof.size() * n_mixprop_trials;
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
        

        //fprintf(stderr, "LL %f c = %f\n", solver.log_likelihood, solver.results[0]);
        contam_prof.clear();
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            //fprintf(stderr, "  %d) %f\n", samp, solver.results_mixcomp[i]);
            contam_prof.insert(make_pair(samp, solver.results_mixcomp[i]));
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, solver.results_mixcomp[solver.results_mixcomp.size()-1]));
        }
        /*
        fprintf(stderr, "LL %f c = %f\n", llmax, c_trials[idxmax]);
        contam_prof.clear();
        for (int i = 0; i < cp_trials[idxmax].size(); ++i){
            int key = i;
            if (i >= n_samples){
                key = -1;
            }
            fprintf(stderr, "  %d) %f\n", key, cp_trials[idxmax][i]);
            contam_prof.insert(make_pair(key, cp_trials[idxmax][i]));
        }
        exit(0);
        */
    }
    else{
        solver.solve();
        contam_prof.clear();
        for (int i = 0; i < idx2samp.size(); ++i){
            int samp = idx2samp[i];
            contam_prof.insert(make_pair(samp, solver.results_mixcomp[i]));
        }
        if (inter_species){
            contam_prof.insert(make_pair(-1, solver.results_mixcomp[solver.results_mixcomp.size()-1]));
        }
    }
    
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
            for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
                if (cp->first == -1){
                    val += cp->second * 0.0;
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
                    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
                        if (cp->first == -1){
                            val += cp->second * adjust_p_err(0.0, e_r, e_a);
                        }
                        else{
                            if (ef_all_avg && cp->first == key.first){
                                val += cp->second * adjust_p_err((double)key.second / 2.0, e_r, e_a);
                            }
                            else if (ef_all_avg && cp->first == key2.first){
                                val += cp->second * adjust_p_err((double)key2.second / 2.0, e_r, e_a);
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

void contamFinder::test_new(double init_c){
    
    if (!contam_prof_initialized){
        // Init contam prof
        contam_prof.clear();
        
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
        
        /*
        if (!inter_species){
            double tot = 0.0;
            map<int, double> props;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end(); ++a){
                if (a->second >= n_samples){
                    pair<int, int> comb = idx_to_hap_comb(a->second, n_samples);
                    if (props.count(comb.first) == 0){
                        props.insert(make_pair(comb.first, 0));
                    }
                    if (props.count(comb.second) == 0){
                        props.insert(make_pair(comb.second, 0));
                    }
                    props[comb.first]++;
                    props[comb.second]++;
                    tot += 2;
                }
                else{
                    if (props.count(a->second) == 0){
                        props.insert(make_pair(a->second, 0));
                    }
                    props[a->second] += 2;
                    tot += 2;
                }
            }
            for (int i = 0; i < idx2samp.size(); ++i){
                int samp = idx2samp[i];
                double min = 0.001;
                if (props.count(samp) == 0){
                    tot += min;
                    props.insert(make_pair(samp, min));
                }    
            }
            for (int i = 0; i < idx2samp.size(); ++i){
                int samp = idx2samp[i];
                contam_prof.insert(make_pair(samp, props[samp] / tot));
            }
        }
        else{
            for (int i = 0; i < idx2samp.size(); ++i){
                int samp = idx2samp[i];
                if (!inter_species){
                    contam_prof.insert(make_pair(samp, 1.0/(double)idx2samp.size()));
                }
                else{
                    contam_prof.insert(make_pair(samp, 1.0/(double)(idx2samp.size()+1)));
                }
            }
            if (inter_species){
                contam_prof.insert(make_pair(-1, 1.0/(double)(idx2samp.size()+1)));
            }
        }
        */
    }
    update_amb_prof_mixture(true, init_c);
}

double contamFinder::est_min_c(){
    map<pair<int, int>, map<pair<int, int>, double> > diffs;
    map<pair<int, int>, map<pair<int, int>, double> > counts;
    
    omit_hets = false;
    bool omit_hets_init = false;

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
                    if (diffs.count(key1) == 0){
                        map<pair<int, int>, double> m;
                        diffs.insert(make_pair(key1, m));
                        counts.insert(make_pair(key1, m));
                    }
                    if (diffs[key1].count(key2) == 0){
                        diffs[key1].insert(make_pair(key2, 0));
                        counts[key1].insert(make_pair(key2, 0));
                    }
                    if (!weighted){
                        diffs[key1][key2] += f;    
                        counts[key1][key2]++;
                    }
                    else{
                        diffs[key1][key2] += f * assn_llr[a->first];
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
                if (diffs.count(key) == 0){
                    map<pair<int, int>, double> m;
                    diffs.insert(make_pair(key, m));
                    counts.insert(make_pair(key, m));
                }
                if (diffs[key].count(nullkey) == 0){
                    diffs[key].insert(make_pair(nullkey, 0));
                    counts[key].insert(make_pair(nullkey, 0));
                }
                if (!weighted){
                    diffs[key][nullkey] += f;
                    counts[key][nullkey]++;
                }
                else{
                    diffs[key][nullkey] += f*assn_llr[a->first];
                    counts[key][nullkey] += assn_llr[a->first];
                }
            }
        }
    }

    double minc_mean = 0.0;
    double minc_weightsum = 0.0;
    for (map<pair<int, int>, map<pair<int, int>, double> >::iterator d = diffs.begin();
        d != diffs.end(); ++d){
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
                    // treat p_c = 0
                    minc = (inf - expec)/(1.0 - expec);
                }   
                else{
                    // treat p_c = 1
                    minc = (inf - expec)/(0.0 - expec);
                }
                //fprintf(stdout, "%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", d->first.first, d->first.second,
                //    d2->first.first, d2->first.second, d2->second / counts[d->first][d2->first], expec, 
                //    minc, counts[d->first][d2->first]);
                if (!omit_hets_init || (d->first.second != 1 && d2->first.second != 1)){
                    minc_mean += minc * counts[d->first][d2->first];
                    minc_weightsum += counts[d->first][d2->first];
                }
            }
        }
    }
    minc_mean /= minc_weightsum;
    return minc_mean;
}
/*
double contamFinder::ll_cell_diff(unsigned long cell, int id1, int id2){
    double c = contam_cell_prior;
    bool is_combo1 = false;
    bool is_combo2 = false;
    pair<int, int> combo1;
    pair<int, int> combo2;
    int shared = -1;
    if (id1 >= n_samples){
        is_combo1 = true;
        combo1 = idx_to_hap_comb(id1, n_samples);
    }
    if (id2 >= n_samples){
        is_combo2 = true;
        combo2 = idx_to_hap_comb(id2, n_samples);
    }
    if (is_combo1 && !is_combo2){
        if (combo1.first == id2){
            shared = combo1.first;
        }
        else if (combo1.second == id2){
            shared = combo2.first;
        }
    }
    else if (is_combo2 && !is_combo1){
        if (combo2.first == id1){
            shared = combo2.first;
        }
        else if (combo2.second == id1){
            shared = combo2.second;
        }
    }
    // Try to keep number of function evaluations consistent across the two IDs
    double ll1 = 0.0;
    double ll2 = 0.0;

    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator x = 
        indv_allelecounts[cell].begin(); x != indv_allelecounts[cell].end(); ++x){
        
        bool pass1 = false;
        if (is_combo1 && is_combo2){
            if (x->first.first == combo1.first || x->first.first == combo2.first){
                pass1 = true;
            }
        }
        else if (!is_combo1 && !is_combo2){
            if (x->first.first == id1 || x->first.first == id2){
                pass1 = true;
            }
        }
        else if (is_combo1 && !is_combo2){
            if (x->first.first == id2 || x->first.first == combo1.first){
                pass1 = true;
            }
        }
        else if (!is_combo1 && is_combo2){
            if (x->first.first == id1 || x->first.first == combo2.first){
                pass1 = true;
            }
        }
        if (pass1){
            for (map<pair<int, int>, pair<float, float> >::iterator y = x->second.begin(); y != x->second.end();
                ++y){
                
                if (!is_combo1 && !is_combo2){
                    if (y->first.first == -1){
                        double ref = y->second.first;
                        double alt = y->second.second;
                        double exp = adjust_p_err((double)x->first.second/2.0, e_r, e_a);
                        double p = (1-c)*exp + c*amb_mu[x->first][y->first];
                        if (x->first.first == id1){
                            ll1 += dbinom(ref+alt, alt, p);    
                        }
                        else{
                            ll2 += dbinom(ref+alt, alt, p);
                        }
                    }
                }
                else if (is_combo1 && is_combo2){
                    if (y->first.first == combo1.second || y->first.first == combo2.second){
                        double ref = y->second.first;
                        double alt = y->second.second;
                        double exp = adjust_p_err((double)(x->first.second + y->first.second)/4.0, e_r, e_a);
                        double p = (1-c)*exp + c*amb_mu[x->first][y->first];
                        if (y->first.first == combo1.second){
                            ll1 += dbinom(ref+alt, alt, p);
                        }
                        if (y->first.first == combo2.second){
                            ll2 += dbinom(ref+alt, alt, p);
                        }
                    }
                }
                else if (is_combo1 && !is_combo2){
                    if (combo1.first == id2){
                        if (y->first.first == combo1.second){
                            double ref = y->second.first;
                            double alt = y->second.second;
                            double exp1 = adjust_p_err((double)(x->first.second + y->first.second)/4.0, e_r, e_a);
                            double exp2 = adjust_p_err((double)x->first.second/2.0, e_r, e_a);
                            double p1 = (1.0-c)*exp1 + c*amb_mu[x->first][y->first];
                            double p2 = (1.0-c)*exp2 + c*amb_mu[x->first][y->first];
                            ll1 += dbinom(ref+alt, alt, p1);
                            ll2 += dbinom(ref+alt, alt, p2);
                        }
                    }
                    else if (combo1.second == id2){
                        if (y->first.first == combo1.second){
                            double ref = y->second.first;
                            double alt = y->second.second;
                            double exp1 = adjust_p_err((double)(x->first.second + y->first.second)/4.0, e_r, e_a);
                            double exp2 = adjust_p_err((double)y->first.second/2.0, e_r, e_a);
                            double p1 = (1.0-c)*exp1 + c*amb_mu[x->first][y->first];
                            double p2 = (1.0-c)*exp2 + c*amb_mu[x->first][y->first];
                            ll1 += dbinom(ref+alt, alt, p1);
                            ll2 += dbinom(ref+alt, alt, p2);
                        }
                    }
                    else{

                    }           
                }
                else if (is_combo2 && !is_combo1){

                }

            }
        }
    }
}
*/

/**
 * Returns whether anything was changed.
 */
bool contamFinder::reclassify_cells(){
    // Assume uniform global contam rate
    double c = contam_cell_prior;
    bool changed = false;
    
    set<unsigned long> cell_rm;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end(); ++a){
        
        // Get a table of log likelihood ratios between every possible
        // pair of identities
        map<int, map<int, double> > llrs;
        llr_table tab(n_samples);
        bool success = populate_llr_table(indv_allelecounts[a->first], llrs, tab, n_samples, allowed_ids, 
            allowed_ids2, 0.5, e_r, e_a, true, c, &amb_mu);
        if (success){
            int a_new;
            double llr_new; 
            tab.get_max(a_new, llr_new);
            if (a_new != -1 && llr_new > 0){
                if (a_new != a->second){
                    changed = true;
                    a->second = a_new;
                    assn_llr[a->first] = llr_new;
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

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != assn.end(); ++a){
        if (a->second >= n_samples){
            // Check against component singlets
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            double ll1 = 0.0;
            double ll2 = 0.0;
            double lldoub = 0.0;
            for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator x = 
                indv_allelecounts[a->first].begin(); x != indv_allelecounts[a->first].end(); ++x){
                for (map<pair<int, int>, pair<float, float> >::iterator y = x->second.begin();
                    y != x->second.end(); ++y){
                    if (x->first.first == combo.first && y->first.first == combo.second){
                        double ref = y->second.first;
                        double alt = y->second.second;
                        double expec1 = (1.0 - c)*((double)x->first.second/2.0) + c*amb_mu[x->first][y->first];
                        double expec2 = (1.0 - c)*((double)y->first.second/2.0) + c*amb_mu[x->first][y->first];
                        double expecd = (1.0 - c)*((double)(x->first.second + y->first.second)/4.0) + 
                            c*amb_mu[x->first][y->first];
                        expec1 = adjust_p_err(expec1, e_r, e_a);
                        expec2 = adjust_p_err(expec2, e_r, e_a);
                        expecd = adjust_p_err(expecd, e_r, e_a);
                        ll1 += dbinom(ref+alt, alt, expec1);
                        ll2 += dbinom(ref+alt, alt, expec2);
                        lldoub += dbinom(ref+alt, alt, expecd);
                    }
                }
            }
            int newid = -1;
            double llr = 0.0;
            if (ll1 >= lldoub){
                newid = combo.first;
                llr = ll1 - lldoub;
            }
            else if (ll2 >= lldoub){
                newid = combo.second;
                llr = ll2 - lldoub;
            }
            if (newid != -1){
                changed = true;
                a->second = newid;
                assn_llr[a->first] = llr;
            }
        }
    }
    return changed;
}

void contamFinder::set_init_c(double c){
    c_init = c;
}

/**
 * Find maximum likelihood estimates of all parameters of interest.
 * Return overall log likelihood.
 */
double contamFinder::fit(){
    //fprintf(stderr, "Computing expected allele matching rates...\n");
    //compute_expected_fracs_all_id();
    if (c_init <= 0){
        c_init = this->est_min_c();
    }
    fprintf(stderr, "Initial contamination rate = %f\n", c_init);
    test_new(c_init);

    //fprintf(stderr, "Get initial parameter estimates...\n"); 
    // Get initial estimates for global contamination rate and all
    // p_c parameters
    //double llinit = solve_params_init();
    //fprintf(stderr, "done\n");
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
        this->est_contam_cells_global();
        
        double dummy = -1;
        double loglik = this->update_amb_prof_mixture(false, dummy, true);

        // Get best estimate of contamination profile, given current
        // estimates of contamination per cell
        //double loglik = this->update_ambient_profile(true);
        //loglik = this->update_amb_prof_mixture(false, dummy);

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
    
    //this->update_ambient_profile(true);
    
    // Once to set prior params
    this->est_contam_cells();
    // Again to use empirical Bayes shrinkage
    this->est_contam_cells();

    // Allow ambient prof to update without considering individuals of origin
    this->update_ambient_profile(false);
    
    // Update mixture components using per-cell contam estimates
    double dummy;
    this->update_amb_prof_mixture(false, dummy, false);
    
    //est_contam_cells();
    fprintf(stderr, "Reclassifying cells...\n");
    bool reclassified = this->reclassify_cells();
    

    /*
    if (reclassified){
        est_contam_cells();
        update_ambient_profile();
        est_contam_cells();
        update_amb_prof_mixture(false, dummy);
    }
    */
    /*
    fprintf(stderr, "before reclass %f\n", compute_ll());
    while (reclassified){

        clear_data();
        compile_data(assn, indv_allelecounts);
        double ll1 = update_ambient_profile();
        fprintf(stderr, "reclass update %f\n", ll1);
        est_contam_cells();
        
        fprintf(stderr, "reclass %f\n", compute_ll());
        
        reclassified = reclassify_cells();
    }
    */
    if (reclassified){
        // Allows log likelihood computation
        clear_data();
        compile_data(assn, indv_allelecounts);
    }
    /*
    if (reclassified){
       update_amb_prof_mixture(false, dummy);
    }
    */
    fprintf(stderr, "est error rates...\n");
    // Check how low the error rates have dropped after modeling contamination
    pair<double, double> err_final = this->est_error_rates(false);
    fprintf(stderr, "Residual error rates:\n");
    fprintf(stderr, "  Reference alleles: %f\n", err_final.first);
    fprintf(stderr, "  Alt alleles: %f\n", err_final.second);
    
    /*
    if (model_mixprop){
        fprintf(stderr, "Computing expected allele matching rates...\n");
        this->compute_expected_fracs_all_id();
        fprintf(stderr, "Modeling ambient RNA as mixture of individuals...\n");
        double mixll = this->model_as_mixture();
    }
    */

    return llprev;
}

