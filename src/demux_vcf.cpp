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
#include <float.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
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
#include "demux_vcf_llr.h"

using std::cout;
using std::endl;
using namespace std;



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
    double error_rate_alt,
    bool use_prior_weights,
    map<int, double>& prior_weights){

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
        
        bool success;
        if (use_prior_weights){
            success = populate_llr_table(x->second, llrs, tab, samples.size(), allowed_assignments,
                allowed_assignments2, doublet_rate, error_rate_ref, error_rate_alt, &prior_weights);
        }
        else{
            success = populate_llr_table(x->second, llrs, tab, samples.size(), allowed_assignments, 
                allowed_assignments2,doublet_rate, error_rate_ref, error_rate_alt);
        }

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
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1){
        p = 1.0-DBL_MIN*1e6;
    }
    double ll = dbinom(n, k, p)/log2(exp(1));
    return ll;
}

double ll_err_persample(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    int ac1 = data_i.at("ac1");
    int ac2 = data_i.at("ac2");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    
    double p;
    if (idx2 == -1){
        double p0 = (double)ac1/2.0;
        double e_r = params[idx1*2];
        double e_a = params[idx1*2 + 1];
        p = p0 - p0*e_a + (1.0 - p0)*e_r;
    }
    else{
        double p0a = (double)ac1/2.0;
        double p0b = (double)ac2/2.0;
        double e_ra = params[idx1*2];
        double e_aa = params[idx1*2+1];
        double e_rb = params[idx2*2];
        double e_ab = params[idx2*2+1];
        p = 0.5*(p0a - p0a*e_aa + (1.0-p0a)*e_ra) + 
            0.5*(p0b - p0b*e_ab + (1.0-p0b)*e_rb);
    }
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1){
        p = 1.0-DBL_MIN*1e6;
    }
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
    double e_r = params[0];
    double e_a = params[1];
    double p0 = data_d.at("exp");
    double p = p0 - p0*e_a + (1.0 - p0)*e_r;
    if (p <= 0){
        p = DBL_MIN*1e6;
    }
    else if (p >= 1.0){
        p = 1.0-DBL_MIN*1e6;
    }
    double ll = dbinom(n, k, p)/log2(exp(1));
    
    double dy_dp = (k-n*p)/(p-p*p);
    double dp_de_a = -p0;
    double dp_de_r = 1.0 - p0;
    results[0] = dy_dp * dp_de_r;
    results[1] = dy_dp * dp_de_a;
}

void dll_err_persample(const vector<double>& params, const map<string, double>& data_d, 
    const map<string, int>& data_i, vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    int ac1 = data_i.at("ac1");
    int ac2 = data_i.at("ac2");
    
    double p;
    if (idx2 == -1){
        double p0 = (double)ac1/2.0;
        double e_r = params[idx1*2];
        double e_a = params[idx1*2+1];
        p = p0 - p0*e_a + (1.0-p0)*e_r;
        
        if (p <= 0){
            p = DBL_MIN*1e6;
        }
        else if (p >= 1.0){
            p = 1.0-DBL_MIN*1e6;
        }
        
        double dy_dp = (k-n*p)/(p-p*p);
        double dp_de_a = -p0;
        double dp_de_r = 1.0 - p0;
        results[idx1*2] += dy_dp * dp_de_r;
        results[idx1*2+1] += dy_dp * dp_de_a;
    }
    else{
        double p0a = (double)ac1/2.0;
        double p0b = (double)ac2/2.0;
        double e_ra = params[idx1*2];
        double e_aa = params[idx2*2+1];
        double e_rb = params[idx2*2];
        double e_ab = params[idx2*2+1];
        p = 0.5*(p0a - p0a*e_aa + (1.0-p0a)*e_ra) + 
            0.5*(p0b - p0b*e_ab + (1.0-p0b)*e_rb);
        if (p <= 0){
            p = DBL_MIN*1e6;
        }
        else if (p >= 1.0){
            p = 1.0-DBL_MIN*1e6;
        }
        
        double dy_dp = (k-n*p)/(p-p*p);
        results[idx1*2] += dy_dp * 0.5*(1.0 - p0a);
        results[idx2*2] += dy_dp * 0.5*(1.0 - p0b);
        results[idx1*2+1] += dy_dp * -0.5*p0a;
        results[idx2*2+1] += dy_dp * -0.5*p0b;
    }
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
    
    solver.set_silent(true);

    double sigma_curr = error_sigma;

    bool success = false;
    while (!success){
        success = true;
        // Sometimes, a math error will be encountered here. If that happens, set the 
        // standard deviation on the prior lower and lower until we can solve with
        // no error.
        try{ 
            solver.solve();
        } 
        catch (const int& err){
            if (err == optimML::OPTIMML_MATH_ERR){
                sigma_curr *= 0.5;
                fprintf(stderr, "Decreasing prior sd to %f...\n", sigma_curr);
                solver.set_prior_param(0, "sigma", sigma_curr);
                solver.set_prior_param(1, "sigma", sigma_curr);
                solver.set_param(0, error_ref);
                solver.set_param(1, error_alt);        
                success = false;
            }
            else{
                fprintf(stderr, "Unknown error encountered while inferring error rates\n");
                exit(1);
            }
        }
    }
    return make_pair(solver.results[0], solver.results[1]);
}

pair<double, double> infer_error_rates_persample(robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
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
    vector<int> ac1;
    vector<int> ac2;
    vector<double> expected;
    vector<double> weights_llr;
    vector<int> idx1;
    vector<int> idx2;

    vector<double> params;
    for (int i = 0; i < n_samples; ++i){
        params.push_back(error_ref);
        params.push_back(error_alt);
    }

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
                            ac1.push_back(x->first.second);
                            ac2.push_back(y->first.second);
                            n.push_back(y->second.first + y->second.second);
                            k.push_back(y->second.second);
                            weights_llr.push_back(weight);
                            idx1.push_back(combo.first);
                            idx2.push_back(combo.second);
                        }
                    }
                }
                else{
                    ac1.push_back(x->first.second);
                    ac2.push_back(-1);
                    n.push_back(x->second[nullkey].first + x->second[nullkey].second);
                    k.push_back(x->second[nullkey].second);
                    weights_llr.push_back(weight);
                    idx1.push_back(a->second);
                    idx2.push_back(-1);
                }
            }
        }
    }

    optimML::multivar_ml_solver solver(params, ll_err_persample, dll_err_persample);
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("ac1", ac1);
    solver.add_data("ac2", ac2);
    solver.add_data("idx1", idx1);
    solver.add_data("idx2", idx2);
    solver.add_weights(weights_llr);
    for (int x = 0; x < params.size(); ++x){
        solver.constrain_01(x);
        solver.add_normal_prior(x, params[x], error_sigma, 0.0, 1.0);
    }
    solver.solve();
    for (int i = 0; i < n_samples; ++i){
        fprintf(stderr, "%s) %f %f\n", samples[i].c_str(), solver.results[i*2], solver.results[i*2+1]);
    }
    exit(0);
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
    /*
    fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    */
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
       //{"index_jump", no_argument, 0, 'j'},
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
    while((ch = getopt_long(argc, argv, "b:v:o:B:i:I:q:D:n:e:E:p:s:fFCSUh", long_options, &option_index )) != -1){
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
            //case 'j':
            //    stream = false;
            //    break;
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
    else if (is_dir(output_prefix)){
        fprintf(stderr, "ERROR: output_prefix %s is a directory, but should be a file \
name prefix.\n", output_prefix.c_str());
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
    bool samples_from_vcf = false;

    // Store data about SNPs
    map<int, var> snpdat;

    int nsnps = 0;
    
    map<pair<int, int>, map<int, float> > conditional_match_fracs;
    map<pair<int, int>, map<int, float> > conditional_match_tots;

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
            read_vcf_samples(vcf_file, samples);         
        }
    }
    else{
        if (vcf_file == ""){
            fprintf(stderr, "ERROR: vcf file is required\n");
            exit(1);
        }
        
        read_vcf_samples(vcf_file, samples);
        
        if (dump_conditional){
            // All we are doing here is computing & writing out conditional match fracs.
            // Get a list of all chromosomes
            bcf_srs_t* sr = bcf_sr_init();    
            bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);

            for (int i = 0; i < bcf_header->n[BCF_DT_CTG]; ++i){
                // Read SNP data from this chromosome
                string hdr_chrom = bcf_hdr_id2name(bcf_header, i);
                map<int, var> snpstmp;
                
                read_vcf_chrom(vcf_file, hdr_chrom, snpstmp, vq); 

                // Add data to conditional_match_fracs
                get_conditional_match_fracs_chrom(snpstmp, conditional_match_fracs, 
                    conditional_match_tots, samples.size());
            }
            bcf_sr_destroy(sr);

            // Normalize conditional_match_fracs
            conditional_match_fracs_normalize(conditional_match_fracs,
                conditional_match_tots, samples.size());

            // Write data.
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            dump_exp_fracs(outf, conditional_match_fracs);
            fclose(outf);
            return 0;
        }

        // If not "dump_conditional", then we will load variants from each chromsome
        // and adjust conditional_match_fracs as we read through the BAM, writing out
        // results at the end.
        
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
    map<int, robin_hood::unordered_map<unsigned long, 
        pair<float, float> > > varcounts_site;
    
    robin_hood::unordered_map<unsigned long, map<int, float> > fracs;
    
    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;
    
    // Do not allow index-jumping; it's incompatible with loading SNPs one chromosome at a time.
    stream = true;

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
        
        int nsnp_processed = 0;
        if (stream){

            // Read through entire BAM file and look for informative SNPs along the way
            // (default setting, appropriate for large numbers of SNPs).
            int curtid = -1;
            
            map<int, var>::iterator cursnp;
            while (reader.next()){
                if (reader.unmapped() || reader.secondary() || reader.qcfail() || reader.dup()){
                    continue;
                }
                if (curtid != reader.tid()){
                    // Started a new chromosome
                    if (curtid != -1){
                        while (cursnp != snpdat.end()){
                            if (varcounts_site.count(cursnp->first) > 0){        
                                dump_vcs_counts(varcounts_site[cursnp->first], 
                                    indv_allelecounts,
                                    cursnp->second, samples.size());
                                varcounts_site.erase(cursnp->first);
                            }
                            ++nsnp_processed;
                            ++cursnp;
                        }
                    }
                    snpdat.clear();
                    char* curchromptr = reader.ref_id();
                    if (curchromptr == NULL){
                        cursnp = snpdat.end();
                    }
                    else{
                        string curchrom = curchromptr;
                        read_vcf_chrom(vcf_file, curchrom, snpdat, vq); 
                        cursnp = snpdat.begin();
                    }
                    curtid = reader.tid();
                    if (!disable_conditional){
                        // Compute contribution of new chrom to conditional match fracs
                        get_conditional_match_fracs_chrom(snpdat, 
                            conditional_match_fracs, conditional_match_tots, samples.size()); 
                    }
                }
                // Advance to position within cur read
                while (cursnp != snpdat.end() && 
                    cursnp->first < reader.reference_start){
                    
                    if (varcounts_site.count(cursnp->first) > 0){
                        dump_vcs_counts(varcounts_site[cursnp->first], 
                            indv_allelecounts, 
                            cursnp->second, samples.size());
                        varcounts_site.erase(cursnp->first);
                    }
                    ++nsnp_processed;
                    snpdat.erase(cursnp++);
                    //++cursnp;
                }
                // Create a second iterator to look ahead for any additional SNPs 
                // within the current read
                map<int, var>::iterator cursnp2 = cursnp;
                while (cursnp2 != snpdat.end() && 
                    cursnp2->first >= reader.reference_start && 
                    cursnp2->first <= reader.reference_end){   
                    process_bam_record(reader, cursnp2->first, cursnp2->second,
                        varcounts_site, cell_barcode, cell_barcodes);
                    ++cursnp2;
                }
                if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                    fprintf(stderr, "Processed %d SNPs\r", nsnp_processed); 
                    last_print = nsnp_processed;
                }
            }
            
            // Handle any final SNPs.
            if (curtid != -1){
                while (cursnp != snpdat.end()){
                    if (varcounts_site.count(cursnp->first) > 0){
                        dump_vcs_counts(varcounts_site[cursnp->first], 
                            indv_allelecounts,
                            cursnp->second, samples.size());
                        varcounts_site.erase(cursnp->first);
                    }
                    ++nsnp_processed;
                    snpdat.erase(cursnp++);
                    //++cursnp;    
                }
            }
        }
        else{
            /*
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
            */
        }
        fprintf(stderr, "Processed %d SNPs\n", nsnp_processed);
        
        // Write the data just compiled to disk.
        string fname = output_prefix + ".counts";
        //FILE* outf = fopen(fname.c_str(), "w");   
        gzFile outf = gzopen(fname.c_str(), "w");
        fprintf(stderr, "Writing allele counts to disk...\n");
        dump_cellcounts(outf, indv_allelecounts, samples);
        fprintf(stderr, "Done\n");
        //fclose(outf);
        gzclose(outf);

        // Finalize conditional match fracs
        if (!disable_conditional){
            
            // Normalize by total number of sites in calculations
            conditional_match_fracs_normalize(conditional_match_fracs, 
                conditional_match_tots, samples.size());
            
            // Write to disk.
            string outname = output_prefix + ".condf";
            FILE* outf = fopen(outname.c_str(), "w");
            dump_exp_fracs(outf, conditional_match_fracs);
            fclose(outf);
        }
    }
        
    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    fprintf(stderr, "Finding likeliest identities of cells...\n");
    map<int, double> prior_weights;
    
    // Note: prior weights are a thing built in here in case we want to influence
    // proportions of different types of IDs based on expectation. 
    // 
    // The thought was that if proportions of certain IDs change drastically before
    // and after inferring error rates, then those might be more present in ambient
    // RNA and we could downweight their expected proportions in the pool.
    //
    // This doesn't work well as of now and isn't being used, but I'm keeping the
    // code around in case we want to use it in the future.

    // Get assignments of cell barcodes
    map<pair<int, int>, double> er_map;
    map<pair<int, int>, double> ea_map; 
    assign_ids(indv_allelecounts, samples, assn, assn_llr, 
        allowed_ids, allowed_ids2, doublet_rate, error_ref, error_alt,
        false, prior_weights);

    robin_hood::unordered_map<unsigned long, int> assncpy = assn;
    
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
        allowed_ids, allowed_ids2, doublet_rate, error_ref_posterior, error_alt_posterior,
        false, prior_weights);
    
    // This function can infer sample-specific error rates, which might be better able to 
    // deal with ambient RNA. One drawback is that we can't infer rates for individuals
    // not yet assigned any cells. Another is that this is a double-edged sword:
    // Good: if there's a ton of ambient RNA from one ID that makes another ID look less
    //   likely, then we account for that and maybe get the correct ID more often.
    // Bad: if we mistakenly ID a bunch of cells to the wrong ID in the first round (because
    //   of ambient RNA contamination), then we are assuming those IDs are often correct
    //   and perhaps making it more likely to choose it on the next round.
    // For now, stick to the original plan, which treats all individuals the same and relies
    //   on the genotype data.    
    //infer_error_rates_persample(indv_allelecounts, samples.size(), assn, assn_llr, error_ref,
    //    error_alt, error_sigma, samples);

    map<int, int> assncount1;
    map<int, int> assncount2;
    int a2tot = 0;
    int a1tot = 0;
    for (int i = 0; i < samples.size(); ++i){
        assncount2.insert(make_pair(i, 1));
        a2tot++;
        assncount1.insert(make_pair(i, 1));
        a1tot++;
        for (int j = i + 1; j < samples.size(); ++j){
            int k = hap_comb_to_idx(i, j, samples.size());
            assncount2.insert(make_pair(k, 1));
            a2tot++;
            assncount1.insert(make_pair(k, 1));
            a1tot++;
        }
    }
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assncpy.begin(); a != 
        assncpy.end(); ++a){
        assncount1[a->second]++;
        a1tot++;
    }
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin(); a != 
        assn.end(); ++a){
        assncount2[a->second]++;
        a2tot++;
    }
    /*
    for (map<int, int>::iterator a = assncount2.begin(); a != assncount2.end(); ++a){
        double p = pbinom(a2tot, assncount2[a->first], (double)assncount1[a->first]/(double)a1tot);
        if (p == 1.0){
            p = 1-1e-6;
        }
        else if (p == 0){
            p = 1e-6;
        }
        //p = log2(assncount2[a->first]) - log2(assncount1[a->first]);
        fprintf(stderr, "%s) %d %d | %d %d | %f\n", idx2name(a->first, samples).c_str(),
            assncount1[a->first], a1tot, assncount2[a->first], a2tot,
            log2(p));
            //1-pbinom(a1tot, assncount2[a->first], (double)assncount1[a->first]/(double)a1tot)); 
        prior_weights.insert(make_pair(a->first, p));
    }
    */

    bool do_again = true;
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
                    error_alt_posterior,
                    false, prior_weights); 
                do_again = false;
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
