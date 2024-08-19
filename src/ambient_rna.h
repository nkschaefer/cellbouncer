#ifndef _CELLID_AMBIENT_RNA_H
#define _CELLID_AMBIENT_RNA_H
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
#include <optimML/multivar.h>
#include "common.h"

// Class to store everything needed for characterizing ambient RNA & 
// estimating per-cell contamination. Cuts down on number of data
// structures to pass to functions.

class contamFinder{
    private:
        
        // Copies of external data structures needed by multiple things
        robin_hood::unordered_map<unsigned long, int> assn;
        robin_hood::unordered_map<unsigned long, double> assn_llr;
        robin_hood::unordered_map<unsigned long, 
                std::map<std::pair<int, int>, 
                    std::map<std::pair<int, int>, 
                        std::pair<float, float> > > > indv_allelecounts;
        std::set<int> allowed_ids; 
        // Map each ID to the sum of all LLRs of all cells assigned to it
        map<int, double> id_llrsum;

        // Vectors of each data value that we will use repeatedly in optimizations 
        std::vector<double> n_all;
        std::vector<double> k_all;
        std::vector<double> p_e_all;
        std::vector<std::pair<int, int> > type1_all;
        std::vector<std::pair<int, int> > type2_all;
        std::vector<double> weight_all;
    
        // Data structures to map different values of interest to indices in the above vectors
        std::map<std::pair<int, int>, std::map<pair<int, int>, std::vector<int> > > expfrac_to_idx;
        std::map<unsigned long, std::vector<int> > cell_to_idx;
        std::map<int, unsigned long> idx_to_cell;
        std::map<std::pair<int, int>, std::vector<int> > sitecomb_type_to_idx;
        
        // How many individuals were in the VCF + set under consideration?
        int n_samples;
        std::vector<int> idx2samp;

        // First key: allele type (individual index, number alt alleles)
        //  Second key: true identity
        //    Value: Expected rate of matching alt alleles at sites of this type
        std::map<std::pair<int, int>, std::map<int, float> > expfracs; 
        

        map<pair<int, int>, map<pair<int, int>, double> > amb_tot;
        
        // Store backup copy of above during each iteration, in case overall log likelihood 
        // decreases 
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> > amb_mu_prev;
        
        // Other backup copies kept for same reason
        robin_hood::unordered_map<unsigned long, double> contam_rate_prev;
        robin_hood::unordered_map<unsigned long, double> contam_rate_prev_se;
        
        // Prior distribution mean for contamination rate per cell
        double contam_cell_prior;
        // Prior distribution sigma for contamination rate per cell
        double contam_cell_prior_se;
       
        // Optionally limit the number of cells used in computing expected alt
        // allele fraction at each allele type given each possible identity
        // (very slow with the entire data set) 
        int max_cells_expfrac;

        // Error rate for ref alleles (rate at which they are misread as alt)
        double e_r;
        // Error rate for alt alleles (rate at which they are misread as ref)
        double e_a;
        
        // Should we model contamination from other species? This assumes that
        // the reference genome is for the current species, that alt alleles are
        // generally derived, and that a foreign species will appear as an excess
        // of reference alleles
        bool inter_species;

        // Should we take the last step of modeling ambient RNA as a mixture of 
        // known individuals?
        bool model_mixprop;
        
        // How many random trials should we do to re-adjust starting proportions
        // when modeling ambient RNA as a mixture of individuals?
        int n_mixprop_trials;

        // Stopping criteria for iteration
        double delta_thresh;
        int maxits;
        int nits;
        
        // Should we consider log likelihood ratio of assignments when inferring
        // contamination?
        bool weighted;

        // Wrangle data
        void compile_data(robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long, 
                std::map<std::pair<int, int>, 
                std::map<std::pair<int, int>, 
                std::pair<float, float> > > >& indv_allelecounts);
        
        void get_reads_expectations(int ident,
            std::map<std::pair<int, int>, 
                std::map<std::pair<int, int>, 
                std::pair<float, float> > >& allelecounts,
            std::vector<double>& n,
            std::vector<double>& k,
            std::vector<double>& p_e,
            std::vector<std::pair<int, int> >& type1,
            std::vector<std::pair<int, int> >& type2);
        
        // Calculate stuff
        double solve_params_init(); 
        void est_contam_cells();
        // Returns log likelihood
        double update_ambient_profile();
        std::pair<double, double> est_error_rates(bool init);
        void compute_expected_fracs_all_id();
        double model_as_mixture();
        
        void test();
        double test_aux(double c);
        bool omit_hets;
        double c_init_global;
        
        double update_amb_prof_mixture(bool est_c, double& global_c); 
        double est_min_c();
        void test_new(double c);

    public:
        
        // Mean alt allele matching frequencies in ambient RNA
        // First key: allele type (individual index, number alt alleles)
        //   Second key: allele type 2 (individual index, number alt alleles, or (-1, -1) for all
        //     Value: Expected rate of matching alt alleles in ambient RNA
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> > amb_mu;
        
        // Contamination rates & their standard errors per cell
        robin_hood::unordered_map<unsigned long, double> contam_rate;
        robin_hood::unordered_map<unsigned long, double> contam_rate_se;
        
        // Most likely mixture of genotypes in ambient RNA
        std::map<int, double> contam_prof;

        // Constructor
        contamFinder(robin_hood::unordered_map<unsigned long, 
                std::map<std::pair<int, int>, 
                    std::map<std::pair<int, int>, 
                        std::pair<float, float> > > >& indv_allelecounts,
            robin_hood::unordered_map<unsigned long, int>& assn,
            robin_hood::unordered_map<unsigned long, double>& assn_llr,
            std::map<std::pair<int, int>, std::map<int, float> >& exp_match_fracs,
            int n_samples);

        // Functions to update parameter values
        void set_error_rates(double e_r, double e_a);
        void model_other_species();
        void model_single_species();
        void set_mixprop_trials(int nt);
        void model_mixture();
        void skip_model_mixture();
        void set_delta(double d);
        void set_maxiter(int i);
        void max_cells_for_expfracs(int mc);
        void use_weights();
        void no_weights();

        // Get log likelihood of data under current parameters
        double compute_ll();

        // Run everything
        // Returns log likelihood
        double fit();

};

#endif
