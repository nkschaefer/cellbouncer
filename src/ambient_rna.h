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
        
        robin_hood::unordered_map<unsigned long, 
                std::map<std::pair<int, int>, 
                    std::map<std::pair<int, int>, 
                        std::pair<float, float> > > > indv_allelecounts;
        std::set<int> allowed_ids; 
        std::set<int> allowed_ids2;

        // Map each ID to the sum of all LLRs of all cells assigned to it
        map<int, double> id_llrsum;
        map<int, double> id_llrsum2;
        map<int, double> id_count;
        double llrtot;
        
        // What is the expected doublet rate (for reclassifying cells)?
        double doublet_rate;

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
        
        // Prior distribution mean for contamination rate per cell
        double contam_cell_prior;
        // Prior distribution sigma for contamination rate per cell
        double contam_cell_prior_se;
       
        // Error rate for ref alleles (rate at which they are misread as alt)
        double e_r;
        // Error rate for alt alleles (rate at which they are misread as ref)
        double e_a;
        
        // Should we model contamination from other species? This assumes that
        // the reference genome is for the current species, that alt alleles are
        // generally derived, and that a foreign species will appear as an excess
        // of reference alleles
        bool inter_species;

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
        
        void clear_data();

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
        void est_contam_cells();
        void est_contam_cells_global();

        // Returns log likelihood
        double update_ambient_profile(bool global_c = false);
        std::pair<double, double> est_error_rates(bool init);
        
        double c_init_global;
        
        void compile_amb_prof_dat(bool solve_for_c, 
            bool use_global_c,
            std::vector<std::vector<double> >& mixfracs,
            std::vector<double>& weights,
            std::vector<double>& n,
            std::vector<double>& k,
            std::vector<double>& p_e,
            std::vector<double>& c);

        double update_amb_prof_mixture(bool est_c, double& global_c, bool use_global_c = false); 
        double est_min_c();
        bool reclassify_cells();
        double init_params(double& c);
        bool contam_prof_initialized;
        bool c_initialized;
        double c_init;

    public:
        
        void set_init_contam_prof(std::map<int, double>& cp);
        void set_init_c(double c);
        
        void set_doublet_rate(double d);

        // Copy of assignments & LLRs of assignments from demux_vcf. These
        // may change if cells are re-assigned after contamination
        // inference
        robin_hood::unordered_map<unsigned long, int> assn;
        robin_hood::unordered_map<unsigned long, double> assn_llr;
        
        // Mean alt allele matching frequencies in ambient RNA
        // First key: allele type (individual index, number alt alleles)
        //   Second key: allele type 2 (individual index, number alt alleles, or (-1, -1) for all
        //     Value: Expected rate of matching alt alleles in ambient RNA
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> > amb_mu;
        
        // Contamination rates & their standard errors per cell
        robin_hood::unordered_map<unsigned long, double> contam_rate;
        robin_hood::unordered_map<unsigned long, double> contam_rate_se;
        robin_hood::unordered_map<unsigned long, double> contam_rate_ll;

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
            int n_samples,
            std::set<int>& allowed_ids,
            std::set<int>& allowed_ids2);

        // Functions to update parameter values
        void set_error_rates(double e_r, double e_a);
        void model_other_species();
        void model_single_species();
        void set_mixprop_trials(int nt);
        void set_delta(double d);
        void set_maxiter(int i);
        void use_weights();
        void no_weights();
        
        // Get variance on contamination profile by bootstrapping it and fitting
        // a Dirichlet distribution
        void bootstrap_amb_prof(int n_boots, std::map<int, double>& contam_prof_cont);

        // Get log likelihood of data under current parameters
        double compute_ll();

        // Run everything
        // Returns log likelihood
        void fit();
};

#endif
