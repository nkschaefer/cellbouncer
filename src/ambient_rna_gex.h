#ifndef _CELLBOUNCER_AMBIENT_RNA_GEX_H
#define _CELLBOUNCER_AMBIENT_RNA_GEX_H
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
#include <random>
#include <mixtureDist/functions.h>
#include <optimML/multivar.h>
#include "common.h"

// Class to store everything needed for characterizing ambient RNA & 
// estimating per-cell contamination. Cuts down on number of data
// structures to pass to functions.

class contam_profiler_gex{
    private:
        
        // Copies of external data structures needed by multiple things
        robin_hood::unordered_map<unsigned long, double>* contam_rate;
        map<int, double>* contam_prof;
        robin_hood::unordered_map<unsigned long, int>* assn;
        int n_samples; 
        robin_hood::unordered_map<unsigned long, int>* clusts;
        int n_clusts;
        
        robin_hood::unordered_map<unsigned long, map<int, long int> >* mtx;
        int n_features;
        
        int nthreads;

        robin_hood::unordered_map<unsigned long, map<int, long int> >* ase_mtx1;
        robin_hood::unordered_map<unsigned long, map<int, long int> >* ase_mtx2;
        
        robin_hood::unordered_map<unsigned long, double> cells_tot;

        bool doublets_as_indvs;
        bool profile_set;
        bool has_ase;
        bool round_counts_mtx;
        
        std::random_device rand_dev;
        std::mt19937 rand_gen;
        std::uniform_real_distribution<> rand_dist;

        bool infer_ase();
        
        bool skip_genes_given;
        std::set<int> skip_genelist;

    public:
        
        contam_profiler_gex(robin_hood::unordered_map<unsigned long, double>& contam_rate,
            map<int, double>& contam_prof,
            robin_hood::unordered_map<unsigned long, int>& assn,
            int n_samples,
            bool doublets_as_indvs = false);

        void set_mtx(robin_hood::unordered_map<unsigned long, map<int, long int> >& mtx,
            int n_features);

        void set_clusts(robin_hood::unordered_map<unsigned long, int>& clusts,
            int n_clusts);

        void set_threads(int n_threads);
        
        void skip_genes(std::set<int>& genelist);

        void set_ase(robin_hood::unordered_map<unsigned long, map<int, long int> >& ase1,
            robin_hood::unordered_map<unsigned long, map<int, long int> >& ase2);
        
        void round_counts();

        // Returns true = worked, or false = did not work
        bool get_profile();
        
        bool decontam();

        // Result data structures
        
        // Ambient RNA expression profile (multinomial parameters)
        std::vector<double> prof_ambient;
        // Endogenous RNA expression profile per cluster (multinomial parameters)
        std::vector<std::vector<double> > prof_clusts;
        
        // Matrix with ambient counts removed
        robin_hood::unordered_map<unsigned long, map<int, double> > mtx_decontam;
        
        // Fraction allele1, allele2, and neither for each gene for ASE data
        vector<vector<double> > ase_fracs;

        // ASE matrices with ambient counts removed
        robin_hood::unordered_map<unsigned long, map<int, double> > mtx_decontam_ase1;
        robin_hood::unordered_map<unsigned long, map<int, double> > mtx_decontam_ase2;
        

};

#endif
