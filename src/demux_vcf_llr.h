#ifndef _DEMUX_VCF_LLR_H
#define _DEMUX_VCF_LLR_H
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
#include <set>
#include <cstdlib>
#include <utility>
#include "common.h"

class llr_table{
    private:
        // Set of pairs rather than map allows duplicate LLRs
        std::map<double, std::vector<std::pair<short, short> > > lookup_llr;
        std::map<double, std::vector<std::pair<short, short> > >::iterator it;
        std::vector<double> maxllr;
        std::vector<double> minllr;

    public:
        std::vector<bool> included;
        int n_indvs;
        
        // Constructor
        llr_table(int x);

        // Destructor
        ~llr_table();
        
        // Print contents of table to stdout.
        void print(std::string& bc_str, std::vector<std::string>& samples);
        
        void print_ranges(std::string& barcode, std::vector<std::string>& samples);
        
        // Add a new item
        void insert(short i1, short i2, double llr);
        
        void disallow(short i);
        
        // Thin the table out, keeping the set number of indvs
        bool del(int n_keep);
        
        // Make a decision from the table.
        void get_max(int& best_idx, double& best_llr);
};

bool populate_llr_table(std::map<std::pair<int, int>,
            std::map<std::pair<int, int>, 
                std::pair<float, float> > >& counts,
    std::map<int, std::map<int, double> >& llrs,
    llr_table& tab,
    int n_samples,
    std::set<int>& allowed_assignments,
    std::set<int>& allowed_assignments2,
    double doublet_rate,
    double error_rate_ref,
    double error_rate_alt,
    bool incl_contam=false,
    double contam_rate=0.0,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >* amb_fracs=NULL,
    std::map<int, double>* prior_weights=NULL);

// Incorporate ref & alt mismatch rates into a probability
double adjust_p_err(double p, double e_r, double e_a);

#endif
