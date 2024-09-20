#ifndef _CELLID_COMMON_H
#define _CELLID_COMMON_H
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
#include <map>
#include <set>
#include <cstdlib>
#include <utility>
#include <htswrapper/bc.h>
/**
 * Contains functions used by more than one program in this
 * repository.
 */

using std::cout;
using std::endl;
using namespace std;

// Print help screen information about the parameter for
// library names that will be appended to cell barcodes
// (used by multiple programs)
void print_libname_help();

// Modifies a string representation of a cell barcode to include unique
// library information before printing.
void mod_bc_libname(string& bc, const string& libname, bool cellranger,
   bool seurat, bool underscore);

// Parse a file mapping cell barcodes to identities
void parse_barcode_map(std::string& fn, 
    std::map<unsigned long, std::string>& bc2hap,
    std::set<std::string>& barcode_groups,
    double llr_cutoff,
    bool keep_doublets);

// Utility functions for translating back and forth between
// identity indices and names, using a standard way of 
// representing doublet identities as integers
short hap_comb_to_idx(short, short, short);
std::pair<short, short> idx_to_hap_comb(short, short);
std::string idx2name(int, std::vector<std::string>&);

// Represent/print a triangular row-major distance matrix
void init_distmat(vector<vector<float> >& distmat, int dim);
void print_distmat(vector<vector<float> >& distmat);
void print_distmat_square(vector<vector<float> >& distmat);

// Transform a log likelihood ratio table of possible cell identities
// into a single identity and its log likelihood ratio
int collapse_llrs(std::map<int, std::map<int, double> >& llrs, double& llr_final);

// Check for over/underrepresentation of counts of specific 
// doublet types in a data set
double doublet_chisq(std::map<int, int>& idcounts, int n_samples);

// Trim the path off of a file name
std::string filename_nopath(std::string& filename);

// Log PDF of binomial distribution
double logbinom(double n, double k, double p);

bool file_exists(std::string name);

// Approximate derivative using slope of neighboring points
void derivative(std::map<double, double>& hist, std::map<double, double>& result, int smooth);

// Find inflection point in a histogram
double find_knee(std::map<double, double>& hist, double min_frac_to_allow);

// Load a market exchange format (MEX) file
void parse_mex(const std::string& barcodesfile,
    const std::string& featuresfile,
    const std::string& matrixfile, 
    robin_hood::unordered_map<unsigned long, map<int, long int> >& counts,
    std::vector<std::string>& labels,
    const std::string& featuretype = "");

void fit_dirichlet(std::vector<double>& mle_fracs,
    std::vector<std::vector<double> >& dirichlet_bootstraps,
    std::vector<double>& conc_param_results);

#endif
