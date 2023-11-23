#ifndef _UNPOOL_COMMON_H
#define _UNPOOL_COMMON_H
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
/**
 * Contains functions used by more than one program in this
 * repository.
 */

using std::cout;
using std::endl;
using namespace std;

void parse_barcode_map(std::string& fn, 
    std::map<std::string, std::string>& bc2hap,
    std::set<std::string>& barcode_groups,
    double llr_cutoff,
    bool keep_doublets);

short hap_comb_to_idx(short, short, short);
std::pair<short, short> idx_to_hap_comb(short, short);
std::string idx2name(int, std::vector<std::string>&);

void init_distmat(vector<vector<float> >& distmat, int dim);
void print_distmat(vector<vector<float> >& distmat);
void print_distmat_square(vector<vector<float> >& distmat);

int collapse_llrs(std::map<int, std::map<int, double> >& llrs, double& llr_final);

#endif
