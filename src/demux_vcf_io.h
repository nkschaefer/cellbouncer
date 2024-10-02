#ifndef _CELLBOUNCER_DEMUX_VCF_IO_H
#define _CELLBOUNCER_DEMUX_VCF_IO_H
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
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"

void parse_idfile(std::string& idfile, 
    std::vector<std::string>& samples,
    std::set<int>& ids_allowed,
    std::set<int>& ids_allowed2,
    bool add_all_doublets);

void write_samples(std::string& filename, 
    std::vector<std::string>& samples);

void load_samples(std::string& filename,
    std::vector<std::string>& samples);

void write_allowed(std::string& filename,
    std::set<int>& allowed,
    std::vector<std::string>& samples);

void load_allowed(std::string& filename,
    std::set<int>& allowed,
    std::vector<std::string>& samples);

void load_counts_from_file(
    robin_hood::unordered_map<unsigned long, std::map<std::pair<int, int>, 
        std::map<std::pair<int, int>, std::pair<float, float> > > >& indv_allelecounts,
    std::vector<std::string>& indvs,   
    std::string& filename,
    std::set<int>& allowed_ids);

void dump_cellcounts(gzFile& out_cell,
    robin_hood::unordered_map<unsigned long, std::map<std::pair<int, int>, 
        std::map<std::pair<int, int>, std::pair<float, float> > > >& indv_allelecounts, 
    std::vector<std::string>& samples);

void load_exp_fracs(std::string& filename,   
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_frac);

void dump_exp_fracs(FILE* exp_frac_out,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_frac);

void dump_assignments(FILE* outf,
    robin_hood::unordered_map<unsigned long, int>& assn_final,
    robin_hood::unordered_map<unsigned long, double>& assn_final_llr,
    std::vector<std::string>& samples,
    std::string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore);

void load_assignments_from_file(std::string& filename,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    std::vector<std::string>& samples);

void write_summary(FILE* outf, 
    std::string& outpre,
    robin_hood::unordered_map<unsigned long, int>& assn,
    std::vector<std::string>& samples,
    double error_ref,
    double error_alt,
    double error_sigma,
    double ref_mm_rate_posterior,
    double alt_mm_rate_posterior,
    std::string& vcf_file,
    int vq_filter,
    double doublet_rate,
    std::map<int, double>& p_ncell,
    std::map<int, double>& p_llr);

void dump_contam_prof(FILE* outf,
    std::map<int, double>& contam_prof,
    std::map<int, double>& contam_prof_conc,
    std::vector<std::string>& samples);

void dump_contam_rates(FILE* outf,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    robin_hood::unordered_map<unsigned long, double>& contam_rate_se,
    std::vector<std::string>& samples,
    std::string& libname,
    bool cellranger,
    bool seurat,
    bool underscore);

void dump_amb_fracs(FILE* outf, 
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double> >& amb_mu);

#endif
