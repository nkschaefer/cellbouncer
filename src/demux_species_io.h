#ifndef _DEMUX_SPECIES_IO_H
#define _DEMUX_SPECIES_IO_H
#include <zlib.h>
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
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "demux_species_io.h"

bool isdigit(char c);
bool filename_base(std::string& fn, std::string& fntrunc);
void print_bc_species_counts(
    robin_hood::unordered_map<unsigned long, std::map<short, int> >& bc_species_counts,
    std::map<short, std::string>& idx2species,
    FILE* countsfile);
void load_from_files(std::string& countsfilename,
    std::string& speciesfilename,
    std::map<short, std::string>& idx2species,
    std::map<std::string, short>& species2idx,
    robin_hood::unordered_map<unsigned long, std::map<short, int> >& bc_species_counts);
void print_assignments(FILE* outf, 
    const std::string& libname,
    robin_hood::unordered_map<unsigned long, short>& bc2species,
    robin_hood::unordered_map<unsigned long, std::pair<unsigned int, unsigned int> >& bc2doublet,
    robin_hood::unordered_map<unsigned long, double>& bc2llr,
    std::map<short, std::string>& idx2species,
    bool use_filter,
    robin_hood::unordered_set<unsigned long>& filter);
void create_library_file(std::vector<std::string>& rna_r1files,
    std::vector<std::string>& atac_r1files,
    std::vector<std::string>& custom_r1files,
    std::vector<std::string>& custom_names,
    std::map<short, std::string>& idx2species,
    const std::string& outdir);



 
 
#endif
