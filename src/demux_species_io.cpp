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
#include <regex>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_species_io.h"

using namespace std;

/**
 * Is the given character a digit between 0 and 9?
 */
bool isdigit(char ch){
    static char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    for (int i = 0; i < 10; ++i){
        if (ch == digits[i]){
            return true;
        }
    }
    return false;
}

/**
 * Truncate a file name to come up with a "base" file name that 
 * can be used to access/write other files in the same group.
 */
bool filename_base(string& filename, string& fnbase){
    
    // Stricter
    static const std::regex fq_regex1("^(.+)_S([0-9]+)_L([0-9]+)_R(1|2|3)(_[0-9]{3})?(\\.[0-9]+)?\\.(fastq|fq)(\\.gz)?$");
    // More permissive
    static const std::regex fq_regex2("^(.+)_R(1|2|3)(_[0-9]{3})?(\\.[0-9]+)?\\.(fastq|fq)(\\.gz)?$");
    smatch matches;
    if (regex_match(filename, matches, fq_regex1)){
        //fnbase = matches[1].str() + "_S" + matches[2].str() + "_L" + matches[3].str();
        fnbase = matches[1].str();
        return true;
    }
    else{
        // Try more permissive matching
        if (regex_match(filename, matches, fq_regex2)){
            fnbase = matches[1].str();
            return true;
        }
    }
    return false;

    /*
    static const std::regex fq_regex("^(.+)(_S[0-9]+)?(_L[0-9]+)?_R(1|2|3)(_[0-9]{3})\.(fastq|fq)(\.gz)$");
    smatch matches;

    if (regex_match(filename, matches, fq_regex)){
        string grp = matches[1].str();
         
        str2bc(matches[3].str().c_str(), as_bitset);
        return as_bitset.to_ulong();
    }
    
    stringstream splitter(filename);
    string chunk;
    vector<string> chunks;
    bool sfound = false;
    bool lfound = false;
    bool rfound = false;
    bool numfound = false;
    
    fnbase = "";    
    while (getline(splitter, chunk, '_')){
        chunks.push_back(chunk);
    }
    
    for (vector<string>::reverse_iterator chunk = chunks.rbegin(); chunk != chunks.rend(); ++chunk){
        if (!numfound){
            string suffix = chunk->substr(chunk->length()-8-1, 9);
            string prefix = chunk->substr(0, chunk->length()-8-1);
            bool pass = true;
            if (suffix == ".fastq.gz"){
                for (int i = 0; i < prefix.length(); ++i){
                    if (!isdigit(prefix[i])){
                        pass = false;
                        break;
                    }
                }
            }
            else{
                pass = false;
            }
            if (pass){
                numfound = true;
            }
            else{
                return false;
            }
        }
        else{
            if (!rfound){
                if ((*chunk)[0] == 'R' && isdigit((*chunk)[1])){
                    rfound = true;
                }
                else{
                    return false;
                }
            }
            else{
                if (!lfound){ 
                    if ((*chunk)[0] == 'L'){
                        bool pass = true;
                        for (int i = 1; i < chunk->length(); ++i){
                            if (!isdigit((*chunk)[i])){
                                pass = false;
                                break;
                            }
                        }
                        if (pass){
                            lfound = true;
                        }
                        else{
                            return false;
                        }
                    }
                }
                else{
                    if (!sfound){
                        if ((*chunk)[0] == 'S'){
                            bool pass = true;
                            for (int i = 1; i < chunk->length(); ++i){
                                if (!isdigit((*chunk)[i])){
                                    pass = false;
                                    break;
                                }
                            }
                            if (pass){
                                sfound = true;
                            }
                            else{
                                return false;
                            }
                        } 
                    }
                    else{
                        // Made it.
                        if (fnbase != ""){
                            fnbase = *chunk + "_" + fnbase;
                        }
                        else{
                            fnbase = *chunk;
                        }
                    }
                }
            }
        } 
    }
    return true;
    */
}

/**
 * Write out all species-specific k-mer counts per barcode in table format.
 */
void print_bc_species_counts(
    robin_hood::unordered_map<unsigned long, map<short, int> >& bc_species_counts,
    map<short, string>& idx2species,
    FILE* countsfile){
    
    for (robin_hood::unordered_map<unsigned long, map<short, int> >::iterator bci = 
        bc_species_counts.begin(); bci != bc_species_counts.end(); ++bci){
        string bcstr = bc2str(bci->first);
        fprintf(countsfile, "%s", bcstr.c_str());
        for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
            int count = 0;
            if (bci->second.count(spec->first) > 0){
                count = bci->second[spec->first];
            }   
            fprintf(countsfile, "\t%d", (int)count);
        }
        fprintf(countsfile, "\n");
    }
}

/**
 * Instead of iterating through the BAM file again, can just load data from a previous run.
 * Requires counts file and species file both generated from previous runs.
 */
void load_from_files(string& countsfilename,
    string& speciesfilename,
    map<short, string>& idx2species,
    map<string, short>& species2idx,
    robin_hood::unordered_map<unsigned long, map<short, int> >& bc_species_counts){
    
    fprintf(stderr, "Loading barcode-species counts\n"); 
    // Load species names
    ifstream speciesfile(speciesfilename);
    string species_idx_str;
    string species_name;
    while (speciesfile >> species_idx_str >> species_name){
        short species_idx = atoi(species_idx_str.c_str());
        idx2species.insert(make_pair(species_idx, species_name));
        species2idx.insert(make_pair(species_name, species_idx));
    }
    
    // Load count data
    string line;
    ifstream countsfile(countsfilename);
    while (getline(countsfile, line)){
        istringstream splitter(line);
        string field;
        short field_idx = 0;
        
        bc bc_bin;
        while(getline(splitter, field, '\t')){
            if (field_idx == 0){
                // Barcode.
                if (field[field.length()-2] == '-' && field[field.length()-1] == '1'){
                    field = field.substr(0, field.length()-2);
                }
                str2bc(field.c_str(), bc_bin);   
                map<short, int> m;
                bc_species_counts.emplace(bc_bin.to_ulong(), m);
            }
            else if (field_idx <= species2idx.size()){
                // Field index == species ID + 1
                float count = atof(field.c_str());
                bc_species_counts[bc_bin.to_ulong()].insert(make_pair(field_idx-1, count));
            }
            else{
                // Optional last column = total # reads
                // Ignore for now
            }
            field_idx++;
        }
    }
    fprintf(stderr, "done\n");
}

/**
 * After assigning barcodes to species, print assignments to an output file that
 * can be reviewed.
 */
void print_assignments(FILE* outf,
    const string& libname, 
    bool cellranger,
    bool seurat,
    bool underscore,
    robin_hood::unordered_map<unsigned long, short>& bc2species,
    robin_hood::unordered_map<unsigned long, pair<unsigned int, unsigned int> >& bc2doublet,
    robin_hood::unordered_map<unsigned long, double>& bc2llr,
    map<short, string>& idx2species,
    bool use_filter,
    robin_hood::unordered_set<unsigned long>& filter){
    
    for (robin_hood::unordered_map<unsigned long, double>::iterator b2llr = bc2llr.begin();
        b2llr != bc2llr.end(); ++b2llr){
        
        if (use_filter && filter.find(b2llr->first) == filter.end()){
            continue;
        }

        char type = 'S';
        string name;
        if (bc2doublet.count(b2llr->first) > 0){
            type = 'D';
            int idx1 = bc2doublet[b2llr->first].first;
            int idx2 = bc2doublet[b2llr->first].second;
            if (idx2species[idx1] < idx2species[idx2]){
                name = idx2species[idx1] + "+" + idx2species[idx2];
            }
            else{
                name = idx2species[idx2] + "+" + idx2species[idx1];
            }
        }
        else if (bc2species.count(b2llr->first) > 0){
            type = 'S';
            name = idx2species[bc2species[b2llr->first]];
        }
        string bc_str = bc2str(b2llr->first);
        mod_bc_libname(bc_str, libname, cellranger, seurat, underscore); 
        fprintf(outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type,
            b2llr->second);
    }
}

/**
 * Create a 10x Genomics-style "Library file" to make it easier to run
 * cellranger-arc, if desired.
 */
void create_library_file(vector<string>& rna_r1files,
    vector<string>& atac_r1files,
    vector<string>& custom_r1files,
    vector<string>& custom_names,
    map<short, string>& idx2species,
    const string& outdir){

    char fullpath[200];
    //if (rna_r1files.size() > 0 && (atac_r1files.size() > 0 || custom_r1files.size() > 0)){
    if (true){  
        // Multiple data types. Create a "library file" for each species
        for (map<short, string>::iterator i2s = idx2species.begin(); i2s != idx2species.end(); ++i2s){
            string fnprefix = "";
            if (outdir == ""){
                //sprintf(&fullpath[0], "");
                fullpath[0] = '\0';
            }
            else{
                realpath(outdir.c_str(), &fullpath[0]);
                fnprefix = fullpath;
                if (fnprefix[fnprefix.length()-1] != '/'){
                    fnprefix += "/";
                }
            }
            string libfilename = fnprefix + i2s->second + ".library";
            if (file_exists(libfilename)){
                // Don't create one.
                // This allows us to, if running in batch mode, run once per file chunk to count,
                // then merge counts, then run once with --dump option to assign species & create
                // library file, then run once per pair of read files to demux each individually
                // (in parallel).
                continue;
            }

            fnprefix += i2s->second;
            FILE* libfile = fopen(libfilename.c_str(), "w");
            fprintf(libfile, "fastqs,sample,library_type\n");
            set<string> atac_base_pre;
            for (int i = 0; i < atac_r1files.size(); ++i){
                // Get base name
                string fn_trunc = filename_nopath(atac_r1files[i]); 
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                if (fn_base.length() < 5 || fn_base.substr(0, 5) != "ATAC_"){
                    fn_base = "ATAC_" + fn_base;
                }
                if (atac_base_pre.find(fn_base) == atac_base_pre.end()){
                    fprintf(libfile, "%s,%s,Chromatin Accessibility\n", fnprefix.c_str(),
                        fn_base.c_str());  
                }
                atac_base_pre.insert(fn_base);
            }
            set<string> rna_base_pre;
            for (int i = 0; i < rna_r1files.size(); ++i){
                string fn_trunc = filename_nopath(rna_r1files[i]);
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                if (fn_base.length() < 4 || fn_base.substr(0, 4) != "GEX_"){
                    fn_base = "GEX_" + fn_base;
                }
                // Ensure only unique entries are written
                if (rna_base_pre.find(fn_base) == rna_base_pre.end()){
                    fprintf(libfile, "%s,%s,Gene Expression\n", fnprefix.c_str(), 
                        fn_base.c_str());
                }
                rna_base_pre.insert(fn_base);
            }
            set<string> custom_base_pre;
            set<string> custom_type_pre;
            for (int i = 0; i < custom_r1files.size(); ++i){
                string fn_trunc = filename_nopath(custom_r1files[i]);
                string fn_base;
                if (!filename_base(fn_trunc, fn_base)){
                    fprintf(stderr, "ERROR: could not parse %s\n", fn_trunc.c_str());
                    exit(1);
                }
                // Determine what to append
                string prefix = custom_names[i] + "_";
                if (fn_base.length() < prefix.length() || 
                    fn_base.substr(0, prefix.length()) != prefix){
                    fn_base = prefix + fn_base;
                }
                string type_lib = custom_names[i];
                if (custom_names[i] == "CRISPR"){
                    type_lib = "CRISPR Guide Capture";
                }
                else if (custom_names[i] == "Ab"){
                    type_lib = "Antibody Capture";
                }
                if (custom_base_pre.find(fn_base) == custom_base_pre.end() ||
                    custom_type_pre.find(type_lib) == custom_type_pre.end()){
                    fprintf(libfile, "%s,%s,%s\n", fnprefix.c_str(), 
                        fn_base.c_str(), type_lib.c_str());
                }
                custom_base_pre.insert(fn_base);
                custom_type_pre.insert(type_lib);
            }
            fclose(libfile);
        }
    }
}
