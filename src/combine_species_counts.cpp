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
#include <map>
#include <set>
#include <cstdlib>
#include <cstdio>
#include <utility>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "combine_species_counts [OPTIONS]\n");
    fprintf(stderr, "After you have run demux_species multiple times with different batch\n");
    fprintf(stderr, "   numbers, run this program to combine all of the batch-specific output\n");
    fprintf(stderr, "   files into merged files. You can then re-run demux_species, which will\n");
    fprintf(stderr, "   load the combined counts.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "   --output_dir -o The output directory used for all batch-specific runs.\n");
    fprintf(stderr, "   --num_chunks -n The number of batches (there should be this many output\n");
    fprintf(stderr, "       files, unless one or more failed.\n");
    exit(code);
}

bool load_species(string& filename, map<short, string>& species_prev){
    map<short, string> species_this;
    ifstream inf(filename.c_str());
    short idx;
    string name;
    while (inf >> idx >> name){
        species_this.insert(make_pair(idx, name));
    }
    if (species_prev.size() > 0){
        for (map<short, string>::iterator sp = species_prev.begin(); sp != species_prev.end(); ++sp){
            if (species_this[sp->first] != sp->second){
                return false;
            }
        }
    }
    else{
        species_prev = species_this;
    }
    return true;
}

void load_counts(string& filename, 
    robin_hood::unordered_map<unsigned long, map<short, int> >& bc_species_counts){
    
    ifstream inf(filename.c_str());
    string line;
    while (getline(inf, line)){
        istringstream splitter(line);
        string elt;
        int idx = 0;
        unsigned long ul;
        while (getline(splitter, elt, '\t')){
            if (idx == 0){
                ul = bc_ul(elt);
                if (bc_species_counts.count(ul) == 0){
                    map<short, int> m;
                    bc_species_counts.emplace(ul, m);
                }
            }
            else{
                int count = atoi(elt.c_str());
                int si = idx-1;
                if (bc_species_counts[ul].count(si) == 0){
                    bc_species_counts[ul].insert(make_pair(si, 0));
                }
                bc_species_counts[ul][si] += count;
            }
            ++idx;
        }
    }
}

bool load_conversion(string& filename, 
    robin_hood::unordered_map<unsigned long, unsigned long>& conv){
    
    ifstream inf(filename.c_str());
    unsigned long c1;
    unsigned long c2;
    while (inf >> c1 >> c2){
        if (conv.count(c1) > 0 && conv[c1] != c2){
            return false;
        }
        else if (conv.count(c1) == 0){
            conv.emplace(c1, c2);
        }
    } 
    return true;
}

int main(int argc, char *argv[]) {    
    
    // Define long-form program options 
    static struct option long_options[] = {
       {"output_directory", required_argument, 0, 'o'},
       {"num_chunks", required_argument, 0, 'n'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string outdir = "";
    int nchunks = -1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:n:h", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                outdir = optarg;
                break;
            case 'n':
                nchunks = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (outdir == ""){
        fprintf(stderr, "ERROR: output_directory / -o is required\n");
        exit(1);
    }
    else if (outdir == "."){
        outdir = "";
    } 
    else if (outdir != "" && outdir[outdir.length() - 1] != '/'){
        outdir += "/";
    }
    if (nchunks <= 0){
        fprintf(stderr, "ERROR: --num_chunks/-n must be positive.\n");
        exit(1);
    }

    robin_hood::unordered_map<unsigned long, map<short, int> > bc_species_counts;
    
    // Multiome data uses different ATAC and RNA-seq barcodes
    // This maps an RNA-seq barcode (which goes into the BAM) to an ATAC-seq barcode
    robin_hood::unordered_map<unsigned long, unsigned long> bc_conversion;
    
    bool has_conversion = false;

    map<short, string> species_names;

    char buf[50];
    
    string counts_out = outdir + "species_counts.txt";
    string species_out = outdir + "species_names.txt";
    string conv_out = outdir + "bcmap.txt";
    
    vector<string> rm;
    
    for (int idx = 1; idx <= nchunks; ++idx){
        sprintf(&buf[0], "%d", idx);
        string bufstr = buf;
        string countsfilename = outdir + "species_counts." + bufstr + ".txt";
        string speciesfilename = outdir + "species_names." + bufstr + ".txt";
        string convfilename = outdir + "bcmap." + bufstr + ".txt";
        if (file_exists(countsfilename)){
            load_counts(countsfilename, bc_species_counts);
            rm.push_back(countsfilename);
            if (file_exists(speciesfilename)){
                bool ok = load_species(speciesfilename, species_names);
                if (!ok){
                    fprintf(stderr, "ERROR: species names from batch %d do not match previous \
species names.\n", idx);
                    exit(1);
                }
                rm.push_back(speciesfilename);
            }
            else{
                fprintf(stderr, "ERROR: missing expected file %s\n", speciesfilename.c_str());
                exit(1);
            }
            if (file_exists(convfilename)){
                bool ok = load_conversion(convfilename, bc_conversion);                
                if (!ok){
                    fprintf(stderr, "ERROR: barcode conversion file for batch %d does not match \
conversions from previous batches.\n", idx);
                    exit(1);
                }
                has_conversion = true;
                rm.push_back(convfilename);
            }
            else if (has_conversion){
                fprintf(stderr, "ERROR: barcode conversion file not found for batch %d\n", idx);
                exit(1);
            }
        }
        else{
            fprintf(stderr, "ERROR: missing expected file %s\n", countsfilename.c_str());
            fprintf(stderr, "One or more jobs likely failed to complete. Please re-run demux_species\n");
            fprintf(stderr, "on all missing batches.\n");
            exit(1);
        }
        ++idx;
    } 
   
    // Now write combined data to disk
    
    FILE* outf_s = fopen(species_out.c_str(), "w");
    for (map<short, string>::iterator s = species_names.begin(); s != 
        species_names.end(); ++s){
        fprintf(outf_s, "%d\t%s\n", s->first, s->second.c_str());
    }
    fclose(outf_s);
    
    FILE* outf_c = fopen(counts_out.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, map<short, int> >::iterator bsc = 
        bc_species_counts.begin(); bsc != bc_species_counts.end(); ++bsc){
        string bc_str = bc2str(bsc->first);
        fprintf(outf_c, "%s", bc_str.c_str());
        for (map<short, int>::iterator m = bsc->second.begin(); m != bsc->second.end(); ++m){
            fprintf(outf_c, "\t%d", m->second);
        }
        fprintf(outf_c, "\n");
    }
    fclose(outf_c);

    if (has_conversion){
        FILE* outf = fopen(conv_out.c_str(), "w");
        for (robin_hood::unordered_map<unsigned long, unsigned long>::iterator c = bc_conversion.begin();
            c != bc_conversion.end(); ++c){
            fprintf(outf, "%ld\t%ld\n", c->first, c->second);
        }
        fclose(outf);    
    }

    // Clean up
    for (int i = 0; i < rm.size(); ++i){
        remove(rm[i].c_str());
    }

    return 0;
}
