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
#include <sys/stat.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <random>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bam.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * This file extracts counts from a composite reference genome and converts them
 * into a table that can be interpreted by demux_species.
 */

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "composite_bam2counts [OPTIONS]\n");
    fprintf(stderr, "Given reads mapped to a composite reference genome, extracts the number of reads\n");
    fprintf(stderr, "  per cell uniquely mapped to each species' reference and outputs a table of counts\n");
    fprintf(stderr, "  and species names, which can then be read and used by demux_species.\n");
    fprintf(stderr, "Assumes that you have prepended or appended a unique identifier to each chromosome/scaffold\n");
    fprintf(stderr, "  name indicating which species it belongs to: for instance, if you have concatenated\n");
    fprintf(stderr, "  the mouse genome mm10 and human genome hg38, you could prepend \"mm10_\" to all mouse\n");
    fprintf(stderr, "  chromosome names and \"hg38_\" to all human chromosome names.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file containing reads mapped to a composite reference.\n");
    fprintf(stderr, "       NOTE: this must be for a single library only: i.e. do not merge BAMs from multiple\n");
    fprintf(stderr, "       10X lanes into a single BAM. Run each separately.\n");
    fprintf(stderr, "    --output_directory Will be created if it does not exist. This program will write files\n");
    fprintf(stderr, "       species_counts.txt and species_names.txt to this directory, which can then be read by\n");
    fprintf(stderr, "       demux_species, passing this argument as the -o option.\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --separator -s What separator was used when prepending/appending unique species IDs\n");
    fprintf(stderr, "       to chromosome/scaffold names? Default = _\n");
    fprintf(stderr, "    --end -e By default, this program assumes unique species IDs are prepended onto\n");
    fprintf(stderr, "       beginnings of chromosome/scaffold names (species ID + separator + sequence name).\n");
    fprintf(stderr, "       With this option set, species IDs are assumed to be appended onto the ends of\n");
    fprintf(stderr, "       chromosome/scaffold names (sequence name + separator + species ID).\n");
    //print_libname_help();
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

bool parse_tid2species(map<int, int>& tid2species,
    vector<string>& species,
    map<string, int>& chrom2tid,
    const string& separator,
    bool suffix){
    
    tid2species.clear();
    species.clear();

    int seplen = separator.length();

    set<string> species_all;

    map<int, string> tid2species_tmp;

    for (map<string, int>::iterator c2t = chrom2tid.begin(); c2t != chrom2tid.end(); ++c2t){
        size_t pos;
        if (suffix){
            pos = c2t->first.rfind(separator);
            if (pos != string::npos){
                string suf = c2t->first.substr(pos+seplen, c2t->first.length()-pos-seplen);
                tid2species_tmp.insert(make_pair(c2t->second, suf));
                species_all.insert(suf);
            }
        }
        else{
            pos = c2t->first.find(separator);
            if (pos != string::npos){
                string pref = c2t->first.substr(0, pos);
                tid2species_tmp.insert(make_pair(c2t->second, pref));
                species_all.insert(pref);
            }
        }    
    }
    
    if (species_all.size() > 0){
        fprintf(stderr, "Found species:\n");
        for (set<string>::iterator s = species_all.begin(); s != species_all.end(); ++s){
            fprintf(stderr, "%s\n", s->c_str());
        }

        map<string, int> species2idx;
        int idx = 0;
        for (set<string>::iterator s = species_all.begin(); s != species_all.end(); ++s){
            species2idx.insert(make_pair(*s, idx));
            species.push_back(*s);
            ++idx;
        }
        for (map<int, string>::iterator ts = tid2species_tmp.begin(); ts != tid2species_tmp.end(); ++ts){
            tid2species.insert(make_pair(ts->first, species2idx[ts->second]));
        }
        return true;
    }
    return false;
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"output_directory", required_argument, 0, 'o'},
       {"separator", required_argument, 0, 's'},
       {"end", no_argument, 0, 'e'},
       //{"libname", required_argument, 0, 'n'},
       //{"cellranger", no_argument, 0, 'C'},
       //{"seurat", no_argument, 0, 'S'},
       //{"underscore", no_argument, 0, 'U'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string output_prefix = "";
    string separator = "_";
    bool suffix = false;

    string libname = "";
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:s:n:eCSUh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'b':
                bamfile = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 's':
                separator = optarg;
                break;
            case 'e':
                suffix = true;
                break;
            case 'n':
                libname = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'S':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: --output_directory/-o required\n");
        exit(1);
    }

    // Init BAM reader
    bam_reader reader = bam_reader();

    if (bamfile.length() == 0){
        fprintf(stderr, "ERROR: bam file (--bam) required\n");
        exit(1);
    }    
    reader.set_file(bamfile);

    // Print progress message every n sites
    long int progress = 50000;

    // retrieve cell barcodes
    reader.set_cb();
    
    // Get a mapping of chromosome names to internal numeric IDs in the 
    // BAM file.
    map<string, int> chrom2tid = reader.get_seq2tid();
    
    // Now map TIDs to species
    map<int, int> tid2species;
    vector<string> species;

    if (!parse_tid2species(tid2species, species, chrom2tid, separator, suffix)){
        fprintf(stderr, "ERROR: no species names found. Did you specify the wrong --separator/-s \
and/or --suffix/-S?\n");
        exit(1);
    }

    long int reads_processed = 0;

    robin_hood::unordered_map<unsigned long, map<int, int> > counts;

    // Read through entire BAM file and look for informative SNPs along the way
    // (default setting, appropriate for large numbers of SNPs).
    
    while (reader.next()){
        
        if (reader.has_cb_z && !reader.unmapped() && !reader.secondary() && !reader.dup()){
            int tid = reader.tid();
            if (tid2species.count(tid) > 0){
                int species_idx = tid2species[tid];
                unsigned long ul = bc_ul(reader.cb_z);
                if (counts.count(ul) == 0){
                    map<int, int> m;
                    counts.emplace(ul, m);
                }
                if (counts[ul].count(species_idx) == 0){
                    counts[ul].insert(make_pair(species_idx, 0));
                }
                counts[ul][species_idx]++;
            }
        }
        
        ++reads_processed;
        if (reads_processed % progress == 0){
            fprintf(stderr, "Processed %ld reads\r", reads_processed);
        }
    }
    fprintf(stderr, "Processed %ld reads\n", reads_processed);
    
    // Write output data
    if (output_prefix[output_prefix.length()-1] == '/'){
        output_prefix = output_prefix.substr(0, output_prefix.length()-1);
    }

    if (!mkdir(output_prefix.c_str(), 0775)){
        // Assume directory already exists
    }
    
    string out_counts_name = output_prefix + "/species_counts.txt";
    string out_names_name = output_prefix + "/species_names.txt";

    FILE* out_counts_f = fopen(out_counts_name.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, map<int, int> >::iterator c = counts.begin();
        c != counts.end(); ++c){
        string bc_str = bc2str(c->first);
        fprintf(out_counts_f, "%s", bc_str.c_str());
        for (int i = 0; i < species.size(); ++i){
            int s_count = 0;
            if (c->second.count(i) > 0){
                s_count = c->second[i];
            }
            fprintf(out_counts_f, "\t%d", s_count);
        }
        fprintf(out_counts_f, "\n");
    }
    fclose(out_counts_f);

    FILE* out_names_f = fopen(out_names_name.c_str(), "w");
    for (int i = 0; i < species.size(); ++i){
        fprintf(out_names_f, "%d\t%s\n", i, species[i].c_str());
    }
    fclose(out_names_f);


    return 0;
}
