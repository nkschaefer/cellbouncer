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
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/sam.h>
#include <zlib.h>
#include <htswrapper/bam.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bam_split_bcs [OPTIONS]\n");
    fprintf(stderr, "Given a BAM file and a file mapping cell barcodes to\n");
    fprintf(stderr, "identities (as output by demux_mt and demux_vcf), splits\n");
    fprintf(stderr, "the BAM file into individual BAM files per inferred\n");
    fprintf(stderr, "individual.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --assignments -a A file of the format output by\n");
    fprintf(stderr, "       demux_mt and demux_vcf, which should be tab\n");
    fprintf(stderr, "       separated, with the following columns:\n");
    fprintf(stderr, "       bc ID S/D LLR. ID is assignment, S/D indicates\n");
    fprintf(stderr, "       singlet or doublet, and LLR is the log likelihood\n");
    fprintf(stderr, "       ratio of best to second best assignment.\n");
    fprintf(stderr, "    --output_prefix -o The prefix for output bam files.\n");
    fprintf(stderr, "       Output files will be in the format \n");
    fprintf(stderr, "       <prefix>.<ID>.bam.\n");
    fprintf(stderr, "    --keep_doublets -d Create an output BAM file for\n");
    fprintf(stderr, "       each doublet combination in the assignment file?\n");
    fprintf(stderr, "       default = only include identified singlets.\n");
    fprintf(stderr, "    --llr -l Cutoff for log likelihood ratio of best to\n");
    fprintf(stderr, "       second best assignment. Cells with assignments\n");
    fprintf(stderr, "       below this threshold will be omitted from output\n");
    fprintf(stderr, "       files. Higher = fewer but more-confidently identified\n");
    fprintf(stderr, "       cells. default = 0 (use all assignments in file).\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"assignments", required_argument, 0, 'a'},
       {"keep_doublets", no_argument, 0, 'd'},
       {"llr", required_argument, 0, 'l'},
       {"output_prefix", required_argument, 0, 'o'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string output_prefix;
    string assignment_file;
    bool keep_doublets = false;
    double llr_thresh = 0;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:a:l:dh", long_options, &option_index )) != -1){
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
            case 'a':
                assignment_file = optarg;
                break;
            case 'd':
                keep_doublets = true;
                break;
            case 'l':
                llr_thresh = atof(optarg);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (bamfile.length() == 0){
        fprintf(stderr, "ERROR: bam file (--bam) required\n");
        exit(1);
    }
    if (assignment_file.length() == 0){
        fprintf(stderr, "ERROR: cell -> individual assignments (-a) required\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix (-o) required\n");
        exit(1);
    }
    
    map<string, string> barcode_map;
    set<string> barcode_groups;
    parse_barcode_map(assignment_file, barcode_map, barcode_groups, llr_thresh, keep_doublets);
    
    fprintf(stderr, "Writing %ld new BAM files\n", barcode_groups.size());
    
    map<string, BGZF* > outs;
    for (set<string>::iterator ind = barcode_groups.begin(); 
        ind != barcode_groups.end(); ++ind){
        // In case of a doublet combination, need to replace any + in
        // a name with an x to keep file names safe.
        string id = *ind;
        size_t ppos = id.find("+");
        if (ppos != string::npos){
            id = id.substr(0, ppos) + "x" + id.substr(ppos+1, id.length()-ppos-1);
        }
        char outbuf[100];
        sprintf(outbuf, "%s.%s.bam", output_prefix.c_str(), id.c_str());
        outs.insert(make_pair(*ind, bgzf_open(&outbuf[0], "w")));
    }   
    
    // Init BAM reader
    bam_reader reader(bamfile);
    
    // Write BAM headers
    for (map<string, BGZF*>::iterator outf = outs.begin(); outf != outs.end(); ++outf){
        // Clear existing read groups
        reader.clear_read_groups_hdr();
        // Add sample-specific read group
        reader.add_read_group_hdr(outf->first, outf->first, "NA", "NA", "Illumina");
        // Write header to file
        reader.write_header(outf->second);
    }
    reader.set_cb();
    
    while(reader.next()){
        if (reader.has_cb_z){
            if (barcode_map.count(reader.cb_z) > 0){
                string indv = barcode_map[reader.cb_z];
                // Add sample-specific read group
                reader.add_read_group_read(indv);
                // Write to specific output file
                reader.write_record(outs[indv]);
            }
        }
        
    }
    
    for (map<string, BGZF*>::iterator out = outs.begin(); out != outs.end(); ++out){
        bgzf_close(out->second);
    }

    return 0;
}
