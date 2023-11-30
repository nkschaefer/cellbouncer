#include <getopt.h>
#include <argp.h>
#include <string>
#include <algorithm>
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
#include <htslib/sam.h>
#include <htswrapper/bam.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bam_indiv_rg [OPTIONS]\n");
    fprintf(stderr, "Given a BAM file and a mapping of cell barcode to group ID,\n");
    fprintf(stderr, "creates read group tags in the BAM with the group ID \n");
    fprintf(stderr, "populating the sample field.\n");
    fprintf(stderr, "This allows for variant calling on this file, with inferred\n");
    fprintf(stderr, " individuals treated as separate samples.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --assignments -a A file (as output by demux_mt or \n");
    fprintf(stderr, "       demux_vcf) that maps cell barcode to individual of \n");
    fprintf(stderr, "       origin. Should be tab separated, with these columns:\n");
    fprintf(stderr, "       barcode ID S/D LLR, where ID is the assignment, S/D \n");
    fprintf(stderr, "       identifies singlets and doublets, and LLR is the log2\n");
    fprintf(stderr, "       likelihood ratio of best to second best assignment.\n");
    fprintf(stderr, "    --outfile -o The name of the output BAM file to create\n");
    fprintf(stderr, "        (default: stdout)\n");
    fprintf(stderr, "    --keep_doublets -d Add read group entries for doublet\n");
    fprintf(stderr, "       combinations (default: ignore doublets)\n");
    fprintf(stderr, "    --llr -l Cutoff for log likelihood ratio of best to next\n");
    fprintf(stderr, "       best assignment. Cells with assignments below this \n");
    fprintf(stderr, "       threshold will be omitted from all output BAM files.\n");
    fprintf(stderr, "       Higher LLR = fewer but higher confidence cells.\n");
    fprintf(stderr, "       default = 0 (use all assignments in file)\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"assignments", required_argument, 0, 'a'},
       {"keep_doublets", no_argument, 0, 'd'},
       {"llr", required_argument, 0, 'l'},
       {"outfile", required_argument, 0, 'o'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string outfile;
    string assignment_file;
    bool keep_doublets = false;
    double llr_thresh = 0.0;;

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
            case 'o':
                outfile = optarg;
                break;
            case 'd':
                keep_doublets = true;
                break;
            case 'l':
                llr_thresh = atof(optarg);
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
        fprintf(stderr, "ERROR: barcode->individual identity assignment \
file (-a) required\n");
        exit(1);
    }
    BGZF* outf;
    if (outfile.length() == 0){
        // Print to stdout
        outf = bgzf_open("-", "w");
    }
    else{
        outf = bgzf_open(outfile.c_str(), "w");
        if (!outf){
            fprintf(stderr, "ERROR opening file %s for writing.\n", outfile.c_str());
            exit(1);
        }
    }
    
    map<unsigned long, string> barcode_map;
    set<string> barcode_groups;
    parse_barcode_map(assignment_file, barcode_map, barcode_groups, 
        llr_thresh, keep_doublets);
    
    fprintf(stderr, "Read %ld identities from file %s\n", barcode_groups.size(), 
        assignment_file.c_str());
    fprintf(stderr, "%ld cell barcodes included\n", barcode_map.size());

    // Init BAM reader
    bam_reader reader(bamfile);
    
    // Remove existing read groups from header
    reader.clear_read_groups_hdr();

    // Add read groups to header
    for (set<string>::iterator sm = barcode_groups.begin(); 
        sm != barcode_groups.end(); ++sm){
        reader.add_read_group_hdr(*sm, *sm, "NA", "NA", "Illumina"); 
    }
    
    // Write BAM header
    reader.write_header(outf);
    
    reader.set_cb();
    
    while(reader.next()){
        // Skip records for which there is no barcode -> individual assignment
        if (reader.has_cb_z){
            bc as_bitset;
            str2bc(reader.cb_z, as_bitset, 16);
            unsigned long as_ulong = as_bitset.to_ulong();
            if (barcode_map.count(as_ulong) > 0){
                reader.add_read_group_read(barcode_map[as_ulong]);
                reader.write_record(outf);
            }
        }
        
    }
    
    bgzf_close(outf);

    return 0;
}
