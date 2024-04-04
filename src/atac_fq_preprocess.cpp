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
#include <cstdlib>
#include <utility>
#include <sys/stat.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bc_scanner.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "atac_fq_preprocess [OPTIONS]\n");
    fprintf(stderr, "Given sc-ATAC-seq FASTQ files from 10X Genomics and a barcode whitelist,\n");
    fprintf(stderr, "checks for whitelisted barcodes, trims them, and inserts them into sequence\n");
    fprintf(stderr, "comments so they can be aligned with an aligner that is capable of inserting\n");
    fprintf(stderr, "sequence comments as BAM tags.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --r1 -1 Read 1 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --r2 -2 Read 2 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --r3 -3 Read 3 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --output_dir -o Directory in which to write output files. Files will have\n");
    fprintf(stderr, "      the same name as input files, but in the output directory.\n");
    fprintf(stderr, "    --whitelist -w Cellranger barcode whitelist\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"r1", required_argument, 0, '1'},
       {"r2", required_argument, 0, '2'},
       {"r3", required_argument, 0, '3'},
       {"output_dir", required_argument, 0, 'o'},
       {"whitelist", required_argument, 0, 'w'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string r1fn = "";
    string r2fn = "";
    string r3fn = "";
    string output_dir = "";
    string wlfn = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "1:2:3:o:w:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case '1':
                r1fn = optarg;
                break;
            case '2':
                r2fn = optarg;
                break;
            case '3':
                r3fn = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'w':
                wlfn = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (r1fn.length() == 0){
        fprintf(stderr, "ERROR: --r1 / -1 required\n");
        exit(1);
    }
    if (r2fn.length() == 0){
        fprintf(stderr, "ERROR: --r2 / -2 required\n");
        exit(1);
    }
    if (r3fn.length() == 0){
        fprintf(stderr, "ERROR: --r3 / -3 required\n");
        exit(1);
    }
    if (output_dir.length() == 0){
        fprintf(stderr, "ERROR: --output_dir / -o required\n");
        exit(1);
    }
    if (wlfn.length() == 0){
        fprintf(stderr, "ERROR: --whitelist / -w required\n");
        exit(1);
    }
    
    if (output_dir[output_dir.length()-1] == '/'){
        output_dir = output_dir.substr(0, output_dir.length()-1);
    }

    // Create output directory if it doesn't exist
    if (!mkdir(output_dir.c_str(), 0775)){
        // Assume directory already exists
    }
    output_dir += "/";

    string r1out = output_dir + filename_nopath(r1fn);
    string r2out = output_dir + filename_nopath(r2fn);
    string r3out = output_dir + filename_nopath(r3fn);

    gzFile outs[3];
    outs[0] = gzopen(r1out.c_str(), "w");
    outs[1] = gzopen(r2out.c_str(), "w");
    outs[2] = gzopen(r3out.c_str(), "w");

    bc_scanner scanner(r1fn, r2fn, r3fn);
    scanner.init_10x_ATAC(wlfn);
    scanner.trim_barcodes(true);

    while(scanner.next()){
        fprintf(stderr, "%s\n", scanner.seq_id);
        fprintf(stderr, "%s\n", scanner.read_f);
        fprintf(stderr, "%s\n", scanner.read_r);
        string bc_str = bc2str(scanner.barcode);
        fprintf(stderr, "%s\n", bc_str.c_str());
        fprintf(stderr, "%s\n", scanner.barcode_read);
        exit(0); 
    }

    gzclose(outs[0]);
    gzclose(outs[1]);
    gzclose(outs[2]);

    return 0;
}
