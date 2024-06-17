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
    fprintf(stderr, "Given sc-ATAC-seq FASTQ files from 10X Genomics (or similar format) and an\n");
    fprintf(stderr, "allowed barcode list, finds valid cell barcodes in reads (within edit distance\n");
    fprintf(stderr, "1), trims barcodes from reads, and inserts barcodes into sequence comments.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs two FASTQ files: forward and reverse reads, with barcodes inserted into\n");
    fprintf(stderr, "both. Users can then map these to a reference genome using an aligner capable\n");
    fprintf(stderr, "of creating BAM tags from sequence comments. For example, minimap2 with the -y\n");
    fprintf(stderr, "flag will have this behavior, and the -a flag instructs it to output SAM.\n"); 
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --r1 -1 Read 1 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --r2 -2 Read 2 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --r3 -3 Read 3 FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --output_dir -o Directory in which to write output files. Files will have\n");
    fprintf(stderr, "      the same name as input files (R1 and R2), but in the output directory.\n");
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

    gzFile outs[2];
    outs[0] = gzopen(r1out.c_str(), "w");
    outs[1] = gzopen(r2out.c_str(), "w");

    bc_scanner scanner(r1fn, r2fn, r3fn);
    scanner.init_10x_ATAC(wlfn);
    scanner.trim_barcodes(true);

    char outbuf[1024];

    while(scanner.next()){
        
        // Print seq ID line
        string bc_str = bc2str(scanner.barcode);
        sprintf(&outbuf[0], "@%s CB:Z:%s\n", scanner.seq_id, bc_str.c_str());
        gzwrite(outs[0], &outbuf[0], 8+scanner.seq_id_len+BC_LENX2/2);
        gzwrite(outs[1], &outbuf[0], 8+scanner.seq_id_len+BC_LENX2/2);
        
        // Print sequences
        gzwrite(outs[0], scanner.read_f, scanner.read_f_len);
        gzwrite(outs[0], "\n+\n", 3);
        gzwrite(outs[1], scanner.read_r, scanner.read_r_len);
        gzwrite(outs[1], "\n+\n", 3);
        
        // Print quality
        gzwrite(outs[0], scanner.read_f_qual, scanner.read_f_len);
        gzwrite(outs[0], "\n", 1);
        gzwrite(outs[1], scanner.read_r_qual, scanner.read_r_len);
        gzwrite(outs[1], "\n", 1);
    
    }

    gzclose(outs[0]);
    gzclose(outs[1]);

    return 0;
}
