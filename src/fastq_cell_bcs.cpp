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
#include <unordered_set>
#include <cstdlib>
#include <utility>
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
    fprintf(stderr, "fastq_cell_bcs [OPTIONS]\n");
    fprintf(stderr, "Given a FASTQ file from 10X genomics and a barcode whitelist, finds\n");
    fprintf(stderr, "all cell barcodes matching reads in the FASTQ and prints to stdout.\n");
    fprintf(stderr, "Must provide R1; assumes 16 bp barcode at the beginning of the read.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --r1 -1 Forward read FASTQ file (can be gzipped)\n");
    fprintf(stderr, "    --whitelist -w Cellranger barcode whitelist\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"r1", required_argument, 0, 'r'},
       {"whitelist", required_argument, 0, 'w'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string r1fn = "";
    string wlfn = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "1:w:h", long_options, &option_index )) != -1){
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
    if (wlfn.length() == 0){
        fprintf(stderr, "ERROR: --whitelist / -w required\n");
        exit(1);
    }
    
    bc_scanner scanner(r1fn);
    scanner.init(wlfn, "", 0, false, false, false);
    
    unordered_set<unsigned long> bc_all;
    while(scanner.next()){
        bc_all.insert(scanner.barcode);
    }
    for (unordered_set<unsigned long>::iterator barcode = bc_all.begin(); barcode != bc_all.end();
        ++barcode){
        string bcstr = bc2str(*barcode);
        fprintf(stdout, "%s\n", bcstr.c_str());
    }
    return 0;
}
