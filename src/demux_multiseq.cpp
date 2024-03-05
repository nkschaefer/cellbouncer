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
#include <utility>
#include <math.h>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/bc_scanner.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_multiseq [OPTIONS]\n");
    fprintf(stderr, "Given reads containing MULTIseq data, counts the occurrences of each\n");
    fprintf(stderr, "  MULTIseq barcode per cell barcode and assigns each cell barcode\n");
    fprintf(stderr, "  an identity based on MULTIseq barcode meaning.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n   ===== REQUIRED =====\n");
    fprintf(stderr, "   --output_prefix -o The prefix for output file names. Will create\n");
    fprintf(stderr, "       a file ending in .counts and another ending in .assignments.\n");
    fprintf(stderr, "   --read1 -1 Forward read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --read2 -2 Reverse read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --whitelist -w Cell barcode whitelist file (i.e. included with 10X Genomics\n");
    fprintf(stderr, "       Cellranger software in cellranger-x.y.z/lib/python/cellranger/barcodes/\n");
    fprintf(stderr, "\n   ===== STRONGLY RECOMMENDED =====\n");
    fprintf(stderr, "   --mapping -m 2-column tab separated file mapping string interpretation of\n");
    fprintf(stderr, "       MULTIseq label (i.e. unique identifier, or treatment name) to MULTIseq\n");
    fprintf(stderr, "       well identifier.\n");
    fprintf(stderr, "\n   ===== OPTIONAL =====\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    fprintf(stderr, "   --barcodes -b A path to a file listing MULTIseq well IDs and corresponding\n");
    fprintf(stderr, "       barcodes, tab separated. If you do not provide one, the default file (in\n");
    fprintf(stderr, "       the data directory) will be used.\n");
    fprintf(stderr, "   --doublet_rate -D What is the prior expected doublet rate?\n");
    fprintf(stderr, "       (OPTIONAL; default = 0.1). Must be a decimal between 0 and 1,\n");
    fprintf(stderr, "       exclusive.\n");
    exit(code);
}

/**
 * Returns barcode length.
 */
int load_well2bc(string& filename, 
    map<unsigned long, string>& bc2well, 
    vector<string>& bclist){

    ifstream inf(filename);
    string bc_str;
    string well;
    bc bc_bitset;
    int bc_len = -1;
    while (inf >> bc_str >> well){
        if (bc_len == -1){
            bc_len = bc_str.length();
        }
        else if (bc_len != bc_str.length()){
            fprintf(stderr, "ERROR: mismatching MULTIseq barcode lengths in file %s:\n", filename.c_str());
            fprintf(stderr, "%d vs %ld\n", bc_len, bc_str.length());
            exit(1);
        }
        str2bc(bc_str.c_str(), bc_bitset, bc_len);
        bc2well.insert(make_pair(bc_bitset.to_ulong(), well));
        bclist.push_back(bc_str);
    }
    return bc_len;
}

void load_well_mapping(string& filename, map<string, string>& well2id){
    ifstream inf(filename);
    string well;
    string id;
    while (inf >> well >> id){
        well2id.insert(make_pair(well, id));
    }
}

void dump_counts(string& filename, 
    robin_hood::unordered_map<unsigned long, map<unsigned long, umi_set > >& bc_ms_umis,
    map<unsigned long, string>& bc2well,
    map<string, string>& well2id){
    
    FILE* f = fopen(filename.c_str(), "w");
    
    // Print header line
    fprintf(f, "bc");
    for (map<unsigned long, string>::iterator bw = bc2well.begin(); bw != bc2well.end(); ++bw){
        fprintf(f, "\t%s", well2id[bw->second].c_str());
    }
    fprintf(f, "\n");

    for (robin_hood::unordered_map<unsigned long, map<unsigned long, umi_set> >::iterator x =
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){
        string bc_str = bc2str(x->first);
        fprintf(f, "%s", bc_str.c_str());
        for (map<unsigned long, string>::iterator bw = bc2well.begin(); bw != bc2well.end(); ++bw){
            if (x->second.count(bw->first) > 0){
                fprintf(f, "\t%d", x->second[bw->first].count());
            }
            else{
                fprintf(f, "\t0");
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
     
int main(int argc, char *argv[]) {    
    
    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"mapping", required_argument, 0, 'm'},
       {"whitelist", required_argument, 0, 'w'},
       {"barcodes", required_argument, 0, 'b'},
       {"doublet_rate", required_argument, 0, 'D'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string output_prefix = "";
    vector<string> read1fn;
    vector<string> read2fn;
    string mapfn = "";
    string wlfn = "";
    string barcodesfn = PROJ_ROOT + "/data/multiseq_indices.txt";
    double doublet_rate = 0.1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:1:2:m:w:b:D:h", 
        long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case '1':
                read1fn.push_back(optarg);
               break;
            case '2':
                read2fn.push_back(optarg);
                break;
            case 'm':
                mapfn = optarg;
                break;
            case 'w':
                wlfn = optarg;
                break;
            case 'b':
                barcodesfn = optarg;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: --output_prefix/-o is required\n");
        exit(1);
    }
    else if (read1fn.size() == 0 || read2fn.size() == 0){    
        fprintf(stderr, "ERROR: input FASTQs for reads (-1 and -2) are required.\n");
        exit(1);
    }
    else if (read1fn.size() != read2fn.size()){
        fprintf(stderr, "ERROR: unequal number of forward (-1) and reverse (-2) reads provided\n");
        exit(1);
    }
    else if (wlfn == ""){
        fprintf(stderr, "ERROR: cell barcode whitelist / -w is required\n");
        exit(1);
    }
    
    // Load MULTIseq barcode -> well ID mapping
    map<unsigned long, string> bc2well;
    vector<string> bclist;
    int ms_bc_len = load_well2bc(barcodesfn, bc2well, bclist);
    
    // Determine appropriate k-mer length for fuzzy matching
    // k must be < (bc_length + 1)/2
    int ms_k;
    if (ms_bc_len/2 % 2 == 0){
        // Even barcode length.
        ms_k = ms_bc_len/2;
    }
    else{
        ms_k = (ms_bc_len-1)/2;
    }

    // Load MULTIseq well -> unique identifier mapping
    map<string, string> well2id;
    if (mapfn == ""){
        fprintf(stderr, "WARNING: no MULTIseq well to identifier mapping provided.\n");
        fprintf(stderr, "Using well names as identifiers.\n");    
        for (map<unsigned long, string>::iterator bw = bc2well.begin(); 
            bw != bc2well.end(); ++bw){
            if (well2id.count(bw->second) == 0){
                well2id.insert(make_pair(bw->second, bw->second));
            }
        }
    }
    else{
        load_well_mapping(mapfn, well2id);
        // Make sure wells specified here have associated valid barcodes.
        set<string> wellnames;
        for (map<unsigned long, string>::iterator bw = bc2well.begin(); bw != 
            bc2well.end(); ++bw){
            wellnames.insert(bw->second);
        }
        for (map<string, string>::iterator wi = well2id.begin(); wi != well2id.end();
            ++wi){
            if (wellnames.find(wi->first) == wellnames.end()){
                fprintf(stderr, "ERROR: well -> ID mapping file contains a well named %s\n", wi->first.c_str());
                fprintf(stderr, "This well ID is not present in the MULTIseq barcode -> well mappings\n");
                exit(1);
            }
        }
        // If the user has omitted one or more wells that turn up often in the data, we will
        // alert the user about this later.
    }
    
    // Initiate MULTIseq barcode whitelist
    bc_whitelist wl_ms(bclist, ms_bc_len, ms_k);
    
    // Data structure to store MULTIseq barcode counts per cell barcode
    // counts come from UMIs
    robin_hood::unordered_map<unsigned long, map<unsigned long, umi_set > > bc_ms_umis;
    
    int umi_len = 10;
    
    // Set up object that will scan each pair of read files
    bc_scanner scanner;
    scanner.init_multiseq_v3(wlfn);

    for (int i = 0; i < read1fn.size(); ++i){
        
        // Initiate object to read through FASTQs and find cell barcodes.
        scanner.add_reads(read1fn[i], read2fn[i]);
        
        while (scanner.next()){
            // At this point, there's a valid cell barcode.
            
            // scanner.barcode_read holds R1
            // scanner.read_f holds R2
            // scanner.umi holds UMI 
            if (scanner.has_umi){
                
                // MULTIseq barcodes are the first 8 bp (or however long supplied barcodes are)
                // of R2, forward orientation
                unsigned long ms_ul;
                if (wl_ms.lookup(scanner.read_f, false, ms_ul)){
                    // Store UMI
                    umi this_umi(scanner.umi, scanner.umi_len);    
                    if (bc_ms_umis.count(scanner.barcode) == 0){
                        map<unsigned long, umi_set> m;
                        bc_ms_umis.emplace(scanner.barcode, m);
                    }
                    if (bc_ms_umis[scanner.barcode].count(ms_ul) == 0){
                        umi_set s(scanner.umi_len);
                        bc_ms_umis[scanner.barcode].insert(make_pair(ms_ul, s));
                    }
                    bc_ms_umis[scanner.barcode].at(ms_ul).add(this_umi);
                }
            }
        } 
    }
    
    // Write counts to disk
    string countsfn = output_prefix + ".counts";
    dump_counts(countsfn, bc_ms_umis, bc2well, well2id);

    return 0;
}
