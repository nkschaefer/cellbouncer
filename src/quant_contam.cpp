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
#include <zlib.h>
#include "robin_hood.h"
#include "common.h"
#include "ambient_rna.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

// ===== Program to profile ambient RNA contamination in cells, =====
//       given output of a demux_vcf run.

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "quant_contam [OPTIONS]\n");
    fprintf(stderr, "Given output of a demux_vcf run, uses the (pre-computed) allele counts\n");
    fprintf(stderr, "to model ambient RNA contamination. Outputs the estimated fraction of\n");
    fprintf(stderr, "each cell's RNA composed of ambient RNA and attempts to find the likeliest\n");
    fprintf(stderr, "mixture of individuals from the VCF that compose the pool of ambient RNA.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --output_prefix -o The output prefix used in a prior run of demux_vcf\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --llr -l Log likelihood ratio cutoff to filter assignments from demux_vcf.\n");
    fprintf(stderr, "       This is the fourth column in the .assignments file. Default = 0 (no filter)\n");
    fprintf(stderr, "    --disable_profile -p In addition to quantifying the ambient RNA contamination\n");
    fprintf(stderr, "        per cell, default behavior is to model the ambient RNA contamination\n");
    fprintf(stderr, "        as a mixture of individuals. This involves learning the rate at which\n");
    fprintf(stderr, "        each individual matches other individuals' homozygous ref, heterozygous,\n");
    fprintf(stderr, "        and homozygous alt alleles, which is slow. This flag disables this step.\n");
    fprintf(stderr, "    --max_cells -c If you do not disable the step (using -p above), speed up the\n");
    fprintf(stderr, "        slow step of modeling by limiting the number of cells assigned to each\n");
    fprintf(stderr, "        individual that will be used in learning alt allele matching fractions.\n");
    fprintf(stderr, "        Specify an integer here (default 50), or set <= 0 to use the whole data set.\n");
    fprintf(stderr, "    --other_species -s If profiling the ambient RNA is enabled (no -p option),\n");
    fprintf(stderr, "        and your data came from a pool of multiple species, demultiplexed and each\n");
    fprintf(stderr, "        mapped to its species-specific reference genome and then demultiplexed by\n");
    fprintf(stderr, "        individual using within-species SNPs, this option models ambient RNA as a\n");
    fprintf(stderr, "        mixture of all individuals in the VCF, plus RNA from other species.\n");
    fprintf(stderr, "    --error_ref -e The underlying, true rate of misreading reference as\n");
    fprintf(stderr, "        alt alleles (should only reflect sequencing error if variant calls\n");
    fprintf(stderr, "        are reliable; default 0.001)\n");
    fprintf(stderr, "    --error_alt -E The underlying, true rate of misreading alt as reference\n");
    fprintf(stderr, "        alleles (should only reflect sequencing error if variant calls are\n");
    fprintf(stderr, "        reliable; default 0.001)\n"); 
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    

    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"disable_profile", no_argument, 0, 'p'},
       {"other_species", required_argument, 0, 's'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"max_cells", required_argument, 0, 'c'},
       {"llr", required_argument, 0, 'l'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string output_prefix = "";
    bool disable_profile = false;
    bool inter_species = false;
    double error_ref = 0.001;
    double error_alt = 0.001;
    int max_cells = 50;
    double llr = 0.0;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:s:e:E:c:l:ph", long_options, &option_index )) != -1){
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
            case 's':
                inter_species = true;
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 'c':
                max_cells = atoi(optarg);
                break;
            case 'p':
                disable_profile = true;
                break;
            case 'l':
                llr = atof(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
     
    // Error check arguments.
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: output_prefix required\n");
        exit(1);
    }
    if (error_ref <= 0 || error_ref >= 1.0 || error_alt <= 0 || error_alt >= 1.0){
        fprintf(stderr, "ERROR: error rates must be between 0 and 1, exclusive.\n");
        exit(1);
    }
    
    // Load files
    string sample_name = output_prefix + ".samples";
    vector<string> samples;
    if (file_exists(sample_name)){
        load_samples(sample_name, samples);
    }
    else{
        fprintf(stderr, "ERROR: no samples file found for %s. Please run demux_vcf with\n", 
            output_prefix.c_str());
        fprintf(stderr, "same output prefix.\n");
        exit(1);
    }

    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    string barcode_group;

    string assn_name = output_prefix + ".assignments";
    if (file_exists(assn_name)){
        fprintf(stderr, "Loading assignments...\n");
        load_assignments_from_file(assn_name, assn, assn_llr, samples, barcode_group);
        if (llr > 0.0){
            // Filter assignments.
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end();){
                if (assn_llr[a->first] <= llr){
                    assn_llr.erase(a->first);
                    // Can't erase while iterating the normal way
                    // https://github.com/martinus/robin-hood-hashing/issues/18
                    a = assn.erase(a);
                }
                else{
                    ++a;
                }
            }
            if (assn.size() == 0){
                fprintf(stderr, "ERROR: LLR filter too high; no assignments left to use.\n");
                exit(1);
            }
        }
    }
    else{
        fprintf(stderr, "ERROR: no assignments found for %s. Please run demux_vcf with same\n", 
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1); 
    }
    
    
    // Store allele counts
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    string counts_name = output_prefix + ".counts";
    if (file_exists(counts_name)){
        load_counts_from_file(indv_allelecounts, samples, counts_name); 
    }
    else{
        fprintf(stderr, "ERROR: no counts found for %s. Please run demux_vcf with same\n",
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1);
    }

    contamFinder cf(indv_allelecounts, assn, assn_llr, samples.size());
    
    cf.set_error_rates(error_ref, error_alt);
    if (inter_species){
        cf.model_other_species();
    } 
    if (disable_profile){
        cf.skip_model_mixture();
    }
    else{
        cf.max_cells_for_expfracs(max_cells);
    }
    cf.fit(); 
        
    // Write contamination profile to disk
    if (!disable_profile){
        string fname = output_prefix + ".contam_prof";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing contamination profile to disk...\n");
        dump_contam_prof(outf, cf.contam_prof, samples);
        fclose(outf);
    }
    
    // Write contamination rate (and standard error) per cell to disk
    {
        string fname = output_prefix + ".contam_rate";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_contam_rates(outf, cf.contam_rate, cf.contam_rate_se, samples,
            barcode_group);
    }
        
    return 0;
}
