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
#include <random>
#include <functional>
#include <utility>
#include <math.h>
#include <float.h>
#include <htslib/vcf.h>
#include <htswrapper/gzreader.h>
#include <zlib.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "filter_vcf_species_fixed [OPTIONS]\n");
    fprintf(stderr, "Given a VCF file containing variatns called across multiple species,\n");
    fprintf(stderr, "along with a file mapping individuals to species, finds sites at which\n");
    fprintf(stderr, "all species are fixed for different alleles, and outputs a subsampled\n");
    fprintf(stderr, "version of the VCF file at those sites.\n");
    fprintf(stderr, "If you do not specify the species an individual belongs to, it will not be\n");
    fprintf(stderr, "used to determine fixed sites, but will still be included in the output.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants.\n");
    fprintf(stderr, "    --species -s A file mapping individual ID to species, one per line,\n");
    fprintf(stderr, "      tab-separated.\n");
    fprintf(stderr, "    --output -o Output file name.\n");
    fprintf(stderr, "    --output_fmt -O Output file format (from HTSLib: v = VCF; z = gzVCF;\n");
    fprintf(stderr, "      b = BCF; etc. Default: z\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

/**
 * Returns true if fixed between species
 * false if not fixed between species
 */
bool proc_bcf_record(bcf1_t* bcf_record,
    bcf_hdr_t* bcf_header,
    int num_samples,
    int num_species,
    vector<int>& sample_species){

    vector<int> species_ref(num_species, 0);
    vector<int> species_alt(num_species, 0);

    // Load ref/alt alleles and other stuff
    // This puts alleles in bcf_record->d.allele[index]
    // Options for parameter 2:
    
    // BCF_UN_STR  1       // up to ALT inclusive
    // BCF_UN_FLT  2       // up to FILTER
    // BCF_UN_INFO 4       // up to INFO
    // BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
    // BCF_UN_FMT  8                           // unpack format and each sample
    // BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
    // BCF_UN_ALL (BCF_UN_SHR|BCF_UN_FMT) // everything
    
    bcf_unpack(bcf_record, BCF_UN_STR);
    bool pass = true;
    // pass = biallelic, no indels, ref/alt both A/C/G/T
    for (int i = 0; i < bcf_record->n_allele; ++i){
        if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
            strcmp(bcf_record->d.allele[i], "C") != 0 &&
            strcmp(bcf_record->d.allele[i], "G") != 0 && 
            strcmp(bcf_record->d.allele[i], "T") != 0){
            pass = false;
            break;
        }
    }
    if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
        pass = false;
    }
    if (pass){
        // Get all available genotypes.
        int32_t* gts = NULL;
        int n_gts = 0;
        int nmiss = 0;
        int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
        if (num_loaded <= 0){
            fprintf(stderr, "ERROR loading genotypes at %s %ld\n", 
                bcf_hdr_id2name(bcf_header, bcf_record->rid), (long int) bcf_record->pos);
            exit(1);
        }
        
        // Assume ploidy = 2
        int ploidy = 2;
        //int ploidy = n_gts / num_samples; 
       
        int nref = 0;
        int nalt = 0;
         
        for (int i = 0; i < num_samples; ++i){
            int32_t* gtptr = gts + i*ploidy;
            int species_idx = sample_species[i];
            if (species_idx >= 0){   
                if (!bcf_gt_is_missing(gtptr[0])){
                    if (bcf_gt_allele(gtptr[0]) == 0){
                        nref++;        
                        species_ref[species_idx]++;
                    }
                    else{
                        nalt++;
                        species_alt[species_idx]++;
                    }
                    if (bcf_gt_allele(gtptr[1]) == 0){
                        nref++;
                        species_ref[species_idx]++;
                    }
                    else{
                        nalt++;
                        species_alt[species_idx]++;
                    }
                }
                if (species_ref[species_idx] > 0 && species_alt[species_idx] > 0){
                    return false;
                }
            }
        }
        // Exclude sites without polymorphism
        if (nref == 0 || nalt == 0){
            return false;
        } 
        // At this point, there is no within-species polymorphism, and there is polymorphism.
        // Therefore, this must be a fixed difference across species.
        return true;
    }
    // Site did not pass filters.
    return false;
}

void parse_species(string speciesfile, map<string, string>& id2species){
    gzreader reader(speciesfile);
    reader.delimited('\t');
    /*
    inf = ifstream(speciesfile);
    string id;
    string species;
    while (inf >> id >> species){
    */
    while(reader.next()){
        string id = reader.fields[0];
        string species = reader.fields[1];
        if (id2species.count(id) > 0 && id2species[id] != species){
            fprintf(stderr, "ERROR: id %s mapped to multiple species in %s\n", 
                id.c_str(), speciesfile.c_str());
            exit(1);
        }
        id2species.insert(make_pair(id, species));
    }
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"species", required_argument, 0, 's'},
       {"output", required_argument, 0, 'o'},
       {"output_fmt", required_argument, 0, 'O'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string outfile = "";
    string speciesfile = "";
    string outfmt = "z";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
   while((ch = getopt_long(argc, argv, "v:o:s:O:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'O':
                outfmt = optarg;
                break;
            case 's':
                speciesfile = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    if (speciesfile == ""){
        fprintf(stderr, "ERROR: species file required\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: output file name required\n");
        exit(1);
    }

    if (outfmt == "z"){
        // Clean up output file name
        if (outfile.rfind(".vcf") == string::npos){
            outfile += ".vcf.gz";
        }
        else if (outfile.rfind(".gz") == string::npos){
            outfile += ".gz";
        }
    }
    
    // Read individual -> species assignments
    map<string, string> id2species;
    parse_species(speciesfile, id2species);
    
    // Store the names of all individuals in the VCF
    vector<string> samples;
    vector<string> samples_species;

    // First pass: count occurrences of each branch, only considering sites without
    // missing genotypes.
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
        if (id2species.count(bcf_header->samples[i]) == 0){
            fprintf(stderr, "WARNING: sample %s not assigned to a species\n", bcf_header->samples[i]);
            fprintf(stderr, "  This individual will be included in output but will not be used to determine \
which SNPs are fixed species differences.\n");
        }
    }
    
    // Write gz-compressed by default
    string writestr = "w" + outfmt;
    htsFile* outf = hts_open(outfile.c_str(), writestr.c_str());
    int write_success = bcf_hdr_write(outf, bcf_header);

    // Convert sample names to species indices
    vector<int> sample_species;
    map<string, int> species2idx;
    int species_idx = 0;
    for (vector<string>::iterator samp = samples.begin(); samp != samples.end(); ++samp){
        if (id2species.count(*samp) == 0){
            sample_species.push_back(-1);
        }
        else{
            string spec = id2species[*samp];
            int idx;
            if (species2idx.count(spec) > 0){
                idx = species2idx[spec];
            }
            else{
                idx = species_idx;
                species2idx.insert(make_pair(spec, idx));
                species_idx++;
            }
            sample_species.push_back(idx);
        }
    }
    int num_species = species_idx;
    
    long int kept = 0;
    
    long int nsnp;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){ 
            bool pass;        
            
            if (proc_bcf_record(bcf_record, bcf_header, num_samples,
                num_species, sample_species)){
                
                // This is a fixed difference.
                int ret = bcf_write1(outf, bcf_header, bcf_record);
                kept++;        
            }
        }
        ++nsnp;
        if (nsnp % 1000 == 0){
            fprintf(stderr, "Processed %ld SNPs\r", nsnp);
        }
    }
    
    hts_close(bcf_reader);
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(outf);
    
    fprintf(stderr, "Processed %ld SNPs\n", nsnp);
    fprintf(stderr, "Kept %ld of %ld SNPs (%.2f%%)\n", kept, nsnp, 100.0*(double)kept/(double)nsnp);

    return 0;
}
