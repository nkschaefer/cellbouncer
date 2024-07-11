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
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bam.h>
#include <optimML/mixcomp.h>
#include "common.h"
#include "demux_vcf_hts.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "bulkprops [OPTIONS]\n");
    fprintf(stderr, "Given bulk sequencing data and a VCF of individuals, attempts to infer\n");
    fprintf(stderr, "the proportion of the bulk data composed of each individual.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants. Only biallelic SNPs \n");
    fprintf(stderr, "       will be considered, and phasing will be ignored.\n");
    fprintf(stderr, "    --output_prefix -o Base name for output files. Will create [output_prefix].bulkprops\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --qual -q Minimum variant quality to consider (default 50)\n");
    fprintf(stderr, "    --ids -i If the VCF file contains individuals that you do not\n");
    fprintf(stderr, "       expect to see in your sample, specify individuals to include here.\n");
    fprintf(stderr, "       These identities will be considered, as well as all possible doublet\n");
    fprintf(stderr, "       combinations of them. Expects a text file, with one individual ID per\n");
    fprintf(stderr, "       line, matching individual names in the VCF.\n");
    fprintf(stderr, "----- I/O options -----\n");
    fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"index_jump", no_argument, 0, 'j'},
       {"ids", required_argument, 0, 'i'},
       {"qual", required_argument, 0, 'q'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string vcf_file = "";
    bool stream = true; // Opposite of index-jumping
    string output_prefix = "";
    int vq = 50;
    string idfile;
    bool idfile_given = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:q:i:jh", long_options, &option_index )) != -1){
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
            case 'v':
                vcf_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'i':
                idfile_given = true;
                idfile = optarg;
                break;
            case 'q':
                vq = atoi(optarg);
                break;
            case 'j':
                stream = false;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (vq < 0){
        fprintf(stderr, "ERROR: variant quality must be a positive integer (or 0 for no filter)\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        exit(1);
    }
    
    // Init BAM reader
    bam_reader reader = bam_reader();

    if (bamfile.length() == 0){
        fprintf(stderr, "ERROR: bam file (--bam) required\n");
        exit(1);
    }    
    reader.set_file(bamfile);

    // Store the names of all individuals in the VCF
    vector<string> samples;
    
    // Store data about SNPs
    map<int, map<int, var> > snpdat;
    
    int nsnps;
    
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: vcf file is required\n");
        exit(1);
    }
    fprintf(stderr, "Loading variants from VCF/BCF...\n");
    nsnps = read_vcf(vcf_file, reader, samples, snpdat, vq, false, false);

    set<int> allowed_ids;
    set<int> allowed_ids2;

    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile.c_str());
        }
    }
    
    // Store ref and alt counts at SNPs
    map<int, map<int, pair<float, float> > > snp_ref_alt;
    
    // Store erronenous (non ref-alt) counts at SNPs
    map<int, map<int, float> > snp_err;

    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;

    // initialize the BAM reader
    reader.set_file(bamfile);
    
    // retrieve cell barcodes
    reader.set_cb();
    
    // Get a mapping of chromosome names to internal numeric IDs in the 
    // BAM file. This is necessary to reconcile how HTSLib might represent
    // the same chromosome name in the BAM vs the VCF file.
    map<string, int> chrom2tid = reader.get_seq2tid();
    map<int, string> tid2chrom;
    for (map<string, int>::iterator ct = chrom2tid.begin(); ct != chrom2tid.end(); ++ct){
        tid2chrom.insert(make_pair(ct->second, ct->first));
    } 
    
    int nsnp_processed = 0;
    if (stream){

        // Read through entire BAM file and look for informative SNPs along the way
        // (default setting, appropriate for large numbers of SNPs).
        int curtid = -1;
        
        map<int, var>::iterator cursnp;
        while (reader.next()){
            if (curtid != reader.tid()){
                // Started a new chromosome
                if (curtid != -1){
                    while (cursnp != snpdat[curtid].end()){
                        ++nsnp_processed;
                        ++cursnp;
                    }
                }
                cursnp = snpdat[reader.tid()].begin();
                curtid = reader.tid();
            }
            // Advance to position within cur read
            while (cursnp != snpdat[reader.tid()].end() && 
                cursnp->first < reader.reference_start){
                ++nsnp_processed;
                ++cursnp;
            }
            // Create a second iterator to look ahead for any additional SNPs 
            // within the current read
            map<int, var>::iterator cursnp2 = cursnp;
            while (cursnp2 != snpdat[reader.tid()].end() && 
                cursnp2->first >= reader.reference_start && 
                cursnp2->first <= reader.reference_end){   
                process_bam_record_bulk(reader, cursnp2->first, cursnp2->second,
                    snp_ref_alt, snp_err);
                ++cursnp2;
            }
            if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, nsnps); 
                last_print = nsnp_processed;
            }
        }
        
    }
    else{
        // Visit each SNP and index-jump in the BAM to it.
        for (map<int, map<int, var> >::iterator curchrom = snpdat.begin();
            curchrom != snpdat.end(); ++curchrom){

            int tid = curchrom->first;
            for (map<int, var>::iterator cursnp = curchrom->second.begin();
                cursnp != curchrom->second.end(); ++cursnp){
                int pos = cursnp->first;
                
                reader.set_query_site(tid, pos); 
                
                while(reader.next()){
                    process_bam_record_bulk(reader, cursnp->first, cursnp->second, 
                        snp_ref_alt, snp_err);
                }
                
                ++nsnp_processed;        
            
                if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                    fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, 
                        nsnps); 
                    last_print = nsnp_processed;
                }
            }
        }
    }
    fprintf(stderr, "Processed %d of %d SNPs\n", nsnp_processed, nsnps);
    
    // Calculate error rate
    int nsnp_data = 0;
    for (map<int, map<int, float> >::iterator s = snp_err.begin(); s != snp_err.end(); ++s){
        for (map<int, float>::iterator s2 = s->second.begin(); s2 != s->second.end(); ++s2){
            nsnp_data++;
        }
    }
    double frac = 1.0/(double)nsnp_data;
    double err_rate = 0.0;
    for (map<int, map<int, float> >::iterator s = snp_err.begin(); s != snp_err.end(); ++s){
        for (map<int, float>::iterator s2 = s->second.begin(); s2 != s->second.end(); ++s2){
            double ref = snp_ref_alt[s->first][s2->first].first;
            double alt = snp_ref_alt[s->first][s2->first].second;
            if (ref + alt + s2->second > 0){
                double ef = s2->second / (s2->second + ref + alt);
                err_rate += frac*ef;
            }
        }
    }
    if (err_rate == 0.0 || isnan(err_rate)){
        err_rate = 0.001;
    }
    fprintf(stderr, "Approximate error rate = %f\n", err_rate);

    // Now solve mix prop problem
    vector<double> n;
    vector<double> k;
    vector<vector<double> > expfracs_all;

    for (map<int, map<int, pair<float, float> > >::iterator s = snp_ref_alt.begin(); s != snp_ref_alt.end();
        ++s){
        for (map<int, pair<float, float> >::iterator s2 = s->second.begin(); s2 != s->second.end(); ++s2){
            
            bool miss = false;
            vector<double> expfracs;

            for (int i = 0; i < samples.size(); ++i){
                if (snpdat[s->first][s2->first].haps_covered[i]){
                    int nalt = 0;
                    if (snpdat[s->first][s2->first].haps1[i]){
                        nalt++;
                    }
                    if (snpdat[s->first][s2->first].haps2[i]){
                        nalt++;
                    }
                    double expfrac = (double)nalt/2.0;
                    // Account for error rate
                    if (expfrac == 0){
                        expfrac += err_rate;
                    }
                    else if (expfrac == 1.0){
                        expfrac -= err_rate;
                    }
                    expfracs.push_back(expfrac);
                }
                else{
                    miss = true;
                    break;
                }
            }

            if (!miss){
                n.push_back(s2->second.first + s2->second.second);
                k.push_back(s2->second.second);
                expfracs_all.push_back(expfracs);          
            } 
        }
    }
    
    optimML::mixcomp_solver solver(expfracs_all, "binom", n, k);
    solver.solve();
    
    string outfn = output_prefix + ".bulkprops";
    FILE* outf = fopen(outfn.c_str(), "w");
    for (int i = 0; i < solver.results.size(); ++i){
        fprintf(outf, "%s\t%f\n", samples[i].c_str(), solver.results[i]);
    }
    fclose(outf);

    return 0;
}
