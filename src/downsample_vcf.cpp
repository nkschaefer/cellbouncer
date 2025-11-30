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
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <optimML/brent.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <iomanip>
#include "common.h"
#include "demux_vcf_io.h"
#include "demux_vcf_hts.h"
#include "downsample_vcf.h"
#include <chrono>

using std::cout;
using std::endl;
using namespace std;


bool operator<(const std::bitset<NBITS>& a, const std::bitset<NBITS>& b) {
    for (int i = NBITS-1; i >= 0; --i) {
        if (a[i] != b[i]) return b[i]; // compare from MSB down
    }
    return false;
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "downsample_vcf [OPTIONS]\n");
    fprintf(stderr, "Given a VCF file containing variants segregating within a panel of\n");
    fprintf(stderr, "deeply divergent genomes (e.g. multiple species on the same reference\n");
    fprintf(stderr, "genome), intelligently downsamples the variant panel to target a desired\n");
    fprintf(stderr, "final number of variants.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This is useful because in these cases, many SNPs can reflect fixed differences\n");
    fprintf(stderr, "between groups; for the sake of demultiplexing, it is not necessary to track\n");
    fprintf(stderr, "all such variants. This program counts SNPs defining clade membership. It then\n");
    fprintf(stderr, "finds a strategy to reach the target output SNP count by letting low-count clades\n");
    fprintf(stderr, "pass through and downsampling high-count clades geometrically (so the highest-count\n");
    fprintf(stderr, "clades are more dramatically downsampled than lower-count clades).\n");
    fprintf(stderr, "This should result in the targeted number of variants, with extremely long\n");
    fprintf(stderr, "branches (e.g. orangutan vs the human/chimp/bonobo clade) downsampled more\n");
    fprintf(stderr, "than shorter branches.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "To save memory, this program currently treats both heterozygous and homozygous alt\n");
    fprintf(stderr, "alleles the same - e.g. two individuals are related at a site if either of them have\n");
    fprintf(stderr, "the same allele at any frequency.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Writes a new gzipped VCF to the specified output file.\n");
    fprintf(stderr, "Writes information about downsampled clades to stdout. This only includes clades that\n");
    fprintf(stderr, "were downsampled.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants.\n");
    fprintf(stderr, "    --num -n Desired final number of SNPs.\n");
    fprintf(stderr, "    --output -o Output VCF file name. Will be gzipped.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

/**
 * Returns true if the site can be used
 * Returns false if the site cannot be used
 */
bool proc_bcf_record(bcf1_t* bcf_record,
    bcf_hdr_t* bcf_header,
    int num_samples,
    bitset<NBITS>& alt,
    //bitset<NBITS>& alt2,
    bitset<NBITS>& alt_flip,
    //bitset<NBITS>& alt2_flip,
    bitset<NBITS>& present,
    bool& pass){

    alt.reset();
    alt_flip.reset();
    present.reset();

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
    pass = true;
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
         
        //bool has_het = false;
        for (int i = 0; i < num_samples; ++i){
            int32_t* gtptr = gts + i*ploidy;
            
            if (!bcf_gt_is_missing(gtptr[0])){
                present.set(i);
                
                alt_flip.set(i);
                
                if (bcf_gt_allele(gtptr[0]) != 0 || bcf_gt_allele(gtptr[1]) != 0){
                    // This individual has at least one alt allele copy.
                    alt.set(i);
                    alt_flip.reset(i);
                    
                    nalt++;
                } 
                else{
                    nref++;
                }
                /*
                if (bcf_gt_allele(gtptr[0]) != 0){
                    nalt++;
                    
                    // alt2 tracks het & hom alt sites.
                    alt2.set(i);
                    alt2_flip.reset(i);

                    if (bcf_gt_allele(gtptr[1]) != 0){
                        nalt++;
                        // hom alt
                        // alt tracks hom alt sites.
                        alt.set(i);
                        alt_flip.reset(i);
                    }
                    else{
                        nref++;
                    }
                }
                else{
                    nref++;
                    if (bcf_gt_allele(gtptr[1]) != 0){
                        nalt++;
                        // het
                        alt2.set(i);
                        alt2_flip.reset(i);
                    }
                    else{
                        nref++;
                    }
                }
                */
            }    
        }
     
        //alt = alt2;

        // Exclude sites without polymorphism
        if (nref == 0 || nalt == 0){
            pass = false;
            return false;
        } 
        return true;
    }
    return false;
}

void count_branch(unordered_map<bitset<NBITS>, double>& branchcounts,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int cladecount,
    int cladecount_flip,
    double count){
   
    if (cladecount < cladecount_flip || (cladecount == cladecount_flip && 
        clade < clade_flip)){
        if (branchcounts.count(clade) == 0){
            branchcounts.insert(make_pair(clade, count));
        }
        else{
            branchcounts[clade] += count;
        }
    }
    else{
        if (branchcounts.count(clade_flip) == 0){
            branchcounts.insert(make_pair(clade_flip, count));
        }
        else{
            branchcounts[clade_flip] += count;
        }
    }
}

void count_branch_missing(unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double>& branchcounts,
    unordered_map<bitset<NBITS>, bitset<NBITS> >& miss2flip,
    bitset<NBITS>& present,
    bitset<NBITS>& clade,
    bitset<NBITS>& clade_flip,
    int cladecount,
    int cladecount_flip,
    double count){
    
    pair<bitset<NBITS>, bitset<NBITS> > key;
    
    if (cladecount < cladecount_flip || (cladecount == cladecount_flip && 
       clade < clade_flip)){
        key = make_pair(present, clade);
        if (miss2flip.count(clade) == 0){
            miss2flip.insert(make_pair(clade, clade_flip));
        }
    }
    else{
        key = make_pair(present, clade_flip);
        if (miss2flip.count(clade_flip) == 0){
            miss2flip.insert(make_pair(clade_flip, clade));
        }
    }
    if (branchcounts.count(key) == 0){
        branchcounts.insert(make_pair(key, count));
    }
    else{
        branchcounts[key] += count;
    }
}

/**
 * Function for Brent's method root finding
 * Finds x that allows sum_i(a_i^x) = n where a is the vector of clade counts
 *  (and i is index into that vector), and n is the desired number of SNPs in
 *  output file
 */
double brentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    double sum_target = (double)data_i.at("target");
    char buf[50];
    string bufstr;
    double sum = 0.0;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        sum += pow(count, param);
    }
    return sum - sum_target;
}

/**
 * Derivative of function for Brent's method root finding
 */
double dbrentfun(double param, const map<string, double>& data_d, const map<string, int>& data_i){
    int num = data_i.at("num");
    double sum_target = (double)data_i.at("target");
    char buf[50];
    string bufstr;
    double d_df = 0.0;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "x%d", i);
        bufstr = buf;
        double count = data_d.at(bufstr);
        d_df += pow(count, param) * log(count);
    }
    return d_df;
}


int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"output", required_argument, 0, 'o'},
       {"num", required_argument, 0, 'n'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string outfile = "";
    int num = -1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:n:o:h", long_options, &option_index )) != -1){
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
            case 'n':
                num = atoi(optarg);
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
    if (num <= 0){
        fprintf(stderr, "ERROR: --num/-n must be positive (Recommend > 1M)\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: --output/-o is required.\n");
        exit(1);
    }
    // Clean up output file name
    if (outfile.rfind(".vcf") == string::npos){
        outfile += ".vcf.gz";
    }
    else if (outfile.rfind(".gz") == string::npos){
        outfile += ".gz";
    }

    // Count branches
    unordered_map<bitset<NBITS>, double> branchcounts;
    unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double> branchcounts_missing;
    
    //unordered_map<bitset<NBITS>, int> branchcounts_het;

    // Store the names of all individuals in the VCF
    vector<string> samples;
    
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
    }
    
    if (samples.size() * 2 > NBITS){
        fprintf(stderr, "ERROR: too many samples in VCF. Please recompile with NBITS = %ld or higher.\n",
            samples.size() * 2);
        exit(1); 
    }

    int nsnp = 0;
    
    bitset<NBITS> alt;
    //bitset<NBITS> alt2;
    bitset<NBITS> alt_flip;
    //bitset<NBITS> alt2_flip;
    bitset<NBITS> present;
    
    unordered_map<bitset<NBITS>, bitset<NBITS> > miss2flip;

    fprintf(stderr, "Counting mutations on branches...\n");
    
    //long int het_site_count = 0;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){ 
            bool pass;        
            if (proc_bcf_record(bcf_record, bcf_header, num_samples,
                alt, alt_flip, present, pass)){  
                
                int ac = alt.count();
                int afc = alt_flip.count();
            
                if (present.count() == num_samples){
                    // Only need to track alt (no hets)
                    count_branch(branchcounts, alt, alt_flip, ac, afc, 1.0);      
                }
                else{
                    // Missing some individuals.
                    pair<bitset<NBITS>, bitset<NBITS> > key;
                    // No need to track hets.
                    count_branch_missing(branchcounts_missing,
                        miss2flip, present, alt, alt_flip, ac, afc, 1.0);      
                }
            }
        }
        ++nsnp;
        if (nsnp % 1000 == 0){
            fprintf(stderr, "Processed %d SNPs\r", nsnp);
        }
    }
    
    hts_close(bcf_reader);
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);

    fprintf(stderr, "Processed %d SNPs\n", nsnp);
    nsnp = 0;
    
    if (branchcounts.size() >= num){
        fprintf(stderr, "ERROR: %ld distinct allele sharing patterns found, but %d SNPs requested.\n", 
            branchcounts.size(), num);
        fprintf(stderr, "  Please set sample size to at least %ld.\n", branchcounts.size());
        exit(1);
    }

    map<int, vector<bitset<NBITS> > > size2clades;
    for (unordered_map<bitset<NBITS>, double>::iterator bc = branchcounts.begin(); bc != 
        branchcounts.end(); ++bc){
        for (int i = 1; i <= bc->first.count(); ++i){
            if (size2clades.count(i) == 0){
                vector<bitset<NBITS> > v;
                size2clades.insert(make_pair(i, v));
            }
            size2clades[i].push_back(bc->first);
        }
    }
    
    fprintf(stderr, "Adjust counts using data from sites with missing genotypes...\n");
    unordered_map<bitset<NBITS>, vector<bitset<NBITS> > > keyuniqparent;
    
    unordered_map<bitset<NBITS>, float> to_add;
    // This will eventually store the probability that the key ((mask, clade)) representing a site type
    // with missing genotypes actually comes from the next map's key (clade) representing a site type
    // with no missing genotypes
    unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, unordered_map<bitset<NBITS>, float> > miss_clade_probs;
    // Determine how to assign clades with missing members to parent clades
    for (unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double>::iterator bcm = branchcounts_missing.begin();
        bcm != branchcounts_missing.end(); ++bcm){
        
        unordered_map<bitset<NBITS>, float> mx;
        miss_clade_probs.insert(make_pair(bcm->first, mx));
        double parent_tot = 0.0;

        int k = bcm->first.second.count();
        if (keyuniqparent.count(bcm->first.second) > 0){
            // Use cached.
            for (vector<bitset<NBITS> >::iterator v = keyuniqparent[bcm->first.second].begin();
                v != keyuniqparent[bcm->first.second].end(); ++v){
                if ((*v & bcm->first.first).count() == k){
                    // Candidate clade
                    double cladecount = branchcounts[*v];
                    parent_tot += cladecount;
                    miss_clade_probs[bcm->first].insert(make_pair(*v, (float)cladecount));
                }
            }
        }
        else{
            for (vector<bitset<NBITS> >::iterator v = size2clades[k].begin();
                v != size2clades[k].end(); ++v){
                if ((*v & bcm->first.second).count() == k){
                    // Parent
                    if (keyuniqparent.count(bcm->first.second) == 0){
                        vector<bitset<NBITS> > vec;
                        keyuniqparent.insert(make_pair(bcm->first.second, vec));
                    }
                    keyuniqparent[bcm->first.second].push_back(*v);
                    if ((*v & bcm->first.first).count() == k){
                        // Candidate clade
                        double cladecount = branchcounts[*v];
                        parent_tot += cladecount;
                        miss_clade_probs[bcm->first].insert(make_pair(*v, (float)cladecount));
                    }
                }
            }
        }

        // Also look up flipped version
        if (miss2flip.count(bcm->first.second) > 0){
            bitset<NBITS> flipped = (miss2flip[bcm->first.second] & bcm->first.first);
            k = flipped.count();
            if (keyuniqparent.count(flipped) > 0){
                // use cached
                for (vector<bitset<NBITS> >::iterator v = keyuniqparent[flipped].begin();
                    v != keyuniqparent[flipped].end(); ++v){
                    if ((*v & bcm->first.first).count() == k){
                        // Candidate clade
                        double cladecount = branchcounts[*v];
                        parent_tot += cladecount;
                        miss_clade_probs[bcm->first].insert(make_pair(*v, (float)cladecount));
                    }
                }
            }
            else{
                 for (vector<bitset<NBITS> >::iterator v = size2clades[k].begin();
                    v != size2clades[k].end(); ++v){
                    if ((*v & flipped).count() == k){
                        // Parent
                        if (keyuniqparent.count(flipped) == 0){
                            vector<bitset<NBITS> > vec;
                            keyuniqparent.insert(make_pair(flipped, vec));
                        }
                        keyuniqparent[flipped].push_back(*v);
                        if ((*v & bcm->first.first).count() == k){
                            // Candidate clade
                            double cladecount = branchcounts[*v];
                            parent_tot += cladecount;
                            miss_clade_probs[bcm->first].insert(make_pair(*v, (float)cladecount));
                        }
                    }
                }   
           } 
        }
        
        // Visit each no-missing-genotypes clade that might be the true state of a clade
        // with missing genotypes
        for (unordered_map<bitset<NBITS>, float>::iterator x = miss_clade_probs[bcm->first].begin(); 
            x != miss_clade_probs[bcm->first].end(); ++x){
            x->second /= (double)parent_tot;
            if (to_add.count(x->first) == 0){
                to_add.insert(make_pair(x->first, 0.0));
            }
            to_add[x->first] += (x->second * (float)bcm->second);
        }
    }

    // Determine number of each clade to subsample.
    unordered_map<bitset<NBITS>, double> downsample_prob;
    
    fprintf(stderr, "Determine number of sites to downsample...\n");

    double allbc2 = 0; 
    vector<bitset<NBITS> > bs_idx;
    vector<pair<double, int> > clcountsort;
    
    for (unordered_map<bitset<NBITS>, double>::iterator bc = branchcounts.begin(); 
        bc != branchcounts.end(); ){
        if (to_add.count(bc->first) > 0){
            bc->second += to_add[bc->first];
            to_add.erase(bc->first);
        }
        allbc2 += bc->second;
        clcountsort.push_back(make_pair(-bc->second, bs_idx.size()));
        bs_idx.push_back(bc->first);
        branchcounts.erase(bc++);
    }
    
    map<string, double> sampstr2prob;
    if (allbc2 > (double)num){
        
        sort(clcountsort.begin(), clcountsort.end());
        
        for (int n_hi = 1; n_hi <= clcountsort.size(); ++n_hi){ 
            double sum_hi = 0;
            
            optimML::brent_solver solver(brentfun, dbrentfun);
            solver.set_root();
            //solver.constrain_01();
            char buf[50];
            vector<vector<double> > vecs;

            for (int i = 0; i < n_hi; ++i){
                vecs.push_back(vector<double>{ -clcountsort[i].first });    
                sprintf(&buf[0], "x%d", i);
                string bufstr = buf;
                bool success = solver.add_data(bufstr, vecs[vecs.size()-1]);
                sum_hi += -clcountsort[i].first;
            }

            double sum_lo = allbc2 - sum_hi;
            if (n_hi == clcountsort.size() || sum_lo < (double)num ){
                double lastcount = -clcountsort[n_hi-1].first;
                double nextcount = -1;
                if (n_hi < clcountsort.size()){
                    nextcount = -clcountsort[n_hi].first;
                }
                
                solver.add_data_fixed("target", (int)round((double)num - sum_lo));
                solver.add_data_fixed("num", n_hi);
                double res = solver.solve(0,1);

                if (res < 1.0){
                    lastcount = pow(lastcount, res);
                    
                    double testsum = 0.0;
                    for (int i = 0; i < n_hi; ++i){
                        testsum += pow(-clcountsort[i].first, res);
                    } 
                    if (lastcount > nextcount){
                        fprintf(stdout, "group1\tgroup2\tcount\tfrac_keep\n");
                        for (int i = 0 ; i < n_hi; ++i){
                            double dsp = pow(-clcountsort[i].first, res) / -clcountsort[i].first;
                            downsample_prob.insert(make_pair(bs_idx[clcountsort[i].second], dsp));
                            set<string> samp1;
                            set<string> samp2;
                            for (int x = 0; x < samples.size(); ++x){
                                if (bs_idx[clcountsort[i].second].test(x)){
                                    samp1.insert(samples[x]);
                                }
                                else{
                                    samp2.insert(samples[x]);
                                }
                            }
                            string sampstr1 = "";
                            string sampstr2 = "";
                            bool first = true;
                            for (set<string>::iterator s = samp1.begin(); s != samp1.end(); ++s){
                                if (!first){
                                    sampstr1 += ",";
                                }
                                first = false;
                                sampstr1 += *s;
                            }
                            first = true;
                            for (set<string>::iterator s = samp2.begin(); s != samp2.end(); ++s){
                                if (!first){
                                    sampstr2 += ",";
                                }
                                first = false;
                                sampstr2 += *s;
                            }
                            fprintf(stdout, "%s\t%s\t%d\t%f\n", sampstr1.c_str(), sampstr2.c_str(),
                                (int)round(-clcountsort[i].first), dsp);
                        }
                        break;
                    }
                }
            }
        }
    }
    else{
        fprintf(stderr, "%f total clade counts; no downsampling needed\n", allbc2);
        return 0;
    }
    
    // Determine probability of keeping clades with missing genotypes
    unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, double> downsample_miss_prob;

    // Eliminate any missing -> full probability mapping for clades that
    // are not going to be downsampled
    for (unordered_map<pair<bitset<NBITS>, bitset<NBITS> >, unordered_map<bitset<NBITS>, float> >::iterator mcp = 
        miss_clade_probs.begin(); mcp != miss_clade_probs.end(); ){
        // Probability of keeping
        double p = 0.0;
        for (unordered_map<bitset<NBITS>, float>::iterator mcp2 = mcp->second.begin(); 
            mcp2 != mcp->second.end(); ){
            if (downsample_prob.count(mcp2->first) > 0){
                p += mcp2->second * downsample_prob[mcp2->first];
            }
            else{
                 p += mcp2->second;
            }
            mcp->second.erase(mcp2++);
        }
        if (p < 1.0){
            downsample_miss_prob.insert(make_pair(mcp->first, p));
        }
        miss_clade_probs.erase(mcp++);
    } 
    if (downsample_prob.size() == 0){
        fprintf(stderr, "No clades were found with more than %d occurrences.\n", num);
        fprintf(stderr, "There is nothing to downsample.\n");
        return 0;
    }
    // Go back through VCF and downsample as necessary.
    int n_rm = 0;
    
    // Re-open input
    bcf_reader = hts_open(vcf_file.c_str(), "r");
    bcf_header = bcf_hdr_read(bcf_reader);
    bcf_record = bcf_init();
    
    // Open output VCF (stdout)
    // Write gz-compressed by default
    htsFile* outf = hts_open(outfile.c_str(), "wz");
    int write_success = bcf_hdr_write(outf, bcf_header);

    // Generate random numbers between 0 and 1
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> rand_dist(0.0, 1.0);
    
    int elim = 0;
    
    fprintf(stderr, "Filtering VCF...\n");

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){ 
            bool pass;
            if (proc_bcf_record(bcf_record, bcf_header, num_samples,
                alt, alt_flip, present, pass)){  
                
                double samp_prob = 1.0;
                int ac = alt.count();
                int afc = alt_flip.count();
                
                if (present.count() == num_samples){
                    // No hets to track.
                    if (ac < afc || (ac == afc && alt < alt_flip)){
                        if (downsample_prob.count(alt) > 0){
                            samp_prob = downsample_prob[alt];
                        }
                    }
                    else{
                        if (downsample_prob.count(alt_flip) > 0){
                            samp_prob = downsample_prob[alt_flip];
                        }
                    }
                    /*
                    else{
                        double samp_prob1 = 1.0;
                        double samp_prob2 = 1.0;
                        if (ac < afc || (ac == afc && alt < alt_flip)){
                            if (downsample_prob.count(alt) > 0){
                                samp_prob1 = downsample_prob[alt];
                            }
                        }
                        else{
                            if (downsample_prob.count(alt_flip) > 0){
                                samp_prob1 = downsample_prob[alt_flip];
                            }
                        }
                        if (a2c < a2fc || (a2c == a2fc && alt2 < alt2_flip)){
                            if (downsample_prob.count(alt2) > 0){
                                samp_prob2 = downsample_prob[alt2];
                            }
                        }
                        else{
                            if (downsample_prob.count(alt2_flip) > 0){
                                samp_prob2 = downsample_prob[alt2_flip];
                            }
                        }
                        samp_prob = 0.5 * samp_prob1 + 0.5 * samp_prob2;
                    }
                    */
                }
                else{
                    pair<bitset<NBITS>, bitset<NBITS> > key;

                    // No hets to track.
                    if (ac < afc || (ac == afc && alt < alt_flip)){
                        key = make_pair(present, alt);
                    }
                    else{
                        key = make_pair(present, alt_flip);
                    }
                    if (downsample_miss_prob.count(key) > 0){
                        samp_prob = downsample_miss_prob[key];
                    }
                    /*
                    else{
                        double samp_prob1 = 1.0;
                        double samp_prob2 = 1.0;
                        if (ac < afc || (ac == afc && alt < alt_flip)){
                            key = make_pair(present, alt);
                        }
                        else{
                            key = make_pair(present, alt_flip);
                        }
                        if (downsample_miss_prob.count(key) > 0){
                            samp_prob1 = downsample_miss_prob[key];
                        }
                        if (a2c < a2fc || (a2c == a2fc && alt2 < alt2_flip)){
                            key = make_pair(present, alt2);
                        }
                        else{
                            key = make_pair(present, alt2_flip);
                        }
                        if (downsample_miss_prob.count(key) > 0){
                            samp_prob2 = downsample_miss_prob[key];
                        }
                        samp_prob = 0.5*samp_prob1 + 0.5*samp_prob2;
                    }
                    */
                }

                bool print_rec = true;
                if (samp_prob < 1.0){
                    double r = rand_dist(gen);
                    if (r <= samp_prob){
                        print_rec = true;
                    }
                    else{
                        print_rec = false;
                    }
                }
                if (print_rec){
                    int ret = bcf_write1(outf, bcf_header, bcf_record);
                }
                else{
                    elim++;
                }
            }
            else{
                // Did not pass (filtering reasons).
                ++elim;
            }
        }
        ++nsnp;
        if (nsnp % 1000 == 0){
            fprintf(stderr, "Processed %d SNPs\r", nsnp);
        }
    }
    
    fprintf(stderr, "Processed %d SNPs\n", nsnp);
    fprintf(stderr, "Eliminated %d SNPs\n", elim);

    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    hts_close(outf);

    return 0;
}
