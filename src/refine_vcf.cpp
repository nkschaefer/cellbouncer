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
#include <float.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/gzreader.h>
#include <mixtureDist/functions.h>
#include <optimML/multivar_ml.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_io.h"
#include "demux_vcf_hts.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "refine_vcf [OPTIONS]\n");
    fprintf(stderr, "Given a BAM file of aligned single-cell sequencing data, a file\n");
    fprintf(stderr, "assigning cell barcodes to individual identities, and a VCF file\n");
    fprintf(stderr, "containing genotypes of individuals, refines the genotypes in the\n");
    fprintf(stderr, "VCF to maximize their likelihood given the cell-individual assignments.\n");
    fprintf(stderr, "Prints a new VCF to stdout.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants. Only biallelic SNPs \n");
    fprintf(stderr, "       will be considered, and phasing will be ignored.\n");
    fprintf(stderr, "    --assignments -a The .assignments file from a CellBouncer program\n");
    fprintf(stderr, "       like demux_vcf\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

void parse_assignments(string& filename, 
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    vector<string>& samples){
    
    map<string, int> samp2idx;
    for (int i = 0; i < samples.size(); ++i){
        samp2idx.insert(make_pair(samples[i], i));
    }

    ifstream inf(filename);
    string bc_str;
    string id;
    char type;
    double llr;

    while (inf >> bc_str >> id >> type >> llr){
        
        unsigned long bckey = bc_ul(bc_str);
        if (type == 'S'){
            if (samp2idx.count(id) == 0){
                fprintf(stderr, "ERROR: id %s not found in VCF\n", id.c_str());
                exit(1);
            }
            int id_i = samp2idx[id];
            assn.emplace(bckey, id_i);
            assn_llr.emplace(bckey, llr);
        }
        else{
            size_t splitpos = id.find("+");
            if (splitpos != string::npos){
                string id1 = id.substr(0, splitpos);
                string id2 = id.substr(splitpos+1, id.length()-splitpos-1);
                if (samp2idx.count(id1) == 0){
                    fprintf(stderr, "ERROR: individual %s in assignments but not VCF\n", id1.c_str());
                    exit(1);
                }
                if (samp2idx.count(id2) == 0){
                    fprintf(stderr, "ERROR: individual %s in assignments but notn VCF\n", id2.c_str());
                    exit(1);
                }
                int id1_i = samp2idx[id1];
                int id2_i = samp2idx[id2];
                int id3_i;
                if (id1_i < id2_i){
                    id3_i = hap_comb_to_idx(id1_i, id2_i, samples.size());
                }
                else if (id1_i == id2_i){
                    fprintf(stderr, "ERROR parsing doublet ID %s\n", id.c_str());
                    exit(1);
                }
                else{
                    id3_i = hap_comb_to_idx(id2_i, id1_i, samples.size());
                }
                assn.emplace(bckey, id3_i);
                assn_llr.emplace(bckey, llr);
            }
            else{
                fprintf(stderr, "ERROR parsing doublet ID %s\n", id.c_str());
                exit(1);
            }
        }
    }
}

double ll_regt(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    double p;
    double target;
    if (idx2 == -1){
        p = params[idx1];
        //p = round(p/0.5)*0.5;
        target = round(p/0.5)*0.5;
    }
    else{
        //double p1 = round(params[idx1]/0.5)*0.5;
        //double p2 = round(params[idx2]/0.5)*0.5;
        //p = 0.5*p1 + 0.5*p2;
        p = 0.5*params[idx1] + 0.5*params[idx2];
        target = 0.5*round(params[idx1]/0.5)*0.5 + 0.5*round(params[idx2]/0.5)*0.5;
    }

    /*
    if (p < 0.001){
        p = 0.001;
    }
    else if (p > 0.999){
        p = 0.999;
    }
    */
    if (p < 1e-8){
        p = 1e-8;
    } 
    else if ( p > 1.0-1e-8){
        p = 1.0-1e-8;
    }
    return logbinom(n,k,p);
    //1000*abs(p-target);
}

void dll_regt(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    int idx1 = data_i.at("idx1");
    int idx2 = data_i.at("idx2");
    
    double p;
    double target;
    if (idx2 == -1){
        p = params[idx1];
        //p = round(p/0.5)*0.5;
        //pen = penalty(p);
        target = round(p/0.5)*0.5;
    }
    else{
        //double p1 = round(params[idx1]/0.5)*0.5;
        //double p2 = round(params[idx2]/0.5)*0.5;
        //p = 0.5*p1 + 0.5*p2;
        p = 0.5*params[idx1] + 0.5*params[idx2];
        //pen = 0.5*penalty(params[idx1]) + 0.5*penalty(params[idx2]);
        target = 0.5*round(p/0.5)*0.5 + 0.5*round(p/0.5)*0.5;
    }
    
    if (p < 1e-8){
        p = 1e-8;
    } 
    else if ( p > 1.0-1e-8){
        p = 1.0-1e-8;
    }
    /*
    if (p < 0.001){
        p = 0.001;
    }
    else if (p > 0.999){
        p = 0.999;
    }
    */

    double dy_dp = (k-n*p)/(p - p*p);
    if (p < target){
        // Increase p = closer = higher LL
        //dy_dp += 100;
    }
    else{
        //dy_dp -= 100;
    }
    if (idx2 == -1){
        results[idx1] += dy_dp;
    }
    else{
        results[idx1] += 0.5*dy_dp;
        results[idx2] += 0.5*dy_dp;
    }
}


/**
 * Re-genotype a SNP and write results to stdout in VCF format.
 */
bool regt_snp(const string& chrom,
    int pos,
    var& v,
    int n_samples,
    map<int, pair<float, float> >& id_dat,
    bool& removed,
    map<int, double>& weights){
    
    // Solve for maximum likelihood alt allele fractions, and then
    // convert these to discrete genotypes
    vector<double> params;
    vector<int> samp_inds;
    map<int, int> samp2param;
    vector<string> gts;

    map<int, double> newalt;
    set<int> gts_pass;
    
    map<int, float> ind_tots;
    map<int, float> ind_alts;
    map<int, int> n_doub_seen;
    
    for (map<int, pair<float, float> >::iterator x = id_dat.begin(); x != id_dat.end(); ++x){
        if (x->second.first + x->second.second > 0){
            if (x->first >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(x->first, n_samples);
                if (ind_tots.count(comb.first) == 0){
                    ind_tots.insert(make_pair(comb.first, x->second.first + x->second.second));
                    ind_alts.insert(make_pair(comb.first, x->second.second));
                    n_doub_seen.insert(make_pair(comb.first, 1));
                }
                else{
                    ind_tots[comb.first] += x->second.first + x->second.second;
                    ind_alts[comb.first] += x->second.second;
                    n_doub_seen[comb.first]++;
                }
                if (ind_tots.count(comb.second) == 0){
                    ind_tots.insert(make_pair(comb.second, x->second.first + x->second.second));
                    ind_alts.insert(make_pair(comb.second, x->second.second));
                    n_doub_seen.insert(make_pair(comb.second, 1));
                }
                else{
                    ind_tots[comb.second] += x->second.first + x->second.second;
                    ind_alts[comb.second] += x->second.second;
                    n_doub_seen[comb.second]++;
                }
            }
            else{
                if (ind_tots.count(x->first) == 0){
                    ind_tots.insert(make_pair(x->first, x->second.first + x->second.second));
                    ind_alts.insert(make_pair(x->first, x->second.second));
                    n_doub_seen.insert(make_pair(x->first, 2));
                }
                else{
                    ind_tots[x->first] += x->second.first + x->second.second;
                    ind_alts[x->first] += x->second.second;
                    n_doub_seen[x->first] += 2;
                }
            }    
        }
    }
    
    for (int i = 0; i < n_samples; ++i){
        //if (v.haps_covered[i]){
        if (v.haps_covered[i] || (n_doub_seen[i] >= 2 && ind_tots[i] >= 1.0)){
            int nalt = 0;
            if (v.haps_covered[i]){
                if (v.haps1[i]){
                    nalt++;
                }
                if (v.haps2[i]){
                    nalt++;
                }
            }
            else{
                double f = round((ind_alts[i]/ind_tots[i])/0.5)*0.5;
                if (f == 0.5){
                    nalt = 1;
                }
                else if (f == 1){
                    nalt = 2;
                }
                else{
                    nalt = 0;
                }
            }
            double frac = (double)nalt/2.0;
            if (frac == 0){
                frac = 0.001;
            }
            else if (frac == 1){
                frac = 0.999;
            }
            samp2param.insert(make_pair(i, params.size()));
            params.push_back(frac);
            samp_inds.push_back(i);
            if (nalt == 0){
                gts.push_back("0/0");
                newalt.insert(make_pair(i, 0.001));
                gts_pass.insert(0);
            }
            else if (nalt == 1){
                gts.push_back("0/1");
                newalt.insert(make_pair(i, 0.5));
                gts_pass.insert(1);
            }
            else{
                gts.push_back("1/1");
                newalt.insert(make_pair(i, 0.999));
                gts_pass.insert(2);
            }
        }
        else{
            gts.push_back("./.");
        }
    }

    vector<double> n;
    vector<double> k;
    vector<double> w;
    vector<int> idx1;
    vector<int> idx2;
    
    double tot = 0.0; 
    double alt_tot = 0.0;
    
    for (map<int, pair<float, float> >::iterator x = id_dat.begin(); x != id_dat.end(); ++x){
        if (x->second.first + x->second.second > 0){
            tot += x->second.first + x->second.second;
            alt_tot += x->second.second;
            if (x->first >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(x->first, n_samples);
                if (v.haps_covered[comb.first] && v.haps_covered[comb.second]){
                    n.push_back(x->second.first + x->second.second);
                    k.push_back(x->second.second);
                    idx1.push_back(samp2param[comb.first]);
                    idx2.push_back(samp2param[comb.second]);
                    w.push_back(weights[x->first]);
                }
            }
            else if (v.haps_covered[x->first]){
                n.push_back(x->second.first + x->second.second);
                k.push_back(x->second.second);
                idx1.push_back(samp2param[x->first]);
                idx2.push_back(-1);
                w.push_back(weights[x->first]);
            }
        }
    }
    
    bool solved = false;
    bool changed = false;

    if (n.size() > 0){    

        optimML::multivar_ml_solver solver(params, ll_regt, dll_regt);
        solver.add_data("n", n);
        solver.add_data("k", k);
        solver.add_data("idx1", idx1);
        solver.add_data("idx2", idx2);
        //solver.add_weights(w);
        for (int i = 0; i < params.size(); ++i){
            solver.constrain_01(i);
        }
        

        try{
            bool success = solver.solve();
            
            if (success){
                
                //double ll = solver.log_likelihood;
                //vector<int> gtresults;
                /*
                for (int i = 0; i < solver.results.size(); ++i){
                    fprintf(stderr, "RESULT %d) %f\n", i, solver.results[i]);
                    
                }
                */
                /*
                for (int i = 0; i < solver.results.size(); ++i){
                    double orig = solver.results[i];
                    if (solver.results[i] < 0.5){
                        solver.set_param(i, 0.001);
                        double ll1 = solver.eval_ll_all();
                        solver.set_param(i, 0.5);
                        double ll2 = solver.eval_ll_all();
                        if (ll1 > ll2){
                            gtresults.push_back(0);
                        }
                        else if (ll1 == ll2){
                            gtresults.push_back(-1);
                        }
                        else{
                            gtresults.push_back(1);
                        }
                    }
                    else if (solver.results[i] > 0.5){
                        solver.set_param(i, 0.5);
                        double ll1 = solver.eval_ll_all();
                        solver.set_param(i, 0.999);
                        double ll2 = solver.eval_ll_all();
                        if (ll1 > ll2){
                            gtresults.push_back(1);
                        }
                        else if (ll1 == ll2){
                            gtresults.push_back(-1);
                        }
                        else{
                            gtresults.push_back(2);
                        }
                    }
                    solver.set_param(i, orig);
                }
                */
                
                solved = true;
                newalt.clear();
                gts_pass.clear();
                
                for (int i = 0; i < samp_inds.size(); ++i){
                    int s = samp_inds[i];
                    //fprintf(stderr, "res %d) %f\n", s, solver.results[i]);                    
                    double rounded = round(solver.results[i]/0.5)*0.5;
                    newalt.insert(make_pair(s, solver.results[i]));
                    if (rounded == 0){
                        if (gts[s] != "0/0"){
                            changed = true;
                        }
                        gts[s] = "0/0";
                        gts_pass.insert(0);
                    }
                    else if (rounded == 0.5){
                        if (gts[s] != "0/1"){
                            changed = true;
                        }
                        gts[s] = "0/1";
                        gts_pass.insert(1);
                    }
                    else if (rounded == 1.0){
                        if (gts[s] != "1/1"){
                            changed = true;
                        }
                        gts[s] = "1/1";
                        gts_pass.insert(2);
                    }
                    /*  
                    if (gtresults[i] == -1){
                        gts[s] = "./.";
                    }
                    else if (gtresults[i] == 0){
                        gts[s] = "0/0";
                        gts_pass.insert(0);
                        newalt.insert(make_pair(s, 0.001));
                    }
                    else if (gtresults[i] == 1){
                        gts[s] = "0/1";
                        gts_pass.insert(1);
                        newalt.insert(make_pair(s, 0.5));
                    }
                    else if (gtresults[i] == 2){
                        gts[s] = "1/1";
                        gts_pass.insert(2);
                        newalt.insert(make_pair(s, 0.999));
                    }
                    */
                    /*
                    double alt = solver.results[i];
                    double d0 = dbinom(2, 0, alt);
                    double d1 = dbinom(2, 1, alt);
                    double d2 = dbinom(2, 2, alt);
                    fprintf(stderr, "idx %d alt %f d %f %f %f\n", i, alt, d0, d1, d2);

                    if (d0 > d1 && d0 > d2){
                        gts[s] = "0/0";
                        newalt.insert(make_pair(s, 0.001));
                        gts_pass.insert(0);
                    }
                    else if (d1 > d0 && d1 > d2){
                        gts[s] = "0/1";
                        newalt.insert(make_pair(s, 0.5));
                        gts_pass.insert(1);
                    }
                    else if (d2 > d0 && d2 > d1){
                        gts[s] = "1/1";
                        newalt.insert(make_pair(s, 0.999));
                        gts_pass.insert(2);
                    }
                    else{
                        gts[s] = "./.";
                    }
                    */
                }
                /*
                if (chrom == "chr1" && pos +1 == 841742){

                    fprintf(stderr, "%f %f\n", id_dat[0].first, id_dat[0].second);
                    fprintf(stderr, "AF %f\n", solver.results[samp2param[0]]);
                    for (map<int, pair<float, float> >::iterator x = id_dat.begin(); x != id_dat.end(); ++x){
                        if (x->first < n_samples){
                            fprintf(stderr, "  %d) %f %f\n", x->first, x->second.first, x->second.second);
                            fprintf(stderr, "  GT: %s\n", gts[x->first].c_str());
                        } 
                        else{
                            pair<int, int> y = idx_to_hap_comb(x->first, n_samples);
                            fprintf(stderr, "  %d+%d) %f %f\n", y.first, y.second, x->second.first, x->second.second);
                            fprintf(stderr, "  GT: %s %s\n", gts[y.first].c_str(), gts[y.second].c_str());
                        }
                    }
                    fprintf(stderr, "solved? %d\n", solved);
                    fprintf(stderr, "gts_pass %ld\n", gts_pass.size());
                    exit(0);
                }
                */
            }
        }
        catch(int x){
            // Math error. Don't re-genotype.

        }
    }
    
    /*
    double ll_all_ref = 0.0;
    double ll_not_all_ref = 0.0;

    for (map<int, pair<float, float> >::iterator x = id_dat.begin(); x != id_dat.end(); ++x){
        if (x->first >= n_samples){
            pair<int, int> comb = idx_to_hap_comb(x->first, n_samples);
            if (newalt.count(comb.first) > 0 && newalt.count(comb.second) > 0){
                ll_all_ref += dbinom(x->second.first + x->second.second, x->second.second, 0.001);
                double alt2 = 0.5*newalt[comb.first] + 0.5*newalt[comb.second];
                ll_not_all_ref += dbinom(x->second.first + x->second.second, x->second.second, alt2);                
            }
        }
        else{
            if (newalt.count(x->first) > 0){
                ll_all_ref += dbinom(x->second.first + x->second.second, x->second.second, 0.001);
                ll_not_all_ref += dbinom(x->second.first + x->second.second, x->second.second, newalt[x->first]);
            }
        }
   }
    */

   // Check whether it's a variable SNP unlikely to just be errors.
   if (gts_pass.size() > 1){
        removed = false;
        /* 
        double qual;
        double prob = pbinom(tot, alt_tot, 0.001);
        if (prob == 1){
            qual = 1e6;
        }
        else{
            qual = -10*log10(1.0-prob);
        }
        */
        /*
         double prob = pow(2,ll_not_all_ref)/(pow(2,ll_not_all_ref) + pow(2,ll_all_ref));
        double qual;
        if (prob == 1){
            qual = 10000;
        }
        else{
            qual = -10*log10(1.0-prob);
        }
        */

        fprintf(stdout, "%s\t%d\t.\t%c\t%c\t%.3f\t.\t.\tGT", chrom.c_str(), pos+1, 
            v.ref, v.alt, v.vq);
        for (int i = 0; i < gts.size(); ++i){
            fprintf(stdout, "\t%s", gts[i].c_str());
        }    
        fprintf(stdout, "\n");
    }
    else{
        removed = true;
        solved = false;
        changed = false;
    }
    return changed;
}

void compute_weights(robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignment_llr,
    map<int, double>& weights){
    
    double weightsum = 0.0;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assignments.begin();
        a != assignments.end(); ++a){
        weightsum += assignment_llr[a->first];
        if (weights.count(a->second) == 0){
            weights.insert(make_pair(a->second, 0.0));
        }
        weights[a->second] += assignment_llr[a->first];
    }
    for (map<int, double>::iterator w = weights.begin(); w != weights.end(); ++w){
        w->second /= weightsum;
    }
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"assignments", required_argument, 0, 'a'},
       {"index_jump", no_argument, 0, 'j'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile = "";
    string vcf_file = "";
    string assnfile = "";
    bool stream = true;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:a:jh", long_options, &option_index )) != -1){
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
            case 'a':
                assnfile = optarg;
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
    if (bamfile == ""){
        fprintf(stderr, "ERROR: BAM file required\n");
        exit(1);
    }    
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    if (assnfile == ""){
        fprintf(stderr, "ERROR: assignments file required\n");
        exit(1);
    }

    
    // Init BAM reader
    bam_reader reader = bam_reader();
    reader.set_file(bamfile);

    // Store the names of all individuals in the VCF
    vector<string> samples;
    
    // Store data about SNPs
    map<int, map<int, var> > snpdat;
    
    int nsnps;
    
    fprintf(stderr, "Loading variants from VCF/BCF...\n");
    nsnps = read_vcf(vcf_file, reader, samples, snpdat, 0.0, false, false);
    
    // Write VCF header to stdout
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);

    // Open output VCF (stdout)
    htsFile* outf = hts_open("-", "w");
    int write_success = bcf_hdr_write(outf, bcf_header);

    fprintf(stderr, "Loading assignments...\n");
    robin_hood::unordered_map<unsigned long, int> assignments;
    robin_hood::unordered_map<unsigned long, double> assignment_llr;
    parse_assignments(assnfile, assignments, assignment_llr, samples);
    
    // Compute weights
    map<int, double> weights;
    compute_weights(assignments, assignment_llr, weights);

    // Data structure to store SNP data
    map<int, map<int, map<int, pair<float, float> > > > snp_id_counts;

    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;

    fprintf(stderr, "Counting alleles in BAM file...\n");

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
    
    int n_updated = 0;
    int n_rm = 0;
    int n_tot = 0;

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
                        if (snpdat[curtid].count(cursnp->first) > 0){        
                            bool rm = false; 
                            bool updated = regt_snp(tid2chrom[curtid], cursnp->first,
                                snpdat[curtid][cursnp->first],
                                samples.size(),
                                snp_id_counts[curtid][cursnp->first],
                                rm,
                                weights);
                            if (updated){
                                n_updated++;
                            }
                            else if (rm){
                                n_rm++;
                            }
                            n_tot++;
                        }
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
                if (snp_id_counts[reader.tid()].count(cursnp->first) > 0){ 
                    bool rm = false;
                    bool updated = regt_snp(tid2chrom[reader.tid()], cursnp->first,
                        cursnp->second,
                        samples.size(),
                        snp_id_counts[reader.tid()][cursnp->first],
                        rm,
                        weights);
                    if (updated){
                        n_updated++;
                    }
                    else if (rm){
                        n_rm++;
                    }
                    n_tot++;
                }
                ++nsnp_processed;
                ++cursnp;
            }
            // Create a second iterator to look ahead for any additional SNPs 
            // within the current read
            map<int, var>::iterator cursnp2 = cursnp;
            while (cursnp2 != snpdat[reader.tid()].end() && 
                cursnp2->first >= reader.reference_start && 
                cursnp2->first <= reader.reference_end){   
                process_bam_record_bysnp(reader, cursnp2->first, cursnp2->second,
                    assignments, snp_id_counts[reader.tid()][cursnp2->first]);
                ++cursnp2;
            }
            if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, nsnps); 
                last_print = nsnp_processed;
            }
        }
        // Handle any final SNPs.
        if (curtid != -1){
            while (cursnp != snpdat[curtid].end()){
                if (snpdat[curtid].count(cursnp->first) > 0){        
                    bool rm = false; 
                    bool updated = regt_snp(tid2chrom[curtid], cursnp->first,
                        snpdat[curtid][cursnp->first],
                        samples.size(),
                        snp_id_counts[curtid][cursnp->first],
                        rm,
                        weights);
                    if (updated){
                        n_updated++;
                    }
                    else if (rm){
                        n_rm++;
                    }
                    n_tot++;
                }
                ++nsnp_processed;
                ++cursnp;    
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
                    process_bam_record_bysnp(reader, cursnp->first, cursnp->second,
                        assignments, snp_id_counts[reader.tid()][cursnp->first]);
                }
                
                bool rm;
                bool updated = regt_snp(tid2chrom[tid], cursnp->first,
                    cursnp->second,
                    samples.size(),
                    snp_id_counts[tid][cursnp->first],
                    rm,
                    weights);
                if (updated){
                    n_updated++;
                }
                else if (rm){
                    n_rm++;
                }
                n_tot++;
                ++nsnp_processed;        
            
                if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                    fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, 
                        nsnps); 
                    last_print = nsnp_processed;
                }
            }
        }
    }
    hts_close(outf);

    fprintf(stderr, "Processed %d of %d SNPs\n", nsnp_processed, nsnps);    
   
    fprintf(stderr, "Updated %d of %d (%.2f%%)\n", n_updated, n_tot, 100*(double)n_updated/(double)n_tot);
    fprintf(stderr, "Removed %d of %d (%.2f%%)\n", n_rm, n_tot, 100*(double)n_rm/(double)n_tot);
    return 0;
}
