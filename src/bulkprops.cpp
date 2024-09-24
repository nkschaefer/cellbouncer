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
#include <random>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bam.h>
#include <optimML/mixcomp.h>
#include <optimML/multivar_ml.h>
#include <mixtureDist/functions.h>
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
    fprintf(stderr, "    --bootstrap -B If you wish to compute the variance of the estimate of\n");
    fprintf(stderr, "       pool composition, set a number of bootstrap replicates to perform here.\n");
    fprintf(stderr, "       After all replicates are computed, a Dirichlet distribution will be fit\n");
    fprintf(stderr, "       and MLE Dirichlet parameters will be written to the output file as a third\n");
    fprintf(stderr, "       column. Default = 100\n");
    fprintf(stderr, "    --qual -q Minimum variant quality to consider (default 50)\n");
    fprintf(stderr, "    --ids -i If the VCF file contains individuals that you do not\n");
    fprintf(stderr, "       expect to see in your sample, specify individuals to include here.\n");
    fprintf(stderr, "       These identities will be considered, as well as all possible doublet\n");
    fprintf(stderr, "       combinations of them. Expects a text file, with one individual ID per\n");
    fprintf(stderr, "       line, matching individual names in the VCF.\n");
    fprintf(stderr, "    --n_trials -N Initial guesses will be randomly shuffled a number of times equal\n");
    fprintf(stderr, "       to this number, in order to avoid reporting a local maximum, Default = 0\n");
    fprintf(stderr, "    --error_rate -e Sequencing error rate (will attempt to calculate from data if \n");
    fprintf(stderr, "       not provided.\n");
    fprintf(stderr, "----- I/O options -----\n");
    fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    fprintf(stderr, "===== ALTERNATIVE RUN MODE =====\n");
    fprintf(stderr, "    --props -p A preexisting file listing mixture proportions (or an .assignments\n");
    fprintf(stderr, "       file from which to calculate them). In this mode, instead of inferring MLE\n");
    fprintf(stderr, "       mixture proportions, it will compute the log likelihood of the read counts\n");
    fprintf(stderr, "       at every SNP given the mixture proportions supplied here, and output in BED\n");
    fprintf(stderr, "       format to stdout.\n");
    fprintf(stderr, "    --genes -g Instead of printing log likelihood for every SNP, aggregate across\n");
    fprintf(stderr, "       genes and output one average LL (weighted by number of reads at each SNP)\n");
    fprintf(stderr, "       per gene. Requires CellRanger/STARsolo-like GN/GX tags to be present in BAM.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

/**
 * Finds MLE of mixture proportions, given counts and expectations
 */
void infer_props_aux(vector<double>& n,
    vector<double>& k,
    vector<vector<double> >& expfracs_all,
    int n_samples,
    int n_trials,
    double& maxll,
    vector<double>& maxprops){

    vector<double> startfrac;
    for (int i = 0; i < n_samples; ++i){
        startfrac.push_back(1.0/(double)n_samples);
    }    
    optimML::mixcomp_solver solver(expfracs_all, "binom", n, k);
    solver.solve();
    
    maxll = solver.log_likelihood;
    maxprops = solver.results;
    for (int i = 0; i < n_trials; ++i){
        solver.add_mixcomp_fracs(startfrac);
        solver.randomize_mixcomps();
        solver.solve();
        if (solver.log_likelihood > maxll){
            maxll = solver.log_likelihood;
            maxprops = solver.results;
        }
    }
}



bool infer_props(map<int, map<int, pair<float, float> > >& snp_ref_alt,
    map<int, map<int, var> >& snpdat,
    double err_rate,
    int n_trials,
    int n_samples,
    double& nreads,
    double& maxll,
    vector<double>& maxprops,
    vector<double>& dirichlet_mle,
    int bootstrap){
    
    nreads = 0;
    
    // Transform data into counts & expectations
    vector<double> n;
    vector<double> k;
    vector<vector<double> > expfracs_all;

    for (map<int, map<int, pair<float, float> > >::iterator s = snp_ref_alt.begin(); s != snp_ref_alt.end();
        ++s){
        for (map<int, pair<float, float> >::iterator s2 = s->second.begin(); s2 != s->second.end(); ++s2){
            
            bool miss = false;
            vector<double> expfracs;

            for (int i = 0; i < n_samples; ++i){
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
                nreads += s2->second.first + s2->second.second;
                n.push_back(s2->second.first + s2->second.second);
                k.push_back(s2->second.second);
                expfracs_all.push_back(expfracs);          
            } 
        }
    }
    
    if (expfracs_all.size() == 0){
        return false;
    }
    
    infer_props_aux(n, k, expfracs_all, n_samples, n_trials, maxll, maxprops);
    
    if (bootstrap > 0){
        // Resample.
        // Init random number generator - use static variables so it only 
        // happens once
        static bool init = false;
        static random_device dev;
        static mt19937 rand_gen;
        static uniform_int_distribution<int> uni_dist;
        if (!init){
            rand_gen = mt19937(dev());
            uni_dist = uniform_int_distribution<int>(0, n.size()-1);
            init = true;
        }
        
        vector<vector<double> > props_bootstrap;
        vector<vector<double> > dirprops;
        for (int i = 0; i < n_samples; ++i){
            vector<double> v;
            dirprops.push_back(v);
        }

        for (int b = 0; b < bootstrap; ++b){
            fprintf(stderr, "Bootstrap sample %d...\r", b+1);

            // Sample with replacement.
            vector<double> n_bootstrap;
            vector<double> k_bootstrap;
            vector<vector<double> > expfrac_bootstrap;
            for (int x = 0; x < n.size(); ++x){
                int r = uni_dist(rand_gen); 
                n_bootstrap.push_back(n[r]);
                k_bootstrap.push_back(k[r]);
                expfrac_bootstrap.push_back(expfracs_all[r]);
            }
            vector<double> p_bootstrap;
            double ll_bootstrap;
            // Solve
            infer_props_aux(n_bootstrap, k_bootstrap, expfrac_bootstrap,
                n_samples, n_trials, ll_bootstrap, p_bootstrap);

            props_bootstrap.push_back(p_bootstrap);
            for (int j = 0; j < n_samples; ++j){
                dirprops[j].push_back(p_bootstrap[j]);
            }
        }
        fprintf(stderr, "\n");
        
        // Find MLE concentration parameters of Dirichlet (in common.cpp)
        fit_dirichlet(maxprops, dirprops, dirichlet_mle);
        
    }
    return true;
}

/**
 * Splits a line (from an input text file) into fields using tabs as separators.
 */
void split_line(string& line, vector<string>& fields){
    istringstream splitter(line);
    string field;
    while(getline(splitter, field, '\t')){
        fields.push_back(field);
    }    
}

/**
 * Parse an input file (assignments or bulkprops) and convert to proportions of 
 * individuals. If assignments, treat as bulk (split up doublet identities).
 */
void parse_props_prev(string& filename, map<int, double>& props, vector<string>& samples){
    map<string, int> name2idx;
    for (int i = 0; i < samples.size(); ++i){
        name2idx.insert(make_pair(samples[i], i));
    }
    ifstream inf(filename);
    string line;
    bool is_assignments;
    bool first = true;

    map<int, int> assncounts;
    int assntot;

    while (getline(inf, line)){
        vector<string> fields;
        split_line(line, fields);
        if (first){
            // Check whether assignments file.
            if (fields.size() == 4 && (fields[2] == "S" || fields[2] == "D")){
                is_assignments = true;
            }
            else{
                is_assignments = false;
            }
            first = false;
        }
        if (is_assignments){
            if (fields[2] == "D"){
                size_t splitpos = fields[1].find("+");
                if (splitpos == string::npos){
                    fprintf(stderr, "ERROR parsing doublet identity %s\n", fields[1].c_str());
                    exit(1);
                }
                string name1 = fields[1].substr(0, splitpos);
                string name2 = fields[1].substr(splitpos+1, fields[1].length()-splitpos-1);
                if (name2idx.count(name1) == 0){
                    fprintf(stderr, "ERROR: name %s not found in VCF data.\n", name1.c_str());
                }
                if (name2idx.count(name2) == 0){
                    fprintf(stderr, "ERROR: name %s not found in VCF data.\n", name2.c_str());
                }
                int id1 = name2idx[name1];
                int id2 = name2idx[name2];
                if (assncounts.count(id1) == 0){
                    assncounts.insert(make_pair(id1, 0));
                }
                if (assncounts.count(id2) == 0){
                    assncounts.insert(make_pair(id2, 0));
                }
                assncounts[id1]++;
                assncounts[id2]++;
                assntot += 2; 
            }
            else{
                if (name2idx.count(fields[1]) == 0){
                    fprintf(stderr, "ERROR: name %s not found in VCF data.\n", fields[1].c_str());
                    exit(1);
                }
                int id = name2idx[fields[1]];
                if (assncounts.count(id) == 0){
                    assncounts.insert(make_pair(id, 0));
                }
                assncounts[id] += 2;
                assntot += 2;
            }
        }
        else{
            double prop = atof(fields[1].c_str());
            if (name2idx.count(fields[0]) == 0){
                fprintf(stderr, "ERROR: name %s not found in VCF data.\n", fields[0].c_str());
                exit(1);
            }
            int id = name2idx[fields[0]];
            props.insert(make_pair(id, prop));
        }
    }
    if (is_assignments){
        // Convert counts into proportions.
        for (map<int, int>::iterator a = assncounts.begin(); a != assncounts.end(); ++a){
            double prop = (double)a->second/(double)assntot;
            props.insert(make_pair(a->first, prop));
        }
    }
}

void compute_ll_snps(map<int, map<int, pair<float, float> > >& snp_ref_alt,
    map<int, map<int, var> >& snpdat,
    map<int, string>& tid2chrom,
    map<int, double>& props,
    double err_rate,
    int n_samples,
    bool genes,
    map<pair<int, int>, set<string> >& snp_gene_ids,
    map<string, string>& gene_id2name){
    
    map<string, double> genesums;
    map<string, double> genecounts;
     
    for (map<int, map<int, pair<float, float> > >::iterator x = snp_ref_alt.begin();
        x != snp_ref_alt.end(); ++x){
        for (map<int, pair<float, float> >::iterator y = x->second.begin();
            y != x->second.end(); ++y){
            bool miss = false;
            // Get expected freqs
            double freqsum = 0.0;
            for (int i = 0; i < n_samples; ++i){
                if (!snpdat[x->first][y->first].haps_covered[i]){
                    miss = true;
                    break;
                }
                else{
                    int nalt = 0;
                    if (snpdat[x->first][y->first].haps1[i]){
                        nalt++;
                    }
                    if (snpdat[x->first][y->first].haps2[i]){
                        nalt++;
                    }
                    double expec = (double)nalt/(double)2.0;
                    if (expec == 0){
                        expec += err_rate;
                    }
                    else if (expec == 1.0){
                        expec -= err_rate;
                    }
                    freqsum += props[i] * expec;
                }
            }
            if (!miss){
                double ref = y->second.first;
                double alt = y->second.second;
                if (ref + alt > 0){
                    double ll = dbinom(ref+alt, alt, freqsum);
                    if (genes){
                        pair<int, int> key = make_pair(x->first, y->first);
                        for (set<string>::iterator gid = snp_gene_ids[key].begin(); gid != 
                            snp_gene_ids[key].end(); ++gid){
                            if (genesums.count(*gid) == 0){
                                genesums.insert(make_pair(*gid, 0));
                                genecounts.insert(make_pair(*gid, 0));
                            }
                            genesums[*gid] += ll*(ref+alt);
                            genecounts[*gid] += ref+alt;
                        }
                    }
                    else{
                        fprintf(stdout, "%s\t%d\t%d\t%f\t%f\n", tid2chrom[x->first].c_str(),
                            y->first, y->first+1, ll, ref+alt);
                    }
                }
            }
        }
    }
    if (genes){
        for (map<string, double>::iterator gs = genesums.begin(); gs != genesums.end(); ++gs){
            double mean = gs->second/genecounts[gs->first];
            string name = gs->first;
            if (gene_id2name.count(gs->first) > 0){
                name = gene_id2name[gs->first];
            }
            fprintf(stdout, "%s\t%s\t%f\t%f\n", gs->first.c_str(), name.c_str(), mean, genecounts[gs->first]);
        }
    }
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"index_jump", no_argument, 0, 'j'},
       {"ids", required_argument, 0, 'i'},
       {"qual", required_argument, 0, 'q'},
       {"n_trials", required_argument, 0, 'N'},
       {"props", required_argument, 0, 'p'},
       {"bootstrap", required_argument, 0, 'B'},
       {"genes", no_argument, 0, 'g'},
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
    int n_trials = 0;
    double err_prior = -1;
    int bootstrap = 100;
    bool bootstrap_given = false;
    string propsfile;
    bool props_given = false;
    bool genes = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:e:v:o:q:i:N:w:B:p:gjh", long_options, &option_index )) != -1){
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
            case 'e':
                err_prior = atof(optarg);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'N':
                n_trials = atoi(optarg);
                break;
            //case 'w':
            //    window = atoi(optarg);
            //    break;
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
            case 'B':
                bootstrap = atoi(optarg);
                bootstrap_given = true;
                break;
            case 'p':
                propsfile = optarg;
                props_given = true;
                break;
            case 'g':
                genes = true;
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
    if (output_prefix.length() == 0 && !props_given){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        exit(1);
    }
    if (n_trials < 0){
        fprintf(stderr, "ERROR: --n_trials/-N must be 0 or greater\n");
        exit(1);
    }
    if (props_given && idfile_given){
        fprintf(stderr, "ERROR: only one of --props/-p and --ids/-i is allowed.\n");
        exit(1);
    }
    if (genes && !props_given){
        fprintf(stderr, "ERROR: --genes/-g argument requires --props/-p.\n");
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
    
    map<int, double> props_prev;
    if (props_given){
        parse_props_prev(propsfile, props_prev, samples);

    }

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

    // Store gene names at SNPs
    map<pair<int, int>, set<string> > snp_gene_ids;
    map<string, string> snp_gene_names;

    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;

    // initialize the BAM reader
    reader.set_file(bamfile);
    
    // retrieve cell barcodes
    reader.set_cb();
    
    if (genes){
        reader.check_genes();
    }

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
                
                if (!genes || (reader.gene_ids.size() > 0 || reader.gene_names.size() > 0)){
                    process_bam_record_bulk(reader, cursnp2->first, cursnp2->second,
                        snp_ref_alt, snp_err, genes, snp_gene_ids, snp_gene_names);
                }

                ++cursnp2;
            }
            if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, nsnps); 
                last_print = nsnp_processed;
            }
        
            if (nsnp_processed == nsnps){
                break;
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
                    if (!genes || (reader.gene_ids.size() > 0 || reader.gene_names.size() > 0)){
                        process_bam_record_bulk(reader, cursnp->first, cursnp->second, 
                            snp_ref_alt, snp_err, genes, snp_gene_ids, snp_gene_names);
                    }
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
    
    if (snp_ref_alt.size() == 0){
        fprintf(stderr, "ERROR: no valid SNPs detected.\n");
        if (genes){
            fprintf(stderr, "Did you run this with --genes on a BAM file that lacks GX/GN tags?\n");
            fprintf(stderr, "make sure these tags were added by an alignment program like CellRanger\n");
            fprintf(stderr, "or STARsolo.\n");
            exit(1);
        }
    }
    // Calculate error rate
    double err_rate = 0.0;
    if (err_prior > 0 && err_prior < 1){
        err_rate = err_prior;
    }
    else{
        int nsnp_data = 0;
        for (map<int, map<int, float> >::iterator s = snp_err.begin(); s != snp_err.end(); ++s){
            for (map<int, float>::iterator s2 = s->second.begin(); s2 != s->second.end(); ++s2){
                nsnp_data++;
            }
        }
        double frac = 1.0/(double)nsnp_data;
        err_rate = 0.0;
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
    }
    if (err_rate == 0.0 || isnan(err_rate)){
        fprintf(stderr, "TRUE\n");
        err_rate = 0.001;
    }
    fprintf(stderr, "Approximate error rate = %f\n", err_rate);
    
    string outfn = output_prefix + ".bulkprops";
    FILE* outf = fopen(outfn.c_str(), "w");
    
    // Window output not used.
    if (false){
        /*
    //if (window > 0){
        fprintf(stdout, "chrom\tstart\tend\tnreads");
        map<string, int> samplesort;
        for (int i = 0; i < samples.size(); ++i){
            samplesort.insert(make_pair(samples[i], i));
        }
        for (map<string, int>::iterator s = samplesort.begin(); s != samplesort.end(); ++s){
            fprintf(stdout, "\t%s", s->first.c_str());
        }
        fprintf(stdout, "\n");

        int winstart = -1;
        map<int, map<int, pair<float, float> > > snpcounts_win;
        for (map<int, map<int, pair<float, float> > >::iterator x = snp_ref_alt.begin();
            x != snp_ref_alt.end(); ++x){
            map<int, pair<float, float> > m;
            snpcounts_win.insert(make_pair(x->first, m));
            for (map<int, pair<float, float> >::iterator y = x->second.begin(); y != 
                x->second.end(); ++y){
                if (y->second.first + y->second.second > 0){
                    snpcounts_win[x->first].insert(make_pair(y->first, y->second));
                    if (winstart == -1){
                        winstart = y->first;
                    }
                    int winend = y->first;
                    if (snpcounts_win[x->first].size() == window){
                        double ll;
                        vector<double> props;
                        double nreads;
                        vector<double> dummy;
                        bool success = infer_props(snpcounts_win, snpdat, 
                            err_rate, n_trials, samples.size(), nreads, ll, props, dummy, bootstrap);
                        if (success){
                            fprintf(stdout, "%s\t%d\t%d\t%f", tid2chrom[x->first].c_str(),
                                winstart, winend+1, nreads);
                            for (map<string, int>::iterator s = samplesort.begin(); s != 
                                samplesort.end(); ++s){
                                fprintf(stdout, "\t%f", props[s->second]);
                            }
                            fprintf(stdout, "\n");
                        }
                        snpcounts_win[x->first].erase(snpcounts_win[x->first].begin()->first);
                        winstart = -1;
                        if (snpcounts_win[x->first].size() > 0){
                            winstart = snpcounts_win[x->first].begin()->first;
                        }

                        //snpcounts_win[x->first].clear();
                        //winstart = -1;
                    }
                }
            }
            // Finish up with this chromosome.
            if (snpcounts_win[x->first].size() > 0){
                int winend = snpcounts_win[x->first].rbegin()->first;
                double ll;
                vector<double> props;
                double nreads;
                vector<double> dummy;
                bool success = infer_props(snpcounts_win, snpdat, err_rate, 
                    n_trials, samples.size(), nreads, ll, props, dummy, bootstrap);
                snpcounts_win.erase(x->first);
                if (success){
                    fprintf(stdout, "%s\t%d\t%d\t%f", tid2chrom[x->first].c_str(),
                        winstart, winend+1, nreads);
                    for (map<string, int>::iterator s = samplesort.begin(); s != 
                        samplesort.end(); ++s){
                        fprintf(stdout, "\t%f", props[s->second]);
                    }
                    fprintf(stdout, "\n");
                }
                winstart = -1;
            }
        }
        */
    }
    else if (props_given){
        // Compute log likelihood of each SNP under given proportions.
        compute_ll_snps(snp_ref_alt, snpdat, tid2chrom, 
            props_prev, err_rate, samples.size(), genes, snp_gene_ids, snp_gene_names);
    }
    else{
        vector<double> props;
        double ll;
        double nreads;
        vector<double> dirichlet_params;
        infer_props(snp_ref_alt, snpdat, err_rate, n_trials, samples.size(), 
            nreads, ll, props, dirichlet_params, bootstrap);
        
        int si = 0;
        map<string, int> samplesort;
        for (int i = 0; i < samples.size(); ++i){
            samplesort.insert(make_pair(samples[i], i));
        }

        for (map<string, int>::iterator s = samplesort.begin(); s != samplesort.end(); ++s){
            if (bootstrap > 0){
                fprintf(outf, "%s\t%f\t%f\n", s->first.c_str(), props[s->second], 
                    dirichlet_params[s->second]);
            }
            else{
                fprintf(outf, "%s\t%f\n", s->first.c_str(), props[s->second]);
            }
        }
    }
    fclose(outf);
    
    return 0;
}
