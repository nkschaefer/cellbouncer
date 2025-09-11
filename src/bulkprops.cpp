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
    fprintf(stderr, "    --num_threads -T Number of parallel threads to use in maximum likelihood\n");
    fprintf(stderr, "       inference (does not apply to BAM reading) (default = none)\n");
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
    //fprintf(stderr, "----- I/O options -----\n");
    //fprintf(stderr, "    --index_jump -j Instead of reading through the entire BAM file \n");
    //fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    //fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    //fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    //fprintf(stderr, "       of SNPs.\n");
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
    vector<double>& maxprops,
    int nthread){

    vector<double> startfrac;
    for (int i = 0; i < n_samples; ++i){
        startfrac.push_back(1.0/(double)n_samples);
    }    
    optimML::mixcomp_solver solver(expfracs_all, "binom", n, k);
    if (nthread > 0){
        // Because of how "mixprops" are handled internally, it's slow to 
        // multithread them (stuff needs to be recalculated for every data
        // pt for every thread). Just multi-thread the solving part, 
        // not the evaluation part.
        solver.set_threads_bfgs(nthread);
    }
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

bool infer_props(vector<double>& n,
    vector<double>& k,
    vector<vector<double> >& expfracs_all,
    int n_trials,
    int n_samples,
    double& maxll,
    vector<double>& maxprops,
    vector<double>& dirichlet_mle,
    int bootstrap,
    int nthread){
    
    if (expfracs_all.size() == 0){
        return false;
    }
    
    infer_props_aux(n, k, expfracs_all, n_samples, n_trials, maxll, maxprops, nthread);
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
                //vector<double> row(n_samples);
                //for (int j = 0; j < n_samples; ++j){
                //    row[j] = expfracs_all[r][j];
                //}
                //expfrac_bootstrap.push_back(row);
                expfrac_bootstrap.push_back(expfracs_all[r]);
            }
            vector<double> p_bootstrap;
            double ll_bootstrap;
            
            // Solve
            infer_props_aux(n_bootstrap, k_bootstrap, expfrac_bootstrap,
                n_samples, n_trials, ll_bootstrap, p_bootstrap, nthread);
            props_bootstrap.push_back(p_bootstrap);
            for (int j = 0; j < n_samples; ++j){
                dirprops[j].push_back(p_bootstrap[j]);
            }
        }
        // Find MLE concentration parameters of Dirichlet (in common.cpp)
        fit_dirichlet(maxprops, dirprops, dirichlet_mle, nthread);
        
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

int compute_ll_snp(pair<float, float>& snp_ref_alt,
    var& snpdat,
    string& chrom,
    int tid, 
    int pos,
    map<int, double>& props,
    double err_rate,
    int n_samples,
    bool genes,
    map<pair<int, int>, set<string> >& snp_gene_ids,
    map<string, string>& gene_id2name,
    map<string, double>& genesums,
    map<string, double>& genecounts){
    
    // Get expected freqs
    double freqsum = 0.0;
    for (int i = 0; i < n_samples; ++i){
        int nalt = 0;
        if (snpdat.haps1[i]){
            nalt++;
        }
        if (snpdat.haps2[i]){
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

    double ref = snp_ref_alt.first;
    double alt = snp_ref_alt.second;
    
    if (ref + alt > 0){
        double ll = dbinom(ref+alt, alt, freqsum);
        if (genes){
            pair<int, int> key = make_pair(tid, pos);
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
            fprintf(stdout, "%s\t%d\t%d\t%f\t%f\n", chrom.c_str(),
                pos, pos+1, ll, ref+alt);
        }
        return 1;
    }
    return 0;
    
}

void summarize_data(pair<float, float>& snp_ref_alt,
    float snp_err_count,
    var& snpdat,
    int n_samples,
    double err_rate,
    double& nreads,
    vector<double>& n,
    vector<double>& k,
    vector<vector<double> >& expfracs_all,
    vector<double>& empirical_err_rates){
   
    if (snp_ref_alt.first + snp_ref_alt.second <= 1){
        return;
    } 

    vector<double> expfracs;
    
    if (err_rate <= 0 || err_rate >= 1){
        double e = snp_err_count/(snp_err_count + snp_ref_alt.first + 
            snp_ref_alt.second);
        empirical_err_rates.push_back(e);
    }

    for (int i = 0; i < n_samples; ++i){
        int nalt = 0;
        // Note: we have already ensured no missing genotypes when
        // loading SNPs
        if (snpdat.haps1[i]){
            nalt++;
        }
        if (snpdat.haps2[i]){
            nalt++;
        }
        double expfrac = (double)nalt/2.0;
        
        // Account for error rate
        // Can only do this if the user supplied an error rate -- if not,
        // will need to do this later after it's calculated.
        if (err_rate > 0 && err_rate < 1){
            if (expfrac == 0){
                expfrac += err_rate;
            }
            else if (expfrac == 1.0){
                expfrac -= err_rate;
            }
        }
        expfracs.push_back(expfrac);
    }

    nreads += snp_ref_alt.first + snp_ref_alt.second;
    n.push_back(snp_ref_alt.first + snp_ref_alt.second);
    k.push_back(snp_ref_alt.second);
    expfracs_all.push_back(expfracs);          
}

int main(int argc, char *argv[]) {    
   
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       //{"index_jump", no_argument, 0, 'j'},
       {"ids", required_argument, 0, 'i'},
       {"qual", required_argument, 0, 'q'},
       {"n_trials", required_argument, 0, 'N'},
       {"props", required_argument, 0, 'p'},
       {"bootstrap", required_argument, 0, 'B'},
       {"genes", no_argument, 0, 'g'},
       {"error_rate", required_argument, 0, 'e'},
       {"num_threads", required_argument, 0, 'T'},
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
    int nthread = 0;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:e:v:o:q:i:N:w:B:p:T:gh", long_options, &option_index )) != -1){
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
            //case 'j':
            //    stream = false;
            //    break;
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
            case 'T':
                nthread = atoi(optarg);
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
    else if (output_prefix.length() > 0 && is_dir(output_prefix)){
        fprintf(stderr, "ERROR: output_prefix %s is a directory, but should be a file \
name prefix.\n", output_prefix.c_str());
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
    if (props_given && (err_prior <= 0 || err_prior >= 1)){
        fprintf(stderr, "ERROR: if providing pre-computed proportions, you must \
also provide a valid --error_rate/-e.\n");
        exit(1);
    } 
    if (nthread <= 1){
        nthread = 0;
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
    map<int, var> snpdat;
    
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: vcf file is required\n");
        exit(1);
    }
    
    // Get samples from VCF
    read_vcf_samples(vcf_file, samples);

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
    map<int, pair<float, float> > snp_ref_alt;
    
    // Store erronenous (non ref-alt) counts at SNPs
    map<int, float> snp_err;

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
    
    vector<double> n;
    vector<double> k;
    vector<vector<double> > expfracs_all;
    vector<double> empirical_err_rates;
    double nreads = 0.0;
    
    map<string, double> genesums;
    map<string, double> genecounts;

    int ll_snps = 0;

    int nsnp_processed = 0;
    if (stream){

        // Read through entire BAM file and look for informative SNPs along the way
        // (default setting, appropriate for large numbers of SNPs).
        int curtid = -1;
        
        map<int, var>::iterator cursnp;
        while (reader.next()){
            
            if (reader.unmapped() || reader.secondary() || reader.supplementary() ||
                reader.qcfail() || reader.dup()){
                continue;
            }

            if (curtid != reader.tid()){
                // Started a new chromosome
                if (curtid != -1){
                    while (cursnp != snpdat.end()){
                        if (props_given){
                            ll_snps += compute_ll_snp(snp_ref_alt[cursnp->first], cursnp->second, 
                                tid2chrom[curtid], curtid, cursnp->first, props_prev, err_prior, 
                                samples.size(), genes, snp_gene_ids, snp_gene_names, genesums,
                                genecounts);
                        }
                        else{
                            summarize_data(snp_ref_alt[cursnp->first], snp_err[cursnp->first], 
                                cursnp->second, samples.size(), err_prior, nreads, n, k, 
                                expfracs_all, empirical_err_rates);
                        }
                        snp_ref_alt.erase(cursnp->first);
                        snp_err.erase(cursnp->first);
                        ++nsnp_processed;
                        snpdat.erase(cursnp++);
                    }
                }
                
                snpdat.clear();
                char* curchromptr = reader.ref_id();
                if (curchromptr != NULL){    
                    string curchrom = curchromptr;
                    // The last argument here is very important: do not load any SNPs
                    // where there are missing genotypes. They are impossible to 
                    // model.
                    read_vcf_chrom(vcf_file, curchrom, snpdat, vq, false);
                    cursnp = snpdat.begin();
                }
                else{
                    cursnp = snpdat.end();
                }
                curtid = reader.tid();
            }
            // Advance to position within cur read
            while (cursnp != snpdat.end() && 
                cursnp->first < reader.reference_start){
                if (props_given){
                    ll_snps += compute_ll_snp(snp_ref_alt[cursnp->first], cursnp->second, 
                        tid2chrom[curtid], curtid, cursnp->first, props_prev, err_prior, 
                        samples.size(), genes, snp_gene_ids, snp_gene_names, genesums,
                        genecounts);
                }
                else{
                    summarize_data(snp_ref_alt[cursnp->first], snp_err[cursnp->first],
                        cursnp->second, samples.size(), err_prior, nreads, n, k, 
                        expfracs_all, empirical_err_rates);
                }
                snp_ref_alt.erase(cursnp->first);
                snp_err.erase(cursnp->first);
                ++nsnp_processed;
                snpdat.erase(cursnp++);
            }
            
            // Create a second iterator to look ahead for any additional SNPs 
            // within the current read
            map<int, var>::iterator cursnp2 = cursnp;
            while (cursnp2 != snpdat.end() && 
                cursnp2->first >= reader.reference_start && 
                cursnp2->first <= reader.reference_end){   
                
                if (!genes || (reader.gene_ids.size() > 0 || reader.gene_names.size() > 0)){
                    process_bam_record_bulk(reader, cursnp2->first, cursnp2->second,
                        snp_ref_alt, snp_err, genes, snp_gene_ids, snp_gene_names);
                }

                ++cursnp2;
            }
            if (nsnp_processed % progress == 0 && nsnp_processed > last_print){
                fprintf(stderr, "Processed %d SNPs\r", nsnp_processed); 
                last_print = nsnp_processed;
            }
        }
        
        if (snp_ref_alt.size() > 0){
            while (cursnp != snpdat.end()){
                if (props_given){
                    ll_snps += compute_ll_snp(snp_ref_alt[cursnp->first], cursnp->second, 
                        tid2chrom[curtid], curtid, cursnp->first, props_prev, err_prior, 
                        samples.size(), genes, snp_gene_ids, snp_gene_names, genesums,
                        genecounts);           
                }
                else{
                    summarize_data(snp_ref_alt[cursnp->first], snp_err[cursnp->first],
                        cursnp->second, samples.size(), err_prior, nreads, n, k, 
                        expfracs_all, empirical_err_rates);
                }
                snp_ref_alt.erase(cursnp->first);
                snp_err.erase(cursnp->first);
                snpdat.erase(cursnp++);
                ++nsnp_processed;
            }
        }
    }
    else{
        /*
         // Index-jumping incompatible with loading SNPs one chromosome at a time.
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
        */
    }
    fprintf(stderr, "Processed %d SNPs\n", nsnp_processed);
    
    if (ll_snps == 0 && n.size() == 0){
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
    bool used_prior_err = true;
    if (err_prior > 0 && err_prior < 1){
        err_rate = err_prior;
    }
    else{
        double frac = 1.0/(double)empirical_err_rates.size();
        for (int i = 0; i < empirical_err_rates.size(); ++i){
            err_rate += frac * empirical_err_rates[i];
        }
        used_prior_err = false;
    }
    
    if (err_rate == 0.0 || isnan(err_rate)){
        // Fall back on minimum error rate
        err_rate = 0.001;
        used_prior_err = false;
    }

    fprintf(stderr, "Approximate error rate = %f\n", err_rate);
    
    if (!used_prior_err){
        // Adjust
        for (int i = 0; i < expfracs_all.size(); ++i){
            for (int j = 0; j < expfracs_all[i].size(); ++j){
                double p = expfracs_all[i][j];
                if (p == 1.0){
                    expfracs_all[i][j] -= err_rate;
                }
                else if (p == 0.0){
                    expfracs_all[i][j] += err_rate;
                }
            }
        }
    }
    
    if (props_given){
        // We've already printed everything, unless summarizing by gene.
        if (genes){
            for (map<string, double>::iterator gs = genesums.begin(); gs != genesums.end(); ++gs){
                double mean = gs->second/genecounts[gs->first];
                string name = gs->first;
                if (snp_gene_names.count(gs->first) > 0){
                    name = snp_gene_names[gs->first];
                }
                fprintf(stdout, "%s\t%s\t%f\t%f\n", gs->first.c_str(), name.c_str(), mean, genecounts[gs->first]);
            }
        }
        // Nothing left to do
        return 0;
    }

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
    else{
        vector<double> props;
        double ll;
        vector<double> dirichlet_params;
        
        infer_props(n, k, expfracs_all, n_trials, samples.size(), 
            ll, props, dirichlet_params, bootstrap, nthread);
        
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
