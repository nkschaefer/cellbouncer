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
#include <regex>
#include <math.h>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/mex.h>
#include "common.h"
#include "ambient_rna.h"
#include "ambient_rna_gex.h"
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
    fprintf(stderr, "    --num_threads -T The number of parallel threads to use when running\n");
    fprintf(stderr, "       optimization problems (default = 1 = no parallelization)\n");
    fprintf(stderr, "    --bootstrap -b Number of bootstrap replicates to run in order to get\n");
    fprintf(stderr, "       variance on mixture proportions of individuals in the ambient RNA\n");
    fprintf(stderr, "       pool. Default = 100.\n");
    fprintf(stderr, "    --doublet_rate -D Expected probability of doublet droplets (for re-IDing\n");
    fprintf(stderr, "       cells). Default = no expectation. Note that this parameter differs from the\n");
    fprintf(stderr, "       one in demux_vcf: in that program, 0.5 = effectively no prior. In this\n");
    fprintf(stderr, "       program, by default no assumption is made about relative frequencies\n");
    fprintf(stderr, "       of different types of singlets and doublets. If you set this parameter,\n");
    fprintf(stderr, "       however, it will compute the overall frequency of each individual in the\n");
    fprintf(stderr, "       data set (as if bulk), and then use this parameter to determine the expected\n");
    fprintf(stderr, "       frequency of each identity (i.e. if ID1 is 5%% of bulk data, ID2 is 10%%, and\n");
    fprintf(stderr, "       D = 0.1, then expected ID1 singlets are 0.9*0.05 + 0.1*0.05*0.05, ID2 singlets\n");
    fprintf(stderr, "       are 0.9*0.1*0.1 + 0.1*0.1*0.1, and ID1+ID2 doublets are 2*0.1*0.05*0.1.\n");
    fprintf(stderr, "       It will then adjust LLRs to encourage identifying the expected proportion of\n");
    fprintf(stderr, "       each identity. This is most useful in high-contamination data sets where\n");
    fprintf(stderr, "       contamination throws off IDs. If you see many more doublets than expected, \n");
    fprintf(stderr, "       set this parameter; if unsure, ignore.\n");
    fprintf(stderr, "    --run_once -r Standard behavior is to iteratively estimate contam profile\n");
    fprintf(stderr, "       and use it to update cell-individual assignments, then repeat until\n");
    fprintf(stderr, "       log likelihood converges. With this option, it will do this process\n");
    fprintf(stderr, "       once and exit.\n");
    fprintf(stderr, "    --ids -i If you limited the individuals to assign when running demux_vcf\n");
    fprintf(stderr, "       (i.e. your VCF contained extra individuals not in the experiment),\n");
    fprintf(stderr, "       provide the filtered list of individuals here. Should be a text file\n");
    fprintf(stderr, "       with one individual name per line, matching individual names in the VCF.\n");
    fprintf(stderr, "    --ids_doublet -I Similar to --ids/-i argument above, but allows control\n");
    fprintf(stderr, "       over which doublet identities are considered. Here, you can specify\n");
    fprintf(stderr, "       individuals and combinations of two individuals to allow. Doublet\n");
    fprintf(stderr, "       combinations not specified in this file will not be considered.\n");
    fprintf(stderr, "       Names of individuals must match those in the VCF, and combinations\n");
    fprintf(stderr, "       of two individuals can be specified by giving both individual names\n");
    fprintf(stderr, "       separated by \"+\", with names in either order.\n");
    fprintf(stderr, "    --dump_freqs -d After inferring the ambient RNA profile, write a file containing\n");
    fprintf(stderr, "        alt allele frequencies at each type of site in ambient RNA. This file will\n");
    fprintf(stderr, "        be called [output_prefix].contam.dat\n");
    fprintf(stderr, "    --llr -l Log likelihood ratio cutoff to filter assignments from demux_vcf.\n");
    fprintf(stderr, "        This is the fourth column in the .assignments file. Default = 0 (no filter)\n");
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
    fprintf(stderr, "    --n_mixprop_trials -N Mixture proportion inference is influenced by initial\n");
    fprintf(stderr, "        guesses. The first time they are inferred, the starting proportions will be\n");
    fprintf(stderr, "        randomly shuffled a number of times equal to this number times the number of\n");
    fprintf(stderr, "        mixture components. Default = 10.\n");
    fprintf(stderr, "    --no_weights -w By default, all observations are weighted by confidence: the log\n");
    fprintf(stderr, "        likelihood ratio of individual ID, divided by the sum of all log likelihood\n");
    fprintf(stderr, "        ratios of assignments of cells to the same individual. This option disables\n");
    fprintf(stderr, "        this weighting and lets all cells contribute equally (although cells with higher\n");
    fprintf(stderr, "        counts will have a stronger influence on the likelihood). You might want to\n");
    fprintf(stderr, "        disable weighting if, for example, you have very unequal numbers of different\n");
    fprintf(stderr, "        individuals and are worried some individual assignments might be mostly noise.\n");
    fprintf(stderr, "        Default behavior would be to give all cells assigned to the noise individual the\n");
    fprintf(stderr, "        same overall weight as all cells assigned to any other individual.\n"); 
    print_libname_help();
    fprintf(stderr, "===== OPTIONAL; FOR INFERRING GENE EXPRESSION =====\n");
    fprintf(stderr, "    --barcodes -B (Optionally gzipped) barcodes file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --features -F (Optionally gzipped) features file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --matrix -M (Optionally gizpped) matrix file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --feature_type -t (OPTIONAL) If --features/-f contains more than one type of data\n");
    fprintf(stderr, "        (i.e. gene expression and feature barcoding), use this to specify which feature\n");
    fprintf(stderr, "        type is RNA-seq (for 10X Genomics, \"Gene Expression\"). By default, includes all\n");
    fprintf(stderr, "        features and does not check.\n");
    fprintf(stderr, "    --clusts -c (RECOMMENDED) cell-cluster assignments computed by another program.\n");
    fprintf(stderr, "    --noround -R By default, decontaminated counts are rounded to the nearest integer, in\n");
    fprintf(stderr, "        a random fashion so that the appropriate number of bulk counts are removed. This\n");
    fprintf(stderr, "        satisfies the requirements of some differential expression tools for integer counts.\n");
    fprintf(stderr, "        Setting this option will instead allow unrounded decimal counts to be output.\n");
    fprintf(stderr, "    --skip_genes -g Provide a list of gene names, one per line, to exclude from considering\n");
    fprintf(stderr, "        as part of the contamination profile. Default = keep all genes.\n");
    fprintf(stderr, "    --skip_genes_regex -G A string of pipe-separated regex-style gene name matching strings\n");
    fprintf(stderr, "        to use to exclude genes from ambient RNA removal. Since ambient RNA is usually detected\n");
    fprintf(stderr, "        from diploid genomic variants, for example, this can be used to exclude mitochondrial\n");
    fprintf(stderr, "        genes (which were not included in the inference). Default = \"^MT-\"\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int parse_clustfile(string& filename, 
    robin_hood::unordered_map<unsigned long, int>& clusts,
    vector<string>& clustnames){
    
    ifstream inf(filename.c_str());
    string bcstr;
    string clustn;
    set<string> clustnames_set;
    robin_hood::unordered_map<unsigned long, string> ctmp;
    while (inf >> bcstr >> clustn){
        clustnames_set.insert(clustn);
        unsigned long ul = bc_ul(bcstr);
        ctmp.emplace(ul, clustn);
    }

    map<string, int> clust2idx;
    int idx = 0;
    for (set<string>::iterator s = clustnames_set.begin(); s != clustnames_set.end(); ++s){
        clust2idx.insert(make_pair(*s, idx));
        clustnames.push_back(*s);
        idx++;
    }
    for (robin_hood::unordered_map<unsigned long, string>::iterator t = ctmp.begin(); 
        t != ctmp.end(); ++t){
        clusts.emplace(t->first, clust2idx[t->second]);
    }
    return clustnames.size();
}

void infer_from_genotypes(string& output_prefix,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    string& idfile,
    bool idfile_given,
    string& idfile_doublet,
    bool idfile_doublet_given,
    vector<string>& samples,
    map<int, double>& contam_prof,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    double doublet_rate,
    int num_threads,
    double error_ref,
    double error_alt,
    bool inter_species,
    int n_mixprop_trials,
    bool weight,
    bool run_once,
    int bootstrap,
    string& libname,
    bool seurat,
    bool cellranger,
    bool underscore){
    
    // Load conditional matching probabilities
    map<pair<int, int>, map<int, float> > exp_match_fracs;
    string expfrac_name = output_prefix + ".condf";
    if (file_exists(expfrac_name)){
        load_exp_fracs(expfrac_name, exp_match_fracs);
    }
    else{
        fprintf(stderr, "ERROR: no conditional matching probability file found for %s.\n", 
            output_prefix.c_str());
        fprintf(stderr, "Please re-run demux_vcf with the same VCF file and output prefix, \
but specify the -F option to create this file. Then re-run this program.\n");
        exit(1); 
    }
    
    // Load filtered ID list(s), if given
    set<int> allowed_ids;
    set<int> allowed_ids2;

    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile.c_str());
        }
    }
    if (idfile_doublet_given){
        parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile_doublet.c_str());
        }
    }

    // Load stored allele counts
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    string counts_name = output_prefix + ".counts";
    if (file_exists(counts_name)){
        fprintf(stderr, "Loading counts...\n");
        load_counts_from_file(indv_allelecounts, samples, counts_name, allowed_ids); 
    }
    else{
        fprintf(stderr, "ERROR: no counts found for %s. Please run demux_vcf with same\n",
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1);
    }
    
    double llprev = 0.0;
    double delta = 999;
    double delta_thresh = 0.1;
    
    map<int, double> contam_prof_conc;
    robin_hood::unordered_map<unsigned long, double> contam_rate_se;
    
    int nits = 0;
    while (delta > delta_thresh){
        fprintf(stderr, "===== ITERATION %d =====\n", nits+1);
        contamFinder cf(indv_allelecounts, assn, assn_llr, exp_match_fracs, samples.size(),
            allowed_ids, allowed_ids2);
        cf.set_doublet_rate(doublet_rate);
        cf.set_num_threads(num_threads);
        if (run_once){
            //cf.no_reassign();
        }

        // Initialize to whatever was the final estimate last time 
        if (nits > 0){
            cf.set_init_contam_prof(contam_prof);
        }
        if (nits > 0){
            double meanc = 0.0;
            double cfrac = 1.0/(double)contam_rate.size();
            for (robin_hood::unordered_map<unsigned long, double>::iterator c = contam_rate.begin();
                c != contam_rate.end(); ++c){
                meanc += cfrac * c->second;
            }
            cf.set_init_c(meanc);
        } 
        // Do standard initialization
        cf.set_error_rates(error_ref, error_alt);
        if (inter_species){
            cf.model_other_species();
        } 
        cf.set_mixprop_trials(n_mixprop_trials);
        if (weight){
            cf.use_weights();
        }
        cf.fit(); 
        
        for (int i = 0; i < samples.size(); ++i){
            fprintf(stderr, "%s) %f\n", samples[i].c_str(), cf.contam_prof[i]);
        }   

        double ll = cf.compute_ll();
        
        if (run_once){
            // Break out of cycle
            assn = cf.assn;
            assn_llr = cf.assn_llr;
            contam_prof = cf.contam_prof;
            contam_rate = cf.contam_rate;
            contam_rate_se = cf.contam_rate_se;
            
            delta = 0;
        }
        else{
            if (llprev == 0 || ll > llprev){
                assn = cf.assn;
                assn_llr = cf.assn_llr;
                contam_prof = cf.contam_prof;
                contam_rate = cf.contam_rate;
                contam_rate_se = cf.contam_rate_se;
            }
            fprintf(stderr, " -- Log likelihood: %f", ll);
            if (llprev != 0){
                delta = ll-llprev;
                fprintf(stderr, " delta = %f\n", delta);
            }
            else{
                fprintf(stderr, "\n");
            }
            llprev = ll;
            nits++;
        }
        if (delta <= delta_thresh && bootstrap > 0){
            // Do bootstrapping
            cf.assn = assn;
            cf.assn_llr = assn_llr;
            cf.contam_prof = contam_prof;
            cf.contam_rate = contam_rate;
            cf.contam_rate_se = contam_rate_se;
            fprintf(stderr, "Computing Dirichlet concentration parameters \
on mixture proportions...\n");
            cf.bootstrap_amb_prof(bootstrap, contam_prof_conc);
        }
    }

    // Write contamination profile to disk
    {
        string fname = output_prefix + ".contam_prof";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing contamination profile to disk...\n");
        dump_contam_prof(outf, contam_prof, contam_prof_conc, samples);
        fclose(outf);
    }
    // Write contamination rate (and standard error) per cell to disk
    {
        string fname = output_prefix + ".contam_rate";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_contam_rates(outf, contam_rate, contam_rate_se, samples,
            libname, cellranger, seurat, underscore);
        fclose(outf);
    }
    /*
    if (dump_freqs){
        // Write alt allele matching frequencies to a file
        string fname = output_prefix + ".contam.dat";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_amb_fracs(outf, cf.amb_mu);
        fclose(outf);
    }
    */
    // Write refined assignments to disk
    // If run_once is set, assume the user is happy with the preliminary, un-refined
    // assignments.
    if (true){
    //if (!run_once){
        string fname = output_prefix + ".decontam.assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_assignments(outf, assn, assn_llr, samples, libname, 
            cellranger, seurat, underscore);
        fclose(outf); 
    }
    else if (file_exists(output_prefix + ".decontam.assignments")){
        string fn = output_prefix + ".decontam.assignments";
        unlink(fn.c_str());
    }
}

void parse_rates(string& filename, robin_hood::unordered_map<unsigned long, double>& contam_rates){
    ifstream inf(filename);
    string bc_str;
    double rate;
    double rate_se;
    while (inf >> bc_str >> rate >> rate_se){
        unsigned long bc_key = bc_ul(bc_str);
        contam_rates.emplace(bc_key, rate);
    }
}

void parse_prof(string& filename, map<int, double>& contam_prof, vector<string>& samples){
    
    // Note: sometimes contam_prof can contain an "other_species" entry, keyed to -1.
    // For this purpose, we are loading contam prof to help guide the GEX profiling.
    // That can't make use of individuals that aren't in the data set, so we will 
    // exclude this individual if it exists and normalize the profile so it sums to 1.

    map<string, int> samp2idx;
    for (int i = 0; i < samples.size(); ++i){
        samp2idx.insert(make_pair(samples[i], i));
    }

    ifstream inf(filename);
    string line;
    double proptot = 0.0;
    while(getline(inf, line)){
        istringstream splitter(line);
        string field;
        int fld_idx = 0;
        int samp_idx = -1;
        double samp_prop = 0.0;
        // Last column, if it exists, is Dirichlet concentration
        // parameters (ignore)
        while(getline(splitter, field, '\t')){
            if (fld_idx == 0){
                if (samp2idx.count(field) > 0){
                    samp_idx = samp2idx[field];
                }
            }
            else if (fld_idx == 1){
                samp_prop = atof(field.c_str());
            }
            ++fld_idx;
        }

        if (samp_idx != -1){
            contam_prof.insert(make_pair(samp_idx, samp_prop));
            proptot += samp_prop;
        }
    }
    for (map<int, double>::iterator p = contam_prof.begin(); p != contam_prof.end();
        ++p){
        p->second /= proptot;
    }
}

void parse_skip_genes(string& skipgenesfile, set<string>& skipgenes){
    ifstream inf(skipgenesfile);
    string line;
    while (inf >> line){
        fprintf(stderr, "Skipping gene (raw expr output): %s\n", line.c_str());
        skipgenes.insert(line);
    }
}

void process_gex_data(string& output_prefix,
    string& barcodesfile,
    string& featuresfile,
    string& matrixfile,
    string& feature_type,
    string& clustfile,
    bool idfile_doublet_given,
    robin_hood::unordered_map<unsigned long, int>& assn,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    map<int, double>& contam_prof,
    bool round,
    int num_threads,
    string& libname,
    bool seurat,
    bool cellranger,
    bool underscore,
    string& skipgenesfile,
    string& skip_genes_regex){
    
    robin_hood::unordered_map<unsigned long, map<int, long int> > mtx;
    vector<string> features;
    fprintf(stderr, "Loading gene expression data...\n");
    bool success = parse_mex(barcodesfile, featuresfile, matrixfile, mtx, features, feature_type);
    if (!success){
        exit(1);
    }
    robin_hood::unordered_map<unsigned long, int> clusts;
    int nclusts = 0;
    vector<string> clustnames;
    if (clustfile != ""){
        nclusts = parse_clustfile(clustfile, clusts, clustnames);
    }
    else{
        fprintf(stderr, "Using cell identities as clusters\n");
        
        if (!idfile_doublet_given){
            clustnames = samples;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                if (a->second < samples.size()){
                    clusts.emplace(a->first, a->second);
                }
            }
            nclusts = samples.size();
        }
        else{
            // A file listing specific allowed doublet identities was given.
            // Treat doublets as unique identities rather than combinations of
            // singlets.

            // This is to make cluster IDs sorted alphabetically.
            set<pair<string, int> > name_set;
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                name_set.insert(make_pair(idx2name(a->second, samples), a->second));
            }
            nclusts = name_set.size();
            int idx = 0;
            map<int, int> idx2idx;
            for (set<pair<string, int> >::iterator n = name_set.begin(); n != name_set.end(); ++n){
                idx2idx.insert(make_pair(n->second, idx));
                clustnames.push_back(n->first);
                ++idx;
            }
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                clusts.emplace(a->first, idx2idx[a->second]);
            }
        }
    }
    contam_profiler_gex contam_profiler(contam_rate, contam_prof, 
        assn, samples.size(), idfile_doublet_given);
    contam_profiler.set_threads(num_threads);
    contam_profiler.set_mtx(mtx, features.size());
    contam_profiler.set_clusts(clusts, nclusts);
    if (round){
        contam_profiler.round_counts();
    }
    if (skipgenesfile != "" || skip_genes_regex != ""){
        set<string> skip_genes_txt;
        if (skipgenesfile != ""){
            parse_skip_genes(skipgenesfile, skip_genes_txt);
        }
        // Also search for gene names using regex
        if (skip_genes_regex != ""){
            const std::regex skip_regex(skip_genes_regex);
            for (int i = 0; i < features.size(); ++i){
                smatch matches;
                if (regex_search(features[i], matches, skip_regex)){
                    fprintf(stderr, "Skipping gene (raw expr output): %s\n", features[i].c_str());
                    skip_genes_txt.insert(features[i]);         
                } 
            }
        }

        if (skip_genes_txt.size() > 0){
            fprintf(stderr, "To avoid skipping genes, omit -g and set option -G \"\"\n");
        }
        // Map to int
        map<string, int> gene2idx;
        for (int i = 0; i < features.size(); ++i){
            gene2idx.insert(make_pair(features[i], i));
        }
        set<int> skipgenes;
        for (set<string>::iterator sgt = skip_genes_txt.begin(); sgt != skip_genes_txt.end(); 
            ++sgt){
            skipgenes.insert(gene2idx[*sgt]);
        }
        contam_profiler.skip_genes(skipgenes);
    }
    
    // Infer ambient RNA profile
    contam_profiler.get_profile();
    
    // Write to disk
    {
        string fname = output_prefix + ".gex_profile";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(outf, "gene\tambient");
        for (int i = 0; i < nclusts; ++i){
            fprintf(outf, "\t%s", clustnames[i].c_str());
        }
        fprintf(outf, "\n");
        for (int i = 0; i < features.size(); ++i){
            fprintf(outf, "%s\t%e", features[i].c_str(), contam_profiler.prof_ambient[i]);
            for (int j = 0; j < nclusts; ++j){
                fprintf(outf, "\t%e", contam_profiler.prof_clusts[j][i]);
            }
            fprintf(outf, "\n");
        }
        fclose(outf);
    }

    // Remove ambient RNA profile
    contam_profiler.decontam();

    // Write cleaned up matrix to disk
    string out_mtx = output_prefix + "_mtx";
    write_mex(out_mtx, contam_profiler.mtx_decontam, 
        features, round, libname, cellranger, seurat, underscore);
}

int main(int argc, char *argv[]) {    

    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"other_species", required_argument, 0, 's'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"llr", required_argument, 0, 'l'},
       {"n_mixprop_trials", required_argument, 0, 'N'},
       {"no_weights", no_argument, 0, 'w'},
       {"dump_freqs", no_argument, 0, 'd'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"run_once", no_argument, 0, 'r'},
       {"bootstrap", required_argument, 0, 'b'},
       {"barcodes", required_argument, 0, 'B'},
       {"features", required_argument, 0, 'F'},
       {"matrix", required_argument, 0, 'M'},
       {"feature_type", required_argument, 0, 't'},
       {"clusts", required_argument, 0, 'c'},
       {"skip_genes", required_argument, 0, 'g'},
       {"skip_genex_regex", required_argument, 0, 'G'},
       {"num_threads", required_argument, 0, 'T'},
       {"noround", no_argument, 0, 'R'},
       {0, 0, 0, 0} 
       
    };
    
    // Set default values
    string output_prefix = "";
    bool inter_species = false;
    double error_ref = 0.001;
    double error_alt = 0.001;
    double llr = 0.0;
    int n_mixprop_trials = 10;
    bool weight = true;
    bool dump_freqs = false;
    string idfile;
    bool idfile_given = false;
    string idfile_doublet;
    bool idfile_doublet_given = false;
    string libname = "";
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;
    bool run_once = false;
    int bootstrap = 100;
    double doublet_rate = -1.0;
    int num_threads = 0;
    string skipgenesfile = "";
    
    string skip_genes_regex = R"(^MT-)";

    string barcodesfile = "";
    string featuresfile = "";
    string matrixfile = "";
    string feature_type = "";
    string clustfile = "";
    bool round = true;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:e:g:G:E:l:N:i:I:n:b:D:B:F:M:t:c:T:RrsCSUdwh", long_options, &option_index )) != -1){
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
            case 'g':
                skipgenesfile = optarg;
                break;
            case 'G':
                skip_genes_regex = optarg;
                break;
            case 'i':
                idfile = optarg;
                idfile_given = true;
                break;
            case 'I':
                idfile_doublet = optarg;
                idfile_doublet_given = true;
                break;
            case 's':
                inter_species = true;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 'l':
                llr = atof(optarg);
                break;
            case 'N':
                n_mixprop_trials = atoi(optarg);
                break;
            case 'w':
                weight = false;
                break;
            case 'd':
                dump_freqs = true;
                break;
            case 'n':
                libname = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'S':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case 'r':
                run_once = true;
                break;
            case 'b':
                bootstrap = atoi(optarg);
                break;
            case 'B':
                barcodesfile = optarg;
                break;
            case 'F':
                featuresfile = optarg;
                break;
            case 'M':
                matrixfile = optarg;
                break;
            case 't':
                feature_type = optarg;
                break;
            case 'c':
                clustfile = optarg;
                break;
            case 'R':
                round = false;
                break;
            case 'T':
                num_threads = atoi(optarg);
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
    if (n_mixprop_trials < 0){
        fprintf(stderr, "ERROR: --n_mixprop_trials must be >= 0\n");
        exit(1);
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i and -I is allowed.\n");
        exit(1);
    }
    if (bootstrap <= 0){
        fprintf(stderr, "WARNING: bootstrapping disabled. Ambient RNA pool proportions will \
be reported without concentration parameters (variance will be unknown).\n");
    }
    if (doublet_rate != -1 && (doublet_rate < 0 || doublet_rate > 1)){
        fprintf(stderr, "ERROR: --doublet_rate/-D must be between 0 and 1, inclusive.\n");
        exit(1);
    }
    if ((barcodesfile != "" || featuresfile != "" || matrixfile != "") &&
        !(barcodesfile != "" && featuresfile != "" && matrixfile != "")){
        fprintf(stderr, "ERROR: if inferring gene expression profile, you must provide all\n");
        fprintf(stderr, "three of --barcodes/-B, --features/-F, and --matrix/-M\n");
        exit(1);
    }
    if (barcodesfile != "" && clustfile == ""){
        fprintf(stderr, "WARNING: inferring expression profile of contamination without \
cluster information. Assuming one default expression profile for each individual (results may \
be less accurate if there is much cell type heterogeneity).\n");
    }
    if (clustfile != "" && barcodesfile == ""){
        fprintf(stderr, "ERROR: --clusters/-c only applicable when loading gene expression data\n");
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
    
    // 1 thread means 0 threads (don't launch extra processes)
    if (num_threads <= 1){
        num_threads = 0;
    }
    
    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    string assn_name = output_prefix + ".assignments";
    if (file_exists(assn_name)){
        fprintf(stderr, "Loading assignments...\n");
        load_assignments_from_file(assn_name, assn, assn_llr, samples);
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
    
    // Obtain contamination profile (from allele matching data)
    string prof_name = output_prefix + ".contam_prof";
    string rate_name = output_prefix + ".contam_rate";

    robin_hood::unordered_map<unsigned long, double> contam_rate;
    map<int, double> contam_prof;
    
    bool load_gex = (barcodesfile != "" && featuresfile != "" && matrixfile != "");
     
    if (file_exists(prof_name) && file_exists(rate_name)){
        if (!load_gex){
            fprintf(stderr, "ERROR: previous run detected, and no gene expression data \
provided. Nothing to do.\n");
            fprintf(stderr, "To repeat the run, remove %s and %s.\n", prof_name.c_str(), 
                rate_name.c_str());
            exit(1);
        }
        fprintf(stderr, "Previous run detected. Loading profile\n");
        fprintf(stderr, "To prevent this behavior, change --output_prefix/-o or remove \
%s and %s.\n", prof_name.c_str(), rate_name.c_str());
        parse_rates(rate_name, contam_rate);
        parse_prof(prof_name, contam_prof, samples);
    }
    else{
        infer_from_genotypes(output_prefix,
            assn,
            assn_llr,
            idfile,
            idfile_given,
            idfile_doublet,
            idfile_doublet_given,
            samples,
            contam_prof,
            contam_rate,
            doublet_rate,
            num_threads,
            error_ref,
            error_alt,
            inter_species,
            n_mixprop_trials,
            weight,
            run_once,
            bootstrap,
            libname,
            seurat,
            cellranger,
            underscore);
    }
    
    if (load_gex){
        process_gex_data(output_prefix,
            barcodesfile,
            featuresfile,
            matrixfile,
            feature_type,
            clustfile,
            idfile_doublet_given,
            assn,
            samples,
            contam_rate,
            contam_prof,
            round,
            num_threads,
            libname,
            seurat,
            cellranger,
            underscore,
            skipgenesfile,
            skip_genes_regex);
    }
    return 0;
}
