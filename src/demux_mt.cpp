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
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <deque>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <zlib.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <mixtureDist/functions.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <optimML/brent.h>
#include "common.h"

/**
 * Given a BAM file of mapped single-cell sequencing data (ATAC or 
 * possibly RNA-seq), attempts to find SNPs on the mitochondrial
 * sequence, cluster these into haplotypes, and then separate
 * cells based on these haplotypes. Does not require the user 
 * to know anything about numbers of individuals or their
 * genotypes a priori. Assumes an approximately equal mixture
 * of all individuals.
 */

using std::cout;
using std::endl;
using namespace std;

typedef bitset<MAX_SITES> hapstr;

int nvars = -1;

// Data structure to store information about a potential variant
// site, before haplotypes are created. Stores all relevant 
// information about the site, minus which cells have each of its
// variants.
struct varsite{
    int pos;
    char allele1;
    char allele2;
    int cov;
    float freq1;
    float freq2;
    
    varsite(){
        this->pos = -1;
        this->allele1 = 'N';
        this->allele2 = 'N';
        this->cov = 0;
        this->freq1 = 0.0;
        this->freq2 = 0.0;
    }; 
    varsite(const varsite& v){
        this->pos = v.pos;
        this->allele1 = v.allele1;
        this->allele2 = v.allele2;
        this->cov = v.cov;
        this->freq1 = v.freq1;
        this->freq2 = v.freq2;
    };
};
// Data structure to store counts at a SNP site in an individual cell.
// counts1 stores reads matching the major allele, counts2 stores
// reads matching the minor allele.
struct var_counts{
    int counts1[MAX_SITES];
    int counts2[MAX_SITES];
    var_counts(){
        for (int i = 0; i < nvars; ++i){
            this->counts1[i] = 0;
            this->counts2[i] = 0;
        }
    };
    var_counts(const var_counts& vc){
        for (int i = 0; i < nvars; ++i){
            this->counts1[i] = vc.counts1[i];
            this->counts2[i] = vc.counts2[i];
        }
    };
};

// Collapsed, easier to compare version of above -- for each cell, 
// store a bitset where every element is a variable site on the
// mitochondrion. mask is 0 or 1 to denote whether or not that
// site was covered/visible in this cell, and vars stores whether
// the site (if covered) matched the major allele (0) or minor
// allele (1).
struct hap{
    hapstr vars;
    hapstr mask;
    
    hap(){};
    hap(const hapstr& m, const hapstr& v){
        this->mask = m;
        this->vars = v;
    };
    hap(const hap& h){
        this->vars = h.vars;
        this->mask = h.mask;
    };
    bool operator==(const hap& h){
        return this->vars == h.vars && this->mask == h.mask;
    };
};

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_mt [OPTIONS]\n");
    fprintf(stderr, "Given a BAM file representing single-cell sequencing data \n");
    fprintf(stderr, "(RNA-seq or ATAC-seq with cell barcode/CB tags present), \n");
    fprintf(stderr, "Attempts to find variable sites on the mitochondrial genome and \n");
    fprintf(stderr, "then cluster these into mitochondrial haplotypes.\n");
    fprintf(stderr, "Attempts to infer the correct number of haplotypes in the mixture\n");
    fprintf(stderr, "(assuming they are close to equally represented).\n");
    fprintf(stderr, "Outputs the list of variable sites used, the haplotypes at these\n");
    fprintf(stderr, "sites inferred for each cluster, the haplotype at each site for\n");
    fprintf(stderr, "every cell in the data set, and the assignment of every cell to its\n");
    fprintf(stderr, "most likely mitochondrial haplotype (optionally including doublet\n");
    fprintf(stderr, "combinations).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Options for all run modes =====\n");
    fprintf(stderr, "   ---------- I/O Options ----------\n");
    fprintf(stderr, "   --output_prefix -o The base file name, including path, for output\n");
    fprintf(stderr, "       files. If inferring haplotypes, will create files with the \n");
    fprintf(stderr, "       extensions .vars, .haps, .cellhaps, .assignments, and .summary.\n");
    fprintf(stderr, "       If assigning cells from pre-computed haplotypes, will only\n");
    fprintf(stderr, "       create .cellhaps, .assignments, and .summary. (REQUIRED)\n");
    fprintf(stderr, "   --barcodes -B To only print cell-haplotype assignments for a subset\n");
    fprintf(stderr, "       of all cell barcodes in the BAM file (i.e. the list of barcodes\n");
    fprintf(stderr, "       deemed to represent legitimate cells), provide a file listing\n");
    fprintf(stderr, "       barcodes in the subset here. If using CellRanger output, this\n");
    fprintf(stderr, "       would be in filtered_feature_bc_matrix/barcodes.tsv or\n");
    fprintf(stderr, "       barcodes.tsv.gz. NOTE: all barcodes in the BAM will still be used\n");
    fprintf(stderr, "       in clustering. To limit which barcodes are used to find variants\n");
    fprintf(stderr, "       and cluster, use the -f option.\n");
    print_libname_help();
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Options for run modes 1 and 2 (clustering & identifying haploytpes) =====\n");
    fprintf(stderr, "   ---------- I/O Options ----------\n");
    fprintf(stderr, "   --dump -d Prints information about allelic state of cell barcodes \n");
    fprintf(stderr, "       at each site (.bchaps file) and quits. If in clustering mode, \n");
    fprintf(stderr, "       this happens after finding variant sites and before clustering\n");
    fprintf(stderr, "       (alternatively, .bchaps will be automatically created after\n");
    fprintf(stderr, "       clustering with only the sites found relevant for clustering).\n");
    fprintf(stderr, "       If loading previously-inferred clusters, prints allelic state \n");
    fprintf(stderr, "       in every cell at all sites in the .vars file. Bypasses the step\n");
    fprintf(stderr, "       of assigning cells to mitochondrial haplotypes.\n");
    fprintf(stderr, "   ---------- General Options ----------\n");
    fprintf(stderr, "   --doublet_rate -D A decimal between 0 and 1 representing the prior\n");
    fprintf(stderr ,"       probability of a cell being a doublet. Set to 0 to disable doublet\n");
    fprintf(stderr, "       identification altogether. Default = 0.5\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Run mode 1: Inferring clusters from a BAM file =====\n");
    fprintf(stderr, "   This mode is used to infer mitochondrial haplotypes from a BAM file, and \n");
    fprintf(stderr, "       then assign cells to the most likely inferred haplotype.\n");
    fprintf(stderr, "   ---------- I/O options ----------\n");
    fprintf(stderr, "   --bam -b The BAM file containing the data to use. (REQUIRED)\n");
    fprintf(stderr, "   --barcodes_filter -f Distinct from the -B option (above), which limits\n");
    fprintf(stderr, "       which barcodes will be assigned mitochondrial haplotypes. This\n");
    fprintf(stderr, "       option limits which barcodes are used to find variants and cluster\n");
    fprintf(stderr, "       haplotypes. This option should only be used if you wish to find sub-\n");
    fprintf(stderr, "       clusters within a specific haplotype from a previous run: pull out all\n");
    fprintf(stderr, "       barcodes assigned to the cluster of interest and provide them here as\n");
    fprintf(stderr, "       a text file, one per line. Optionally, can be gzipped.\n"); 
    fprintf(stderr, "   --mito -m The name of the mitochondrial sequence in the reference\n");
    fprintf(stderr, "       genome (OPTIONAL; default: chrM)\n");
    fprintf(stderr, "   ---------- Filtering options ----------\n");
    fprintf(stderr, "   --no_cov_filt -c Do not remove variable sites with low coverage across all\n");
    fprintf(stderr, "       cells. Default behavior is to compute a histogram of sites passing filter\n");
    fprintf(stderr, "       with increasing coverage cutoffs, then find the knee point (defined as \n");
    fprintf(stderr, "       the point of maximum curvature), ensuring at least 25%% of all sites are\n");
    fprintf(stderr, "       included. Disabling this might be appropriate for low-coverage data sets,\n");
    fprintf(stderr, "       especially scRNA-seq, where coverage is expected to vary from site to site.\n");
    fprintf(stderr, "   --mapq -q The minimum map quality filter (OPTIONAL; default 20)\n");
    fprintf(stderr, "   --baseq -Q The minimum base quality filter (OPTIONAL; default 20)\n");
    fprintf(stderr, "   ---------- General options ----------\n"); 
    fprintf(stderr, "   --nclustmax -N (OPTIONAL) If you know how many individuals should\n");
    fprintf(stderr, "       be in the mixture, you can specify that number here, and only\n");
    fprintf(stderr, "       up to that many haplotypes will be inferred. It's possible that\n");
    fprintf(stderr, "       fewer than that number of haplotypes will be found, however.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Run mode 2: Assigning cells to previously-inferred haplotypes =====\n");
    fprintf(stderr, "   This mode is used when you have already inferred mitochondrial haplotypes (using\n");
    fprintf(stderr, "       run mode 1), and want to assign cells (i.e. from a different BAM file) to\n");
    fprintf(stderr, "       the most likely haplotype, out of the set of haplotypes you already inferred.\n");
    fprintf(stderr, "   ---------- I/O Options ----------\n");
    fprintf(stderr, "   --haps -H Cluster haplotypes from a previous run. Should be that\n");
    fprintf(stderr, "       run's [output_prefix].haps. (REQUIRED)\n");
    fprintf(stderr, "   --vars -v Variants used to build cluster haplotypes in previous \n");
    fprintf(stderr, "       run. Should be that run's [output_prefix].vars. (REQUIRED)\n");
    fprintf(stderr, "   --ids -i (OPTIONAL). A file listing names of IDs to assign, one per line, in the\n");
    fprintf(stderr, "       same order as the haplotypes in the --haps file. If you do not provide this\n");
    fprintf(stderr, "       file, individual names will be 0-based numeric IDs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Run mode 3: Inferring mixing proportions in cells with multiple haplotypes =====\n");
    fprintf(stderr, "   This is a special use case: mitochondrial haplotypes are already known, some\n");
    fprintf(stderr, "       cells are known to contain more than one mitochondrial haplotype (i.e. \n");
    fprintf(stderr, "       tetraploid composite cell lines), and you want to infer the mixing proportion\n");
    fprintf(stderr, "       of each mitochondrial haplotype within each cell. Inferred mixing proportions\n");
    fprintf(stderr, "       will be written to [output_prefix].props.\n");
    fprintf(stderr, "   ---------- I/O Options ----------\n");
    fprintf(stderr, "   --haps -H Cluster haplotypes from a previous run. Should be that\n");
    fprintf(stderr, "       run's [output_prefix].haps. (REQUIRED)\n");
    fprintf(stderr, "   --vars -v Variants used to build cluster haplotypes in previous \n");
    fprintf(stderr, "       run. Should be that run's [output_prefix].vars. (REQUIRED)\n");
    fprintf(stderr, "   --assignments -a A file mapping cell barcode to identity (could be an assignment \n");
    fprintf(stderr, "       file from demux_mt or demux_vcf, but only the first two columns are required.\n");
    fprintf(stderr, "       Only cells with multiple mitochondrial haplotypes (can be designated in either\n");
    fprintf(stderr, "       order and separated by + or x) will be considered. (REQUIRED)\n");
    exit(code);
}

// ===== Functions for pileup, to find SNPs on the mitochondrion =====

// Info to pass to function (below) that will be called by pileup function
struct infile_t {
    const char* fname;
    samFile* fp;
    sam_hdr_t* fp_hdr;
    hts_itr_t* itr;
    hts_idx_t* idx;
    
    // Constructor
    infile_t(const char* fname, const char* mitoname){
        this->fname = fname;
        
        // Load/create BAM index
        this->idx = hts_idx_load(fname, HTS_FMT_BAI);
        
        // Retrieve region of interest from BAM index
        
        // First need to get TID of mitochondrion
        this->fp = sam_open(fname, "r");
        this->fp_hdr = sam_hdr_read(this->fp);
        int mito_tid = -1;
        hts_pos_t mito_len;
        for (int i = 0; i < this->fp_hdr->n_targets; ++i){
            if (strcmp(this->fp_hdr->target_name[i], mitoname) == 0){
                mito_tid = i;
                mito_len = sam_hdr_tid2len(this->fp_hdr, i);
                break;
            }
        } 
        if (mito_tid == -1){
            fprintf(stderr, "ERROR: mitochondrial sequence %s not found \
in BAM header\n", mitoname);
            exit(1);
        }

        this->itr = sam_itr_queryi(this->idx, mito_tid, 0, mito_len);  
    };
    
    // Destructor
    ~infile_t(){
        sam_hdr_destroy(this->fp_hdr);
        sam_close(this->fp);
        hts_idx_destroy(this->idx);
    }; 

};

struct infile_bc_wrapper{
    infile_t* infile;
    set<unsigned long>* bcs;
};

// Iterator function for pileup
static int readaln(void *data, bam1_t *b){
    // Retrieve data in usable format
    infile_t *g = (infile_t*)data;
    int ret;
    while (1){
        ret = sam_itr_next(g->fp, g->itr, b); 
        if (ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        break;
    }
    return ret;
}

// same as above, but including barcode whitelist
static int readaln_bcs(void *data, bam1_t *b){
    // Retrieve data in usable format
    infile_bc_wrapper* wrap = (infile_bc_wrapper*)data;
    infile_t* g = wrap->infile;
    set<unsigned long>* bc_whitelist = wrap->bcs;
    int ret;
    while (1){
        ret = sam_itr_next(g->fp, g->itr, b); 
        if (ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        // try to retrieve cell barcode
        uint8_t* bc_bin = bam_aux_get(b, "CB");
        if (bc_bin != NULL){
            // Convert to string
            char* bc_char = bam_aux2Z(bc_bin);
            // Convert to C++ style string
            string bc_str = bc_char;
            int dashpos = bc_str.find("-");
            if (dashpos != string::npos){
                bc_str = bc_str.substr(0, dashpos);
            }
            bc as_bitset;
            str2bc(bc_str.c_str(), as_bitset);
            unsigned long as_ulong = as_bitset.to_ulong();
            if (bc_whitelist->find(as_ulong) == bc_whitelist->end()){
                // Cell barcode not in whitelist
                continue;
            }
        }
        else{
            // No cell barcode
            continue;
        }
        break;
    }
    return ret;
}


/**
 * Given a list of numbers representing sizes of clusters, finds the optimal way
 * to partition these sizes into true clusters and erroneous clusters. Assumes 
 * that true clusters have sizes drawn from the same Poisson distribution, and that
 * erroneous clusters have sizes drawn from a different Poisson distribution with
 * lower mean.
 *
 * Returns the inferred number of true clusters.
 */
int choose_nclust(vector<int>& clust_sizes){

    // Sort values high to low
    int nclust_tot = 0;
    vector<int> clust_sorted;
    for (vector<int>::iterator c = clust_sizes.begin(); c != clust_sizes.end(); 
        ++c){
        clust_sorted.push_back(-*c);
        nclust_tot += *c;
    }
    sort(clust_sorted.begin(), clust_sorted.end());
    
    // What is the maximum possible number of clusters?
    int nclust_max = clust_sorted.size();
    
    // What is the optimal number of clusters? 
    int nclust_best = -1;
    double loglik_max = 0.0;
    
    int nclust_min = 1;
    for (int nclust = nclust_min; nclust <= nclust_max; ++nclust){
        
        // Get mean size of valid clusters under this model
        double meansize1 = 0.0;
        for (int i = 0; i < nclust; ++i){
            meansize1 += (double)-clust_sorted[i];
        }
        meansize1 /= (double)nclust;
        
        // Get mean size of invalid clusters under this model
        double meansize2 = 0.0;
        for (int j = nclust; j < clust_sorted.size(); ++j){
            meansize2 += (double)-clust_sorted[j];
        }
        meansize2 /= (double)(clust_sorted.size()-nclust);
        
        double loglik = 0.0;
        int ncell_this = 0;
        for (int i = 0; i < nclust; ++i){
            loglik += dpois(-clust_sorted[i], meansize1);
            ncell_this += -clust_sorted[i];
        }
        
        for (int i = nclust; i < clust_sorted.size(); ++i){
            loglik += dpois(-clust_sorted[i], meansize2);
        }
        
        //fprintf(stderr, "nclust %d LL %f means %f %f\n", nclust, 
        //    loglik, meansize1, meansize2);

        if (nclust_best == -1 || loglik > loglik_max){
            nclust_best = nclust;
            loglik_max = loglik;
        }
        else if (nclust_best != -1 && loglik < loglik_max){
            break;
        }
    }
    
    return nclust_best;
}

void write_vars(string& mito_chrom, 
    string& output_prefix, 
    deque<varsite>& vars, 
    hapstr& mask_global, 
    bool use_mask){
    
    string vars_out = output_prefix + ".vars";
    FILE* vars_outf = fopen(vars_out.c_str(), "w");
    
    // Just print the non-blacklisted variable sites.
    int idx = 0;
    for (deque<varsite>::iterator v = vars.begin(); v != vars.end(); ++v){
        if (!use_mask || mask_global.test(idx)){
            fprintf(vars_outf, "%s\t%d\t%c\t%c\t%.4f\t%.4f\t%d\n", mito_chrom.c_str(), v->pos + 1, 
                v->allele1, v->allele2, v->freq1, v->freq2, v->cov); 
        }
        ++idx;
    }

    fclose(vars_outf);

}

void collapse_sites(hapstr& mask_global,
    int nvars,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_mask,
    map<int, int>& orig_to_collapsed,
    map<int, set<int> >& collapsed_to_orig,
    vector<pair<long int, int> >& clsort,
    bool keep_all_vars,
    double one){
    
    double zero = 1.0-one;

    // What sites appear to be erroneous and should be removed?
    // (NOTE: maybe dangerous to do this here)
    set<int> rm_err;

    for (int i = 0; i < clsort.size()-1; ++i){
        int idx1 = clsort[i].second;
        
        // If it's already been collapsed, stop looking
        if (orig_to_collapsed.count(idx1) > 0){
            continue;
        }
        else if (rm_err.find(idx1) != rm_err.end()){
            continue;
        }
        else{
            // Put self in self group
            orig_to_collapsed.insert(make_pair(idx1, idx1));
            if (collapsed_to_orig.count(idx1) == 0){
                set<int> s;
                collapsed_to_orig.insert(make_pair(idx1, s));
            }
            collapsed_to_orig[idx1].insert(idx1);
        }

        for (int j = i + 1; j < clsort.size(); ++j){
            int idx2 = clsort[j].second;
            
            // How many members of A are covered in B?
            int cladesize_A = 0;
            // How many members of B are covered in A?
            int cladesize_B = 0;
            // How many members are common to both clades?
            int cladesize_both = 0;
            // How many unique to A?
            int cladesize_AminusB = 0;
            // How many unique to B?
            int cladesize_BminusA = 0;

            for (robin_hood::unordered_set<unsigned long>::iterator a = 
                clades[idx1].begin(); a != clades[idx1].end(); ++a){
                if (clades_mask[idx2].find(*a) != clades_mask[idx2].end()){
                    cladesize_A++;
                    if (clades[idx2].find(*a) != clades[idx2].end()){
                        cladesize_both++;
                    }
                    else{
                        cladesize_AminusB++;
                    }
                }
            }
            for (robin_hood::unordered_set<unsigned long>::iterator b = 
                clades[idx2].begin(); b != clades[idx2].end(); ++b){
                if (clades_mask[idx1].find(*b) != clades_mask[idx1].end()){
                    cladesize_B++;
                    if (clades[idx1].find(*b) == clades[idx1].end()){
                        cladesize_BminusA++;
                    }
                }
            }
            
            // First, see if B looks like an error from A's perspective.
            double llAB_zero = dbinom(cladesize_A, cladesize_B, zero);
            if (llAB_zero > dbinom(cladesize_A, cladesize_B, 0.5) &&
                llAB_zero > dbinom(cladesize_A, cladesize_B, one)){
                if (!keep_all_vars){
                    // Remove site? (risky)
                    //mask_global.reset(idx2);
                    //rm_err.insert(idx2);
                }
            }
            else{
                double llAboth_one = dbinom(cladesize_A, cladesize_both, one);
                double llBboth_one = dbinom(cladesize_B, cladesize_both, one);
                
                // Two sites are same haplotype if A&B == A == B
                if (llAboth_one > dbinom(cladesize_A, cladesize_both, 0.5) &&
                    llAboth_one > dbinom(cladesize_A, cladesize_both, zero) &&
                    llBboth_one > dbinom(cladesize_B, cladesize_both, 0.5) &&
                    llBboth_one > dbinom(cladesize_B, cladesize_both, zero)){
                    
                    orig_to_collapsed.insert(make_pair(idx2, idx1));
                    collapsed_to_orig[idx1].insert(idx2);
                }
            }
        }
    }
    fprintf(stderr, "%ld collapsed site groups remaining\n", collapsed_to_orig.size());
}

/**
 * Given counts of alleles at variant sites per barcode, does several
 * things:
 *
 * Looks for sites that appear to be erroneous (in cells where the
 * minor allele is more common than the major allele, the sum of
 * log likelihoods that the allele is 50% frequent exceeds the
 * sum of log likelihoods that the allele is 100% frequent)
 *
 * Creates summarized, bitset representations of cell haplotypes
 * (instead of reads covering the minor and major allele,
 * stores 0 or 1 for presence or absence of the minor allele)
 *
 * Creates sets of cell barcodes (hashed as unsigned long) per
 * variant site, one for all cells with the minor allele,
 * one for all cells with the major allele, and one containing
 * both (all cells with coverage at that site)
 *
 * Collapses sites likely to tag the same haplotype into groups - 
 * this information is stored in collapsed_to_orig and 
 * orig_to_collapsed
 *
 * Combines the sets of cell barcodes per site (above) by group
 * of linked sites.
 */
void process_var_counts(robin_hood::unordered_map<unsigned long, var_counts>& hap_counter,
    robin_hood::unordered_map<unsigned long, hap>& bc2hap,
    bool keep_all_vars,
    bool skip_clustering,
    hapstr& mask_global,
    int nvars,
    bool barcodes_given,
    map<int, robin_hood::unordered_set<unsigned long> >& site_minor,
    map<int, robin_hood::unordered_set<unsigned long> >& site_major,
    map<int, robin_hood::unordered_set<unsigned long> >& site_mask,
    map<int, int>& orig_to_collapsed,
    map<int, set<int> >& collapsed_to_orig,
    vector<pair<long int, int> >& clsort,
    double one){

    double zero = 1.0-one;

    // Set up global mask -- include all sites by default
    mask_global.reset();
    for (int i = 0; i < nvars; ++i){
        mask_global.set(i);
    }
   
    // Remove sites that appear to be errors. This is determined by, in cells
    // where a site is more likely to be the minor than major allele, whether
    // the site is more likely to be heterozygous than homozygous for the minor
    // allele. 
    map<int, double> site_hetsum;
    map<int, double> site_homsum;
    
    for (int i = 0; i < nvars; ++i){
        if (mask_global.test(i)){
            site_hetsum.insert(make_pair(i, 0.0));
            site_homsum.insert(make_pair(i, 0.0));
            
            if (!skip_clustering){ 
                // need to cluster
                robin_hood::unordered_set<unsigned long> s;
                site_major.insert(make_pair(i, s));
                site_minor.insert(make_pair(i, s));
                site_mask.insert(make_pair(i, s));
            }
        }
    }
    
    // Create groups of cells with major/minor alleles at each site 
    // Simultaneously build cell haplotypes
    for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
        hap_counter.begin(); hc != hap_counter.end(); ++hc){
        hap h;
        for (int i = 0; i < nvars; ++i){
            if (mask_global.test(i)){
                int maj = hc->second.counts1[i];
                int min = hc->second.counts2[i];
                if (maj + min > 0){
                    double dmaj = dbinom(maj+min, min, zero);
                    double dhet = dbinom(maj+min, min, 0.5);
                    double dmin = dbinom(maj+min, min, one);
                    if (!keep_all_vars && dmin > dmaj){
                        site_homsum[i] += dmin;
                        site_hetsum[i] += dhet;
                    }
                    
                   if (dmaj > dmin && dmaj > dhet){
                        h.mask.set(i);
                        if (!skip_clustering){
                            site_major[i].insert(hc->first);
                            site_mask[i].insert(hc->first);
                        }
                   } 
                   else if (dmin > dmaj && dmin > dhet){
                        h.mask.set(i);
                        h.vars.set(i);
                        if (!skip_clustering){
                            site_minor[i].insert(hc->first);
                            site_mask[i].insert(hc->first);
                        }
                   }
                }
            }
        }
        bc2hap.emplace(hc->first, h);
    }
    
    if (!skip_clustering){      
        int sites_removed = 0;
        vector<pair<long int, int> > sitesort;
        for (map<int, double>::iterator s = site_homsum.begin(); s != site_homsum.end();
            ++s){
            if (!keep_all_vars && (
                site_homsum[s->first] == 0.0 || site_hetsum[s->first] >= 
                site_homsum[s->first])){
                // Remove
                mask_global.reset(s->first);
                site_major.erase(s->first);
                site_minor.erase(s->first);
                site_mask.erase(s->first);
                ++sites_removed;
            }
            else{
                sitesort.push_back(make_pair(-site_minor[s->first].size(), s->first));
            }
        }    
        if (!keep_all_vars){
            fprintf(stderr, "%d sites removed\n", sites_removed);
        }
        sort(sitesort.begin(), sitesort.end());

        // Now collapse sites & dump cell groups together.
        // Be more lenient with definition of 0 and 1 probability here -- for
        // identifying sites in cells, we want to be conservative, but for collapsing
        // sets of sites, we want to be more liberal.
        collapse_sites(mask_global, nvars, site_minor, site_major, site_mask, 
            orig_to_collapsed, collapsed_to_orig, sitesort, keep_all_vars, 0.999);
    
        // For each set of collapsed sites, dump cells together.
        for (map<int, int>::iterator oc = orig_to_collapsed.begin();
            oc != orig_to_collapsed.end(); ++oc){
            for (robin_hood::unordered_set<unsigned long>::iterator cell = 
                site_major[oc->first].begin(); cell != site_major[oc->first].end(); ++cell){
                /*
                if (site_mask[oc->second].find(*cell) == site_mask[oc->second].end()){
                    site_mask[oc->second].insert(*cell);
                    site_major[oc->second].insert(*cell);
                }
                */
                site_mask[oc->second].insert(*cell);
                if (site_minor[oc->second].find(*cell) == site_minor[oc->second].end()){
                    site_major[oc->second].insert(*cell);
                }
            }
            for (robin_hood::unordered_set<unsigned long>::iterator cell = 
                site_minor[oc->first].begin(); cell != site_minor[oc->first].end(); ++cell){
                /*
                if (site_mask[oc->second].find(*cell) == site_mask[oc->second].end()){
                    site_mask[oc->second].insert(*cell);
                    site_minor[oc->second].insert(*cell);
                }
                */
                site_mask[oc->second].insert(*cell);
                if (site_major[oc->second].find(*cell) == site_major[oc->second].end()){
                    site_minor[oc->second].insert(*cell);
                }
            }
        }
    
        // Finally, remove any uninformative sites (100% major or 100% minor allele)
        for (map<int, robin_hood::unordered_set<unsigned long> >::iterator sm = 
            site_minor.begin(); sm != site_minor.end(); ){
            
            if (sm->second.size() == 0 || site_major[sm->first].size() == 0){
                mask_global.reset(sm->first);
                site_major.erase(sm->first);
                site_mask.erase(sm->first);
                site_minor.erase(sm++);
            }
            else{
                clsort.push_back(make_pair(-sm->second.size(), sm->first));
                ++sm;
            }
        }
        sort(clsort.begin(), clsort.end());
    }
    fprintf(stderr, "%ld cell barcodes in data set\n", bc2hap.size());
}

/**
 * Assign barcodes of cells to a mitochondrial haplotype.
 */
void assign_bcs(robin_hood::unordered_map<unsigned long, var_counts>& hap_counter, 
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    vector<hap>& haps_final, 
    hapstr& mask_global,
    int nvars,
    double doublet_rate,
    bool use_filter,
    robin_hood::unordered_set<unsigned long>& cell_filter,
    double one){

    double zero = 1.0-one;
    
    vector<int> all_model_idx;
    for (int i = 0; i < haps_final.size(); ++i){
        if (doublet_rate < 1.0){
            all_model_idx.push_back(i);
        }
        if (doublet_rate > 0.0){
            for (int j = i + 1; j < haps_final.size(); ++j){
                int k = hap_comb_to_idx(i, j, haps_final.size());
                all_model_idx.push_back(k);
            }
        }
    }
    sort(all_model_idx.begin(), all_model_idx.end());
    
    int progress = 1000;
    int n_assigned = 0;

    for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
        hap_counter.begin(); hc != hap_counter.end(); ++hc){
        
        if (use_filter && cell_filter.find(hc->first) == cell_filter.end()){
            continue;
        }

        map<int, map<int, double> > llrs;
        for (int i = 0; i < all_model_idx.size()-1; ++i){
            int idx1 = all_model_idx[i];
            if (llrs.count(idx1) == 0){
                map<int, double> m;
                llrs.insert(make_pair(idx1, m));
            }
            for (int j = i + 1; j < all_model_idx.size(); ++j){
                int idx2 = all_model_idx[j];
                if (llrs[idx1].count(idx2) == 0){
                    // Start with prior
                    // Don't bother with prior if both singlet or
                    // both doublet -- since dealing only with
                    // likelihood ratios
                    double init_llr = 0.0;
                    if (idx1 < haps_final.size() && 
                        idx2 >= haps_final.size()){
                        init_llr = log2(1.0-doublet_rate) - 
                            log2(doublet_rate);
                    }
                    else if (idx2 < haps_final.size() &&
                        idx1 >= haps_final.size()){
                        init_llr = log2(doublet_rate) - 
                            log2(1.0-doublet_rate);
                    }
                    llrs[idx1].insert(make_pair(idx2, init_llr));
                }
            }
        }
        
        // Use 0.5 for every mixture proportion
        map<int, double> mixprops;
        for (int i = 0; i < all_model_idx.size(); ++i){
            int idx = all_model_idx[i];
            if (idx >= haps_final.size()){
                mixprops.insert(make_pair(idx, 0.5));
            }
        }

        int tot_reads = 0;
        int missing_sites = 0;
        bool skip_bc = false;
        
        for (int site = 0; site < nvars; ++site){
            if (mask_global[site]){
                
                // Retrieve major/minor allele counts
                int count1 = hc->second.counts1[site];
                int count2 = hc->second.counts2[site];
                
                if (count1+count2 > 0){                
                    // var set in hapstr == count2 is expected
                    for (int i = 0; i < all_model_idx.size()-1; ++i){
                        int idx1 = all_model_idx[i];
                        bool covered_idx1 = true;
                        // expected minor allele fraction this individual
                        double exp_frac1;
                        if (idx1 >= haps_final.size()){
                            pair<int, int> comb = idx_to_hap_comb(idx1, haps_final.size());
                            if (!haps_final[comb.first].mask[site] || 
                                !haps_final[comb.second].mask[site]){
                                covered_idx1 = false;
                            }
                            else{
                                int nmin1a = 0;
                                int nmin1b = 0;
                                if (haps_final[comb.first].vars[site]){
                                    nmin1a++;
                                }
                                if (haps_final[comb.second].vars[site]){
                                    nmin1b++;
                                }
                                double exp_frac1a = (double)nmin1a;
                                double exp_frac1b = (double)nmin1b;
                                exp_frac1 = mixprops[idx1]* exp_frac1a + (1.0-mixprops[idx1])*exp_frac1b;
                            }
                        }
                        else{
                            int nmin1 = 0;
                            if (!haps_final[idx1].mask[site]){
                                covered_idx1 = false;
                            }
                            else if (haps_final[idx1].vars[site]){
                                nmin1++;
                            }
                            exp_frac1 = (double)nmin1;
                        }
                        if (!covered_idx1){
                            continue;
                        }

                        for (int j = i + 1; j < all_model_idx.size(); ++j){
                            int idx2 = all_model_idx[j];
                            bool covered_idx2 = true;
                            double exp_frac2;
                            if (idx2 >= haps_final.size()){
                                pair<int, int> comb = idx_to_hap_comb(idx2, haps_final.size());
                                if (!haps_final[comb.first].mask[site] ||
                                    !haps_final[comb.second].mask[site]){
                                    covered_idx2 = false;
                                }
                                else{
                                    int nmin2a = 0;
                                    int nmin2b = 0;
                                    if (haps_final[comb.first].vars[site]){
                                        nmin2a++;
                                    }
                                    if (haps_final[comb.second].vars[site]){
                                        nmin2b++;
                                    }
                                    double exp_frac2a = (double)nmin2a;
                                    double exp_frac2b = (double)nmin2b;
                                    exp_frac2 = mixprops[idx2]*exp_frac2a + (1.0-mixprops[idx2])*exp_frac2b;
                                }
                            }
                            else{
                                int nmin2 = 0;
                                if (!haps_final[idx2].mask[site]){
                                    covered_idx2 = false;
                                }
                                else if (haps_final[idx2].vars[site]){
                                    nmin2++;
                                }
                                exp_frac2 = (double)nmin2;
                            }
                            
                            if (!covered_idx2){
                                continue;
                            }

                            if (exp_frac1 != exp_frac2){
                                
                                if (exp_frac1 == 0){
                                    exp_frac1 = zero;
                                }
                                else if (exp_frac1 == 1){
                                    exp_frac1 = one;
                                }
                                if (exp_frac2 == 0){
                                    exp_frac2 = zero;
                                }
                                else if (exp_frac2 == 1){
                                    exp_frac2 = one;
                                }
                                
                                double ll1 = dbinom(count1+count2, count2, exp_frac1);
                                double ll2 = dbinom(count1+count2, count2, exp_frac2);
                                
                                llrs[idx1][idx2] += (ll1-ll2);
                            }
                        }
                    }
                }
            }
        }
        
        double llr;
        int best_assignment = collapse_llrs(llrs, llr);
        if (best_assignment != -1 && llr > 0){
            assignments.emplace(hc->first, best_assignment);
            assignments_llr.emplace(hc->first, llr);
        }       
        ++n_assigned;
        if (!use_filter && n_assigned % progress == 0){
            fprintf(stderr, "%d cells assigned\r", n_assigned);
        }
    }
    if (!use_filter){
        fprintf(stderr, "%d cells assigned\n", n_assigned);
    }
}

/**
 * If a prior run created a set of variants, load them from that file here.
 * No filtering will be done on this set.
 */
string load_vars_from_file(string& filename, deque<varsite>& vars){
    ifstream f(filename);
    if (!f.good()){
        fprintf(stderr, "ERROR opening file %s\n", filename.c_str());
    }
    string line;
    string chrname = "";
    while (getline(f, line)){
        stringstream splitter(line);
        string field;
        string chrom;
        
        int pos;
        char allele1;
        char allele2;
        int fld_idx = 0;
        while (getline(splitter, field, '\t')){
            if (fld_idx == 0){
                chrom = field;    
                if (chrname != "" && chrname != chrom){
                    fprintf(stderr, "ERROR: chrname %s does not match %s\n", 
                    chrname.c_str(), chrom.c_str());
                    exit(1);
                }
                chrname = chrom;
            }
            else if (fld_idx == 1){
                // Convert from 1-based to 0-based
                pos = atoi(field.c_str())-1;
            }
            else if (fld_idx == 2){
                allele1 = field[0];
            }
            else if (fld_idx == 3){
                allele2 = field[0];
                
                // Check for duplicates
                bool add_var = true;
                if (vars.size() > 0){
                    if (vars.back().pos == pos){
                        fprintf(stderr, "WARNING: vars file contains duplicate \
entries for position %d. Using the first occurrence in the file.\n", pos);   
                        add_var = false;
                    }
                    else if (vars.back().pos > pos){
                        fprintf(stderr, "ERROR: non-sorted variant sites detected \
in vars file.\n");
                        exit(1);
                    }
                }
                if (add_var){
                    // Store data
                    varsite v;
                    v.pos = pos;
                    v.allele1 = allele1;
                    v.allele2 = allele2;
                    v.cov = 0;
                    v.freq1 = 0.0;
                    v.freq2 = 0.0;
                    //vars.insert(make_pair(pos, v));
                    vars.push_back(v);
                }
            }
            else{
                break;
            }
            fld_idx++;
        }
    }
    if (vars.size() > MAX_SITES){
        fprintf(stderr, "ERROR: loaded %ld variants from %s, but demux_mt was only \
compiled to allow up to %d variable sites.\n", vars.size(), filename.c_str(), MAX_SITES);
        fprintf(stderr, "Please re-compile using MAX_SITES=%ld or higher.\n", vars.size());
        exit(1);
    }
    return chrname;
}

/**
 * Find a (preliminary) set of variant sites on the mitochondrial sequence
 * using a BAM file. This will later be filtered further by 
 * allele frequency found across all cell barcodes.
 */
void find_vars_in_bam(string& bamfile, 
    string& mito_chrom, 
    int minmapq, 
    int minbaseq, 
    deque<varsite>& vars,
    bool has_bc_whitelist,
    set<unsigned long>& bc_whitelist,
    bool cov_filt){
    
    bam_mplp_t plp;
    
    infile_t infile(bamfile.c_str(), mito_chrom.c_str());
    
    infile_bc_wrapper wrap;
    wrap.infile = &infile;
    wrap.bcs = &bc_whitelist;

    if (has_bc_whitelist){
        // Important to use the correct pileup function, so we include the
        // cell barcode whitelist in the search (we only want sites that are
        // variable in the chosen set of cells)
        void* data[1];
        data[0] = (void*)&wrap;
        plp = bam_mplp_init(1, readaln_bcs, data);
    }
    else{
        void* data[1];
        data[0] = (void*)&infile;
        plp = bam_mplp_init(1, readaln, data);
    }
    if (!plp){
        fprintf(stderr, "ERROR opening BAM file %s as pileup\n", bamfile.c_str());
        exit(1);
    }
    bam_mplp_init_overlaps(plp);
    
    // Chromosome index - names in infile.fp_hdr->target_name array
    int tid = 0;
    // 0-based chromosome position
    int pos = 0;
    // Number of reads at position?
    int n = 0;
    
    // Store all alleles at each site (consider only SNPs)
    char alleles[5];
    int allelecounts[5];
    int n_alleles = 0;
    char alleles_filtered[5];
    int n_alleles_filtered = 0;
    int max_alleles = 4;
    const bam_pileup1_t* p;
    int ret;
    
    // Store initial set of variants, which will then be filtered for
    // coverage based on the median coverage across all found variants
    map<int, varsite> vars_unfiltered;
    
    vector<int> covsort;

    while ((ret = bam_mplp_auto(plp, &tid, &pos, &n, &p)) > 0){
        if (tid < 0){
            break;
        }
        if (tid >= infile.fp_hdr->n_targets){
            fprintf(stderr, "bam_mplp_auto returned tid %d greater than \
n_targets %d\n", tid, infile.fp_hdr->n_targets);
            exit(1);
        }

        // Process
        if (n > 0){
            n_alleles = 0;
            int cov = 0;
            int n_skip = 0;
            for (int i = 0; i < n; i++, p++){
                uint8_t* seq = bam_get_seq(p->b);
                uint8_t* qual = bam_get_qual(p->b);
                uint16_t flag = p->b->core.flag;
                int del_len, is_rev = bam_is_rev(p->b);
                if (!p->is_del && p->indel == 0 && !p->is_refskip){
                    // check map quality
                    if (p->b->core.qual >= minmapq){
                        unsigned char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
                        c = toupper(c);
                        unsigned char qualchr = qual[p->qpos] + 33;
                        if ((c == 'A' || c == 'C' || c == 'G' || c == 'T') && 
                            qualchr >= minbaseq){
                            cov++;
                            if (n_alleles == 0){
                                alleles[0] = c;
                                allelecounts[0] = 1;
                                n_alleles++;
                            }
                            else{
                                bool found = false;
                                for (int z = 0; z < n_alleles; ++z){
                                    if (alleles[z] == c){
                                        allelecounts[z]++;
                                        found = true;
                                        break;
                                    }
                                }
                                if (!found){
                                    alleles[n_alleles] = c;
                                    allelecounts[n_alleles] = 1;
                                    n_alleles++;
                                }
                            }
                            if (n_alleles >= max_alleles){
                                break; // no room left in array
                            }
                        }
                    }
                }
                else{
                    n_skip++;
                }
            }
            
            // If there are more than 2 alleles here, take the two highest-
            // frequency ones. If there's only one allele, it can't be a 
            // variable site.
            if (n_alleles >= 2){

                char allele1;
                char allele2;
                float freq1;
                float freq2;
                int num_pass = 0;
                
                vector<pair<double, char> > allelesort;
                for (int i = 0; i < n_alleles; ++i){
                    double freq = (double)allelecounts[i] / (double)cov;
                    allelesort.push_back(make_pair(-freq, alleles[i]));
                } 
                sort(allelesort.begin(), allelesort.end());
                
                allele1 = allelesort[0].second;
                freq1 = -allelesort[0].first;
                allele2 = allelesort[1].second;
                freq2 = -allelesort[1].first;

                varsite v;
                v.pos = pos;
                v.cov = cov;
                v.allele1 = allele1;
                v.allele2 = allele2;
                v.freq1 = freq1;
                v.freq2 = freq2;
                vars_unfiltered.insert(make_pair(pos, v));
                
                covsort.push_back(cov);
            }
        }
    }
    
    double cov_thresh = 0.0; 
    if (cov_filt){
        map<double, double> sitehist;
        sort(covsort.begin(), covsort.end()); 
        int cprev = -1;
        for (int i = 0; i < covsort.size(); ++i){
            if (covsort[i] != cprev){
                // All cov values greater than this one would survive
                // a filter set at this value.
                sitehist.insert(make_pair((double)covsort[i], (double)covsort.size()-i));
            }
            cprev = covsort[i];
        }
        cov_thresh = find_knee(sitehist, 0.25);
        fprintf(stderr, "Coverage threshold: %f\n", cov_thresh);
    }
    if (n < 0){
        fprintf(stderr, "bam_mplp_auto failed for %s\n", infile.fname);
        exit(1);
    }
    bam_mplp_destroy(plp);
    
    vector<pair<double, int> > vs_sort;
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != 
        vars_unfiltered.end(); ++v){
        if ((double)v->second.cov >= cov_thresh){
            vs_sort.push_back(make_pair(-(double)v->second.freq2, v->first));
        }
    }
    sort(vs_sort.begin(), vs_sort.end());
    set<int> pass_sites;
    for (int i = 0; i < MAX_SITES; ++i){
        if (i > vs_sort.size()-1){
            break;
        }
        pass_sites.insert(vs_sort[i].second);
    }
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != 
        vars_unfiltered.end(); ++v){
        if (pass_sites.find(v->first) != pass_sites.end()){
            // include
            vars.push_back(v->second);
        }
    }
}

/**
 * Given a set of variant sites, counts reads covering each allele of
 * each variant site tied to each barcode in the BAM file, across
 * the mitochondrial sequence.
 */
void count_vars_barcodes(string& bamfile, 
    string& mito_chrom, 
    int minmapq, 
    deque<varsite>& vars, 
    bool has_bc_whitelist, 
    set<unsigned long> & bc_whitelist, 
    robin_hood::unordered_map<unsigned long, var_counts>& hap_counter){
    
    int vars_idx = 0;
    
    bam_reader reader(bamfile);
    reader.set_cb();
    bool success = reader.set_query_region(mito_chrom.c_str(), -1, -1);
    if (!success){
        fprintf(stderr, "ERROR: sequence %s not found in BAM file.\n", mito_chrom.c_str());
        exit(1);
    }
    while (reader.next()){
        if (!reader.unmapped() && !reader.secondary() && reader.has_cb_z 
            && reader.mapq >= minmapq){

            // Get hashable version of barcode.
            string bc_str = reader.cb_z;
            size_t dashpos = bc_str.find('-');
            if (dashpos != string::npos){
                bc_str = bc_str.substr(0, dashpos);
            }
            bc bc_hashable;

            if (str2bc(bc_str.c_str(), bc_hashable) && (
                !has_bc_whitelist || bc_whitelist.find(bc_hashable.to_ulong()) != 
                bc_whitelist.end())){ 
                 
                // Remove any variants already done with
                while (vars.size() > 0 && vars.front().pos < reader.reference_start){
                    vars.pop_front();
                    vars_idx++;
                }

                if (vars.size() > 0 && vars.front().pos >= reader.reference_start && 
                    vars.front().pos <= reader.reference_end){
                    int vars_idx2 = 0;
                    for (deque<varsite>::iterator var = vars.begin(); var != vars.end(); 
                        ++var){
                        if (var->pos > reader.reference_end){
                            break;
                        }
                        else{
                            // Look for variant in read.
                            char base = reader.get_base_at(var->pos + 1);
                            
                            // Make sure an entry for this barcode exists.
                            if (hap_counter.count(bc_hashable.to_ulong()) == 0){
                                var_counts v;
                                hap_counter.emplace(bc_hashable.to_ulong(), v);
                            }
                            if (base == var->allele1){
                                hap_counter[bc_hashable.to_ulong()].counts1[vars_idx + vars_idx2]++;
                            }
                            else if (base == var->allele2){
                                hap_counter[bc_hashable.to_ulong()].counts2[vars_idx + vars_idx2]++;
                            }
                        }
                        ++vars_idx2;
                    }
                }
            } 
        } 
    }
}

/**
 * Returns chosen number of clusters and LLR sum of assignments,
 * disallowing doublet assignments.
 */
pair<int, float> infer_clusters(hapstr& mask_global, 
    robin_hood::unordered_map<unsigned long, hap>& haplotypes, 
    int nvars,
    vector<pair<long int, int> >& clades_sort,
    map<int, set<int> >& collapsed_to_orig,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    vector<hap>& haps_final,
    bool exact_matches_only,
    int nclust_max,
    robin_hood::unordered_set<unsigned long>& cellset,
    robin_hood::unordered_map<unsigned long, var_counts>& hap_counter,
    double one){
    
    hapstr mask;

    vector<hapstr> hapsites;
    vector<set<unsigned long> > hapgroups;
    vector<set<unsigned long> > hapgroups_not;
     
    vector<hapstr> hapsites_prev;
    vector<set<unsigned long> > hapgroups_prev;
    vector<set<unsigned long> > hapgroups_not_prev;

    int nsites_included = 0;
    vector<int> sites_order;
    vector<vector<hap> > haps_final_order;
    vector<hapstr> mask_order;
    vector<hap> haps_final_prev;
    
    int nclust_prev = -1;

    // Try including each site in decreasing order of confidence that
    // it is a true SNP (proxy for this: MAF)
    for (int i = 0; i < clades_sort.size(); ++i){
        
        int sitekey = clades_sort[i].second;
        ++nsites_included;
        
        // First site. Need to build new haplotypes.
        if (i == 0){

            // h1 == haplotype with minor allele at site
            hapstr h1;
            // h2 == haplotype with major allele at site
            hapstr h2;
            
            // Use all sites in the collapsed group.
            for (set<int>::iterator site = collapsed_to_orig[sitekey].begin();
                site != collapsed_to_orig[sitekey].end(); ++site){
                mask.set(*site);
                h1.set(*site);
            }
            
            // Store which (not-collapsed) sites are tied to each haplotype
            hapsites.push_back(h1);
            hapsites.push_back(h2);
            
            // Now create a clade of cells that possess the minor allele 
            // at this site 
            set<unsigned long> cl;
            for (robin_hood::unordered_set<unsigned long>::iterator member = 
                clades[sitekey].begin(); member != clades[sitekey].end(); ++member){
                cl.insert(*member);
            }
            // Create a clade of cells that possess the major allele
            // at this site
            set<unsigned long> cl_not;
            for (robin_hood::unordered_set<unsigned long>::iterator member = 
                clades_not[sitekey].begin(); member != clades_not[sitekey].end(); ++member){
                cl_not.insert(*member);
            }
            
            // Store all current haplotypes, to be modified with
            // the next site
            hapgroups.push_back(cl);
            hapgroups.push_back(cl_not);
            hapgroups_not.push_back(cl_not);
            hapgroups_not.push_back(cl);
                
        }
        else{
            
            // Not the first site: dig up preexisting haplotypes and 
            // modify them to each include the major or minor allele
            // at this site, thus increasing the number of possible
            // haplotypes. 
            
            // Add all sites in this collapsed group to the set
            // of included sites
            for (set<int>::iterator site = collapsed_to_orig[sitekey].begin();
                site != collapsed_to_orig[sitekey].end(); ++site){
                mask.set(*site);
            }

            // Build the new set of haplotypes that include
            // this site.         
            vector<hapstr> hapsites_new;
            vector<set<unsigned long> > hapgroups_new;
            vector<set<unsigned long> > hapgroups_not_new;
            
            for (int j = 0 ; j < hapsites.size(); ++j){
                
                // Each existing group must get two new versions - for this site
                // Existing haplotype with minor allele at this site 
                set<unsigned long> grp1;
                // Existing haplotype with major allele at this site
                set<unsigned long> grp2;
                
                for (robin_hood::unordered_set<unsigned long>::iterator member = 
                    clades[sitekey].begin(); member != clades[sitekey].end(); ++member){
                    if (hapgroups[j].find(*member) != hapgroups[j].end()){
                        grp1.insert(*member);
                    }
                }
                
                for (robin_hood::unordered_set<unsigned long>::iterator member = 
                    clades_not[sitekey].begin(); member != clades_not[sitekey].end(); ++member){
                    if (hapgroups[j].find(*member) != hapgroups[j].end()){
                        grp2.insert(*member);
                    }  
                }
                
                // In addition to groups of cells belonging to both new haplotypes
                // (grp1 and grp2 above), need to store sets of major/minor alleles
                // belonging to both new haplotypes.
                hapstr h1 = hapsites[j];
                hapstr h2 = hapsites[j];
                
                for (set<int>::iterator site = collapsed_to_orig[sitekey].begin();
                    site != collapsed_to_orig[sitekey].end(); ++site){
                    h1.set(*site);
                }
                
                if (grp1.size() > 0){
                    hapsites_new.push_back(h1);
                    hapgroups_new.push_back(grp1);
                    set<unsigned long> hnot = hapgroups_not[j];
                    for (robin_hood::unordered_set<unsigned long>::iterator member = 
                        clades_not[sitekey].begin(); member != clades_not[sitekey].end(); ++member){
                        hnot.insert(*member);    
                    }
                    hapgroups_not_new.push_back(hnot);
                }
                if (grp2.size() > 0){
                    hapsites_new.push_back(h2);
                    hapgroups_new.push_back(grp2);
                    set<unsigned long> hnot = hapgroups_not[j];
                    for (robin_hood::unordered_set<unsigned long>::iterator member = 
                        clades[sitekey].begin(); member != clades[sitekey].end(); ++member){
                        hnot.insert(*member);    
                    }
                    hapgroups_not_new.push_back(hnot);
                }
            }
            
            // Store previous versions of all haplotypes and their members,
            // because the algorithm only tells us when we've gone too far
            // and added too many sites (in that case, we will want to select
            // the set of haplotypes we had before including the last site).
                
            hapsites_prev.clear();
            hapsites_prev = hapsites;
            hapgroups_prev.clear();
            hapgroups_prev = hapgroups;
            hapgroups_not_prev.clear();
            hapgroups_not_prev = hapgroups_not;
            
            // Set currently-tracked haplotypes, and their sites and members,
            // to the new versions updated to include the new site.
             
            hapsites.clear();
            hapsites = hapsites_new;
            hapgroups.clear();
            hapgroups = hapgroups_new;
            hapgroups_not.clear();
            hapgroups_not = hapgroups_not_new;
    
        }

        // Now determine likeliest number of haplotypes, given the current
        // set of sites & clades.
        
        // Compute the size of each currently-tracked haplotype    
        vector<int> sizes;
        // Store sizes of haplotypes mapped to haplotype indices
        vector<pair<int, int> > sizepairs;
        
        for (int i = 0; i < hapgroups.size(); ++i){
            if (hapgroups[i].size() > 0){
                sizes.push_back((int)hapgroups[i].size());
                sizepairs.push_back(make_pair((int)hapgroups[i].size(), i));
            }
        }
        
        int nclust = choose_nclust(sizes);
        if (nclust > nclust_prev){
            fprintf(stderr, "Found %d haplotypes with %d site groups included\n", 
                nclust, nsites_included);
        }

        // Should we stop adding sites here?
    
        // At first, we just want to get the set of sites that results in the 
        // maximum number of clusters (we will deal with the possibility that
        // we included too many + erroneous sites later).
        bool terminate = false;
        
        if (nclust_max != -1 && nclust > nclust_max){
            // Terminate, if the user has specified a set number of clusters
            // and we're about to exceed it.
            terminate = true;
            // Avoid reporting clusters beyond the user-specified number
            nclust = nclust_max;             
        }
        else if (nclust < nclust_prev){
            // Terminate, if the number of clusters dropped because we included
            // too many sites.
            terminate = true;
        }
        
        if (terminate){
            
            // Don't include the last site.
            mask.reset(sitekey);
            
            // Prepare to return everything from the previous iteration
            hapgroups.clear();
            hapgroups = hapgroups_prev;
            hapgroups_not.clear();
            hapgroups_not = hapgroups_not_prev;
            hapsites.clear();
            hapsites = hapsites_prev;
        
            break;
        }
        else{
            
            nclust_prev = nclust;

            haps_final.clear();
            sort(sizepairs.begin(), sizepairs.end());
            
            for (int i = sizepairs.size()-1; i >= sizepairs.size()-nclust; --i){
                hap h;
                h.mask = mask;
                h.vars = hapsites[sizepairs[i].second];
                
                // Remove sites missing in the majority of cells?
                vector<int> site_maj;
                vector<int> site_min;
                vector<int> site_miss;
                for (int x = 0; x < nvars; ++x){
                    if (mask_global[x] && mask[x]){
                        site_maj.push_back(0);
                        site_min.push_back(0);
                        site_miss.push_back(0);
                    }
                }
                for (set<unsigned long>::iterator cell = 
                    hapgroups[sizepairs[i].second].begin();
                    cell != hapgroups[sizepairs[i].second].end(); ++cell){
                    int site_idx = 0;
                    for (int x = 0; x < nvars; ++x){
                        if (mask_global[x] && mask[x]){
                            if (haplotypes[*cell].mask[x]){
                                if (haplotypes[*cell].vars[x]){
                                    site_min[site_idx]++;
                                }
                                else{
                                    site_maj[site_idx]++;
                                }
                            }
                            else{
                                site_miss[site_idx]++;
                            }
                            ++site_idx;
                        }
                    }
                }
                int site_idx = 0;
                for (int x = 0; x < nvars; ++x){
                    if (mask_global[x] && mask[x]){
                        int maj = site_maj[site_idx];
                        int min = site_min[site_idx];
                        int miss = site_miss[site_idx];
                        double llmaj = dbinom(maj+min, maj, one);
                        double llmin = dbinom(maj+min, min, one);
                        double llmiss = dbinom(maj+min+miss, miss, one);
                        if (llmiss > llmaj && llmiss > llmin){
                            h.mask.reset(x);
                        }
                        ++site_idx;
                    }
                }

                haps_final.push_back(h);
                //haps_final.push_back(hapsites[sizepairs[i].second]);
                /*
                fprintf(stderr, "\tsize %d [ ", sizepairs[i].first);
                for (int x = 0; x < nvars; ++x){
                    if (mask[x] && hapsites[sizepairs[i].second][x]){
                        fprintf(stderr, "%d ", x);
                    }
                }
                fprintf(stderr, "]\n");
                */
            }
            
            // Store haps_final and mask, if haps_final is different from the 
            // previous version
            bool hf_changed = false;
            if (haps_final.size() != haps_final_prev.size()){
                hf_changed = true;
            }         
            else{
                for (int i = 0; i < haps_final.size(); ++i){
                    if (haps_final_prev[i].mask != haps_final[i].mask ||
                        haps_final_prev[i].vars != haps_final[i].vars){
                        hf_changed = true;
                        break;
                    }
                }
            }
            if (hf_changed){
                sites_order.push_back(sitekey);        
                haps_final_order.push_back(haps_final);
                mask_order.push_back(mask);
            }
            
            if (!exact_matches_only){
                // Add in newly allowable cells to each haplotype
                for (robin_hood::unordered_map<unsigned long, hap>::iterator hap = 
                    haplotypes.begin(); hap != haplotypes.end(); ++hap){
                    if ((mask & hap->second.mask).count() < mask.count()){
                        int match_idx = -1;
                        int match_count = 0;
                        for (int i = 0; i < hapsites.size(); ++i){
                            if (hapgroups_not[i].find(hap->first) == hapgroups_not[i].end()){
                                if ((mask & hap->second.mask & hap->second.vars) == 
                                    (mask & hap->second.mask & hapsites[i])){
                                    match_idx = i;
                                    match_count++;
                                }
                            }
                        }
                        if (match_count == 1){
                            hapgroups[match_idx].insert(hap->first);
                        }
                    }
                }
            }
        }
    }
    
    fprintf(stderr, "Iteratively removing site groups to test for improved\n");
    fprintf(stderr, "   assignment LLRs...\n"); 

    // One last step: we chose the set of included sites to maximize
    // the number of haplotypes. This guarantees we included enough
    // sites to see the true number of haplotypes, but it also 
    // risks having included too many sites and are seeing small numbers
    // of erroneous haplotypes. 

    // To combat this, go through all valid cells (defined as having coverage
    // at the maximum-frequency site) and choose the best fit haplotype
    // based on actual read counts. Erroneous sites will have low read counts
    // and not affect the totals much.

    // Then reinfer the number of clusters from this, and finally, remove
    // from all haplotypes any sites that are not informative for distinguishing
    // between chosen clusters.
    
    robin_hood::unordered_map<unsigned long, int> assn;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    double llrsum = 0.0;
    assign_bcs(hap_counter, assn, assn_llr, haps_final, mask, nvars, 0.0, true, 
        cellset, one);

    //map<int, int> grpsizes;
    for (robin_hood::unordered_map<unsigned long, double>::iterator al = assn_llr.begin();
        al != assn_llr.end(); ++al){
        llrsum += al->second;
        //if (grpsizes.count(assn[al->first]) == 0){
        //    grpsizes.insert(make_pair(assn[al->first], 0));
        //}
        //grpsizes[assn[al->first]]++;
    }
    /*
    vector<int> sizevec;
    for (map<int, int>::iterator gs = grpsizes.begin(); gs != grpsizes.end(); ++gs){
        sizevec.push_back(gs->second);
    }
    int nclust_from_assn = choose_nclust(sizevec);
    fprintf(stderr, "size from assn: %d\n", nclust_from_assn);
    */
    double llrsum_orig = llrsum;

    int n_rm = 0; 
    // Try iteratively deleting sites backward until we maximize llrsum
    if (haps_final_order.size() > 2){
        for (int site_idx = haps_final_order.size()-2; site_idx >= 0; --site_idx){
            if (haps_final_order[site_idx].size() == 1){
                // Don't consider single haplotype case.
                haps_final = haps_final_order[site_idx+1];
                mask = mask_order[site_idx+1];
                break;
            }
            else{
                ++n_rm;
                assn.clear();
                assn_llr.clear();
                double llrsum_new = 0;
                assign_bcs(hap_counter, assn, assn_llr, haps_final_order[site_idx],
                    mask_order[site_idx], nvars, 0.0, true, cellset, one);
                //grpsizes.clear();
                //sizevec.clear();
                for (robin_hood::unordered_map<unsigned long, double>::iterator al = 
                    assn_llr.begin(); al != assn_llr.end(); ++al){
                    llrsum_new += al->second;
                    //if (grpsizes.count(assn[al->first]) == 0){
                    //    grpsizes.insert(make_pair(assn[al->first], 0));
                    //}
                    //grpsizes[assn[al->first]]++;
                }
                /*
                for (map<int, int>::iterator gs = grpsizes.begin(); gs != grpsizes.end();
                    ++gs){
                    sizevec.push_back(gs->second);
                }
                nclust_from_assn = choose_nclust(sizevec);
                fprintf(stderr, "nclust from assn: %d\n", nclust_from_assn);
                */
                if (llrsum_new < llrsum){
                    // Finished
                    haps_final = haps_final_order[site_idx+1];
                    mask = mask_order[site_idx+1];
                    break;
                }
                else{
                    if (haps_final_order[site_idx].size() != 
                        haps_final_order[site_idx+1].size()){
                        fprintf(stderr, "%ld haplotypes after removing %d site groups\n", 
                            haps_final_order[site_idx].size(), n_rm);
                        fprintf(stderr, "\timprovement: %f\n", llrsum_new-llrsum_orig);
                    }
                    llrsum = llrsum_new;
                }
            }
        } 
    }

    fprintf(stderr, "Final clustering: %ld haplotypes\n", haps_final.size());
    int nhaps = (int)haps_final.size();
    double llr = llrsum;
    
    // Make sure we don't look at sites that weren't included in final haplotypes.
    mask_global &= mask;

    return make_pair(nhaps, llrsum);
}


void parse_idsfile(string& idsfilename, vector<string>& ids){
    ifstream infile(idsfilename.c_str());
    if (!infile.good()){
        fprintf(stderr, "ERROR opening file %s\n", idsfilename.c_str());
        exit(1);
    }
    string id;
    while (infile >> id ){
        ids.push_back(id);
    }
}

void write_bchaps(string& haps_out,
    int nvars,
    hapstr& mask_global,
    robin_hood::unordered_map<unsigned long, hap>& haplotypes,
    string& bc_group,
    bool cellranger,
    bool seurat,
    bool underscore){
    
    // Open output file
    FILE* haps_outf = fopen(haps_out.c_str(), "w");
    
    // Write header
    fprintf(haps_outf, "bc\t");
    bool firstprint = true;
    int varidx = 0;
    for (int i = 0; i < nvars; ++i){
        if (mask_global.test(i)){
            if (!firstprint){
                fprintf(haps_outf, "\t");
            }
            fprintf(haps_outf, "%d", i);
            ++varidx;
            firstprint = false;
        }
    }
    fprintf(haps_outf, "\n");
    
    // Visit and print information for each cell barcode
    for (robin_hood::unordered_map<unsigned long, hap>::iterator h = 
        haplotypes.begin(); h != haplotypes.end(); ++h){
        
        bc bcbits(h->first);
        string bcstr = bc2str(bcbits);
        mod_bc_libname(bcstr, bc_group, cellranger, seurat, underscore);

        if ((h->second.mask & mask_global).count() > 0){
            firstprint = true;
            fprintf(haps_outf, "%s\t", bcstr.c_str());
        
            for (int x = 0; x < nvars; ++x){
                if (mask_global[x]){
                    if (!firstprint){
                        fprintf(haps_outf, "\t");
                    }
                    firstprint = false;
                    if (h->second.mask[x]){
                        if (h->second.vars[x]){
                            fprintf(haps_outf, "1");
                        }
                        else{
                            fprintf(haps_outf, "0");
                        }
                    }
                    else{
                        fprintf(haps_outf, "NA");
                    }
                }
            }
            fprintf(haps_outf, "\n");
        }
    }
    fclose(haps_outf);
}

void write_clusthaps(string& clusthapsfile,
    hapstr& mask_global,
    vector<hap>& clusthaps){
    
    FILE* clusthaps_f = fopen(clusthapsfile.c_str(), "w");
    for (int i = 0; i < clusthaps.size(); ++i){
        for (int x = 0; x < nvars; ++x){
            if (mask_global[x]){
                if (clusthaps[i].mask.test(x)){
                    if (clusthaps[i].vars.test(x)){
                        fprintf(clusthaps_f, "1");
                    }
                    else{
                        fprintf(clusthaps_f, "0");
                    }
                }
                else{
                    fprintf(clusthaps_f, "-");
                }
            }
        }
        fprintf(clusthaps_f, "\n");
    }
    fclose(clusthaps_f);
}

void load_clusthaps(string& hapsfilename,
    hapstr& mask_global,
    vector<hap>& clusthaps){

    ifstream infile(hapsfilename.c_str());
    if (!infile.good()){
        fprintf(stderr, "ERROR opening file %s\n", hapsfilename.c_str());
        exit(1);
    }
    string hapstring;
    bool first = true;
    while (infile >> hapstring ){
        if (first){
            first = false;
            nvars = hapstring.length();
            for (int i = 0; i < nvars; ++i){
                mask_global.set(i);
            }
        }
        hap h;
        for (int i = 0; i < nvars; ++i){
            if (hapstring[i] == '0'){
                h.mask.set(i);
            }
            else if (hapstring[i] == '1'){
                h.mask.set(i);
                h.vars.set(i);
            }
        }
        clusthaps.push_back(h);
    }
    /*
    // Check for any positions where there are no differences among haplotypes
    // and blacklist them.
    for (int i = 0; i < nvars; ++i){
        if (mask_global[i]){
            int nmaj = 0;
            int nmin = 0;
            for (int x = 0; x < clusthaps.size(); ++x){
                if (clusthaps[x].mask[i]){
                    if (clusthaps[x].vars[i]){
                        nmin++;
                    }
                    else{
                        nmaj++;
                    }
                }
            }
            if (nmaj == 0 || nmin == 0){
                //fprintf(stderr, "remove site %d\n", i);
                // Uninformative site.
                mask_global.reset(i);
            }
        }
    }
    */
}

/**
 * Print barcode-identity assignments to disk.
 */
void write_assignments(string& assn_out,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore,
    vector<string>& clust_ids,
    int nhaps,
    map<int, int>& id_counter,
    int& tot_cells,
    int& doub_cells){
    
    tot_cells = 0;
    doub_cells = 0;
    
    FILE* assn_outf = fopen(assn_out.c_str(), "w");
    char namebuf[100];
    for (robin_hood::unordered_map<unsigned long, int>::iterator assn = 
        assignments.begin(); assn != assignments.end(); ++assn){
        tot_cells++;
        string name;
        char s_d;
        if (clust_ids.size() > 0){
            name = idx2name(assn->second, clust_ids);
            if (assn->second >= nhaps){
                s_d = 'D';
                doub_cells++;
            }
            else{
                s_d = 'S';
            }
        }
        else{
            if (assn->second < nhaps){
                sprintf(&namebuf[0], "%d", assn->second);
                s_d = 'S';
            }
            else{
                pair<int, int> combo = idx_to_hap_comb(assn->second, 
                    nhaps);
                sprintf(&namebuf[0], "%d+%d", combo.first, combo.second);
                s_d = 'D';
                doub_cells++;
            }
            name = namebuf;
        }
        if (id_counter.count(assn->second) == 0){
            id_counter.insert(make_pair(assn->second, 0));
        }
        id_counter[assn->second]++;
        bc as_bitset(assn->first);
        string bc_str = bc2str(as_bitset);
        mod_bc_libname(bc_str, barcode_group, cellranger, seurat, underscore);
        fprintf(assn_outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(),
            name.c_str(), s_d, assignments_llr[assn->first]);
    } 
    fclose(assn_outf);
}

void write_statsfile(string& statsfilename,
    string& output_prefix,
    int nclust_model,
    double llrsum_model,
    int tot_cells,
    int doub_cells,
    map<int, int>& id_counter,
    double doublet_rate,
    string& hapsfilename,
    string& varsfilename,
    vector<string>& samples){

    FILE* statsfilef = fopen(statsfilename.c_str(), "w");
    
    fprintf(statsfilef, "%s\tparam.doublet_rate\t%f\n", output_prefix.c_str(),
        doublet_rate);

    if (nclust_model != -1){
        // Clusters were inferred.
        fprintf(statsfilef, "%s\tnum_clusters\t%d\n", output_prefix.c_str(), 
            nclust_model);
        fprintf(statsfilef, "%s\tllr_sum\t%f\n", output_prefix.c_str(), llrsum_model);
    }
    else{
        // Clusters were loaded
        fprintf(statsfilef, "%s\tparam.hapsfile\t%s\n", output_prefix.c_str(),
            hapsfilename.c_str());
        fprintf(statsfilef, "%s\tparam.varsfile\t%s\n", output_prefix.c_str(),
            varsfilename.c_str());
    }
    fprintf(statsfilef, "%s\ttot_cells\t%d\n", output_prefix.c_str(), tot_cells);
    if (doublet_rate > 0 && doublet_rate < 1){
        // Prevent NaN from division by zero
        if (tot_cells == 0){
            tot_cells = 1;
        }
        fprintf(statsfilef, "%s\tdoublets\t%d\n", output_prefix.c_str(), 
            doub_cells);
        fprintf(statsfilef, "%s\tfrac_doublets\t%.3f\n", output_prefix.c_str(), 
            (float)doub_cells/(float)tot_cells);
        
        // Compute chi-squared p value for expected vs actual counts of each doublet type
        double chisq_p = doublet_chisq(id_counter, samples.size());
        if (chisq_p >= 0){
            fprintf(statsfilef, "%s\tdoublet_chisq.p\t%f\n", output_prefix.c_str(),
                chisq_p);
        }
    }
     
    vector<pair<int, int> > idcounts_sorted;
    for (map<int, int>::iterator idc = id_counter.begin(); idc != 
        id_counter.end(); ++idc){
        idcounts_sorted.push_back(make_pair(-idc->second, idc->first));
    }
    sort(idcounts_sorted.begin(), idcounts_sorted.end());
    for (vector<pair<int, int> >::iterator idcs = idcounts_sorted.begin(); 
        idcs != idcounts_sorted.end(); ++idcs){
        fprintf(statsfilef, "%s\t%s\t%d\n", output_prefix.c_str(), 
            idx2name(idcs->second, samples).c_str(), -idcs->first);   
    }
    fclose(statsfilef);
}

/**
 * If inferring mixing proportions, need to load assignments of cells to mitochondrial
 * IDs.
 *
 */
void parse_assignments(string& assnfile,
    vector<string>& hapnames,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    bool doublets_only){
    
    map<string, int> id2idx;
    for (int i = 0; i < hapnames.size(); ++i){
        id2idx.insert(make_pair(hapnames[i], i));
    }

    ifstream infile(assnfile.c_str());
    if (!infile.good()){
        fprintf(stderr, "ERROR opening file %s\n", assnfile.c_str());
        exit(1);
    }
    string bc_str;
    string id_str;
    char s_d;
    double llr;
    
    int doubcount = 0;

    while (infile >> bc_str >> id_str >> s_d >> llr){
        // Process barcode string
        size_t sep_pos = bc_str.find("-");
        if (sep_pos != string::npos){
            bc_str = bc_str.substr(0, sep_pos);
        }
        unsigned long as_ul = bc_ul(bc_str);
        if (s_d == 'D'){
            size_t sep_pos = id_str.find("+");
            if (sep_pos == string::npos){
                sep_pos = id_str.find("x");
                if (sep_pos == string::npos){
                    sep_pos = id_str.find("X");
                }
            }
            if (sep_pos == string::npos){
                fprintf(stderr, "ERROR: could not process ID string %s as a doublet\n", 
                    id_str.c_str());
                exit(1);
            }
            else{
                string id1 = id_str.substr(0, sep_pos);
                string id2 = id_str.substr(sep_pos+1, id_str.length()-sep_pos-1);
                if (id2idx.count(id1) == 0){
                    fprintf(stderr, "ERROR: id %s not found in haplotypes. Did you forget \
to pass the .ids file (-i)?\n", id1.c_str());
                    exit(1);
                }
                if (id2idx.count(id2) == 0){
                    fprintf(stderr, "ERROR: id %s not found in haplotypes. Did you forget \
to pass the .ids file (-i)?\n", id2.c_str());
                    exit(1);
                }
                int idx1 = id2idx[id1];
                int idx2 = id2idx[id2];
                
                int idx3;
                bool is_doub = true;
                if (idx1 < idx2){
                    idx3 = hap_comb_to_idx(idx1, idx2, hapnames.size());
                }
                else if (idx2 < idx1){
                    idx3 = hap_comb_to_idx(idx2, idx1, hapnames.size());
                }
                else{
                    // They're the same.
                    idx3 = idx1;
                    is_doub = false;
                }
                if (!doublets_only || is_doub){
                    assignments.emplace(as_ul, idx3);
                }
            }
        }
        else if (s_d == 'S' && !doublets_only){
            if (id2idx.count(id_str) == 0){
                fprintf(stderr, "ERROR: id %s not found in haplotypes. Did you forget to \
pass the .ids file (-i)?\n", id_str.c_str());
                exit(1);
            }
            int idx = id2idx[id_str];
            assignments.emplace(as_ul, idx);
        }
    }
}

double y_mixprop(double p, const map<string, double>& data_d, const map<string, int>& data_i){
    int m1 = data_i.at("match1");
    int m2 = data_i.at("match2");
    double k = (double)m1;
    double n = (double)(m1+m2);
    return logbinom(n, k, p);
}

double dy_dp_mixprop(double p, const map<string, double>& data_d, const map<string, int>& data_i){
    int m1 = data_i.at("match1");
    int m2 = data_i.at("match2");
    double k = (double)m1;
    double n = (double)(m1+m2);
    return (k-n*p)/(p - p*p);
}

double d2y_dp2_mixprop(double p, const map<string, double>& data_d, const map<string, int>& data_i){
    int m1 = data_i.at("match1");
    int m2 = data_i.at("match2");
    double k = (double)m1;
    double n = (double)(m1+m2);
    return (k*(2*p-1) - n*p*p)/(pow(p-1,2)*p*p);
}

void infer_mixprops(robin_hood::unordered_map<unsigned long, var_counts>& hap_counter,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    vector<hap>& clusthaps,
    hapstr& mask_global,
    int nvars,
    robin_hood::unordered_map<unsigned long, double>& mixprops_mean,
    robin_hood::unordered_map<unsigned long, double>& mixprops_sd,
    map<int, double>& id_mixprop_mean,
    map<int, double>& id_mixprop_sd,
    vector<string>& clust_ids){

    int n_samples = (int)clusthaps.size();
    
    map<int, vector<int> > id_match1;
    map<int, vector<int> > id_match2;

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assignments.begin();
        a != assignments.end(); ++a){
        
        if (hap_counter.count(a->first) > 0 && a->second >= n_samples){
            
            // Doublet.
            pair<int, int> combo = idx_to_hap_comb(a->second, n_samples);
            
            int match1 = 0;
            int match2 = 0;

            for (int x = 0; x < nvars; ++x){
                
                if (mask_global[x] && clusthaps[combo.first].mask[x] &&
                    clusthaps[combo.second].mask[x] && 
                    clusthaps[combo.first].vars[x] != clusthaps[combo.second].vars[x]){
                    
                    // We can use this site.
                    int nmaj = hap_counter[a->first].counts1[x];
                    int nmin = hap_counter[a->first].counts2[x];
                    if (nmaj + nmin > 0){
                        if (clusthaps[combo.first].vars[x] && !clusthaps[combo.second].vars[x]){
                            // minor allele is indv 1
                            match1 += nmin;
                            match2 += nmaj;
                        }
                        else if (clusthaps[combo.second].vars[x] && !clusthaps[combo.first].vars[x]){
                            // minor allele is indv 2
                            match1 += nmaj;
                            match2 += nmin;
                        }
                    }
                }
            }
            
            if (match1 + match2 > 0){
                string bcstr = bc2str(a->first);
                fprintf(stdout, "%s\t%s\t%s\t%d\t%d\n", bcstr.c_str(), clust_ids[combo.first].c_str(),
                    clust_ids[combo.second].c_str(), match1, match2);

                /*        
                double mean = (double)match1/(double)(match1+match2);
                double var = (double)(match1*match2)/(pow((double)match1+match2,2)*(double)(match1+match2+1));
                string n = idx2name(a->second, clust_ids);
                pair<int, int> combo = idx_to_hap_comb(a->second, clust_ids.size());
                string n1 = idx2name(combo.first, clust_ids);
                string n2 = idx2name(combo.second, clust_ids);
                fprintf(stdout, "%s\t%s\t%f\t%f\n", n1.c_str(), n2.c_str(), mean, sqrt(var));
                fprintf(stdout, "%s\t%s\t%f\t%f\n", n2.c_str(), n1.c_str(), mean, sqrt(var));
                continue;

                vector<int> match1v = {match1};
                vector<int> match2v = {match2};
                optimML::brent_solver solver(y_mixprop, dy_dp_mixprop, d2y_dp2_mixprop);
                solver.add_data("match1", match1v);
                solver.add_data("match2", match2v);
                solver.constrain_01();
                double p = solver.solve(0,1);
                mixprops_mean.emplace(a->first, p);
                mixprops_sd.emplace(a->first, solver.se); 
                */
                // Store for whole-ID tests
                if (id_match1.count(a->second) == 0){
                    vector<int> v;
                    id_match1.insert(make_pair(a->second, v));
                }
                id_match1[a->second].push_back(match1);
                if (id_match2.count(a->second) == 0){
                    vector<int> v;
                    id_match2.insert(make_pair(a->second, v));
                }
                id_match2[a->second].push_back(match2);
            }
        }
    }

    for (map<int, vector<int> >::iterator im = id_match1.begin(); im != id_match1.end();
        ++im){
         
        optimML::brent_solver solver(y_mixprop, dy_dp_mixprop, d2y_dp2_mixprop);
        solver.add_data("match1", im->second);
        solver.add_data("match2", id_match2[im->first]);
        solver.constrain_01();
        //solver.add_beta_prior(10000, 10000);
        double p = solver.solve(0,1);
	    pair<int, int> comb = idx_to_hap_comb(im->first, clust_ids.size());
        //fprintf(stdout, "%s\t%s\t%f\t%f\n", clust_ids[comb.first].c_str(), clust_ids[comb.second].c_str(), p, solver.se);
        id_mixprop_mean.insert(make_pair(im->first, p));
        id_mixprop_sd.insert(make_pair(im->first, solver.se));
        
    }
}

void write_mixprops(FILE* outf,
    string& barcode_group,
    bool cellranger,
    bool seurat,
    bool underscore,
    robin_hood::unordered_map<unsigned long, double>& mixprops_mean,
    robin_hood::unordered_map<unsigned long, double>& mixprops_sd,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    vector<string>& clust_ids){
    
    for (robin_hood::unordered_map<unsigned long, double>::iterator mp = mixprops_mean.begin();
        mp != mixprops_mean.end(); ++mp){
        
        string bc_str = bc2str(mp->first);
        mod_bc_libname(bc_str, barcode_group, cellranger, seurat, underscore);

        int assn = assignments[mp->first];
        pair<int, int> combo = idx_to_hap_comb(assn, clust_ids.size());
        
        string name1 = idx2name(combo.first, clust_ids);
        string name2 = idx2name(combo.second, clust_ids);
        
        fprintf(outf, "%s\t%s\t%s\t%f\t%f\n", bc_str.c_str(), name1.c_str(), name2.c_str(),
            mp->second, mixprops_sd[mp->first]);
        fprintf(outf, "%s\t%s\t%s\t%f\t%f\n", bc_str.c_str(), name2.c_str(), name1.c_str(),
            mixprops_sd[mp->first], mp->second);

    }

}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"barcodes_filter", required_argument, 0, 'f'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"mapq", required_argument, 0, 'q'},
       {"baseq", required_argument, 0, 'Q'},
       {"nclust_max", required_argument, 0, 'N'},
       {"mito", required_argument, 0, 'm'},
       {"dump", no_argument, 0, 'd'},
       {"vars", required_argument, 0, 'v'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"haps", required_argument, 0, 'H'},
       {"ids", required_argument, 0, 'i'},
       {"no_cov_filt", no_argument, 0, 'c'},
       {"assignments", required_argument, 0, 'a'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string output_prefix;
    int minmapq = 20;
    int minbaseq = 20;
    int nclust = -1;
    bool dump = false;
    string mito_chrom = "chrM";
    string varsfile;
    bool varsfile_given = false;
    string barcodesfilename;
    string barcodes_filter_filename;
    string hapsfilename;
    bool hapsfile_given = false;
    double doublet_rate = 0.5;
    bool ids_given = false;
    string idsfile;
    string barcode_group = "";
    bool mixing_proportions = false;
    string assnfile = "";
    bool cov_filt = true;
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:n:B:f:g:q:Q:N:m:v:H:i:D:a:CSUcdh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'n':
                barcode_group = optarg;
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
            case 'H':
                hapsfilename = optarg;
                hapsfile_given = true;
                break;
            case 'i':
                idsfile = optarg;
                ids_given = true;
                break;
            case 'b':
                bamfile = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'B':
                barcodesfilename = optarg;
                break;
            case 'f':
                barcodes_filter_filename = optarg;
                break;
            case 'q':
                minmapq = atoi(optarg);
                break;
            case 'Q':
                minbaseq = atoi(optarg);
                break;
            case 'c':
                cov_filt = false;
                break;
            case 'm':
                mito_chrom = optarg;
                break;
            case 'v':
                varsfile = optarg;
                varsfile_given = true;
                break;
            case 'd':
                dump = true;
                break;
            case 'a':
                mixing_proportions = true;
                assnfile = optarg;
                break;
            case 'N':
                nclust = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }

    // Error check arguments.
    if (bamfile.length() == 0){
        fprintf(stderr, "ERROR: bam file is required\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix is required\n");
        exit(1);
    }
    if (hapsfile_given && !varsfile_given){
        fprintf(stderr, "ERROR: If --haps / -H is given, --vars / -v is also required\n");
        exit(1);
    }
    if (doublet_rate < 0 || doublet_rate > 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1, inclusive\n");
        exit(1);
    }
    if (mixing_proportions && (!hapsfile_given || !varsfile_given)){
        fprintf(stderr, "ERROR: if inferring mixing proportions, -H, -v, and -a are all required\n");
        exit(1);
    }
    if (mixing_proportions && dump){
        fprintf(stderr, "ERROR: cannot infer mixing proportions and also dump counts\n");
        exit(1);
    }
    if (nclust == 0 || (nclust < 0 && nclust != -1)){
        fprintf(stderr, "ERROR: maximum number of clusters must either be a positive number or\n");
        fprintf(stderr, "-1 (for no limit)\n");
        exit(1);
    }
    // If haplotypes are already constructed, we can treat a barcode list as a filter.
    // This is because we don't need to cluster/toss out variants anyway, so no need
    // loading the full set of barcodes.
    if (hapsfile_given && varsfile_given && barcodesfilename.length() > 0 && 
        barcodes_filter_filename.length() == 0){
        barcodes_filter_filename = barcodesfilename;
        barcodesfilename = "";
    }

    set<unsigned long> bc_whitelist;
    bool has_bc_whitelist = false;
    if (barcodes_filter_filename.length() > 0){
        parse_barcode_file(barcodes_filter_filename, bc_whitelist);
        if (bc_whitelist.size() > 0){
            has_bc_whitelist = true;
        }
    }
    set<unsigned long> bc_filter_assn;
    bool has_bc_filter_assn = false;
    if (barcodesfilename.length() > 0){
        parse_barcode_file(barcodesfilename, bc_filter_assn);
        if (bc_filter_assn.size() > 0){
            has_bc_filter_assn = true;
        }
    }
    
    // How should we represent zero/one probability in binomial dist?
    double one = 0.999;

    // Get an ordered list of variable sites on the mitochondrial sequence
    deque<varsite> vars;
    map<int, varsite> vars_unfiltered;
    
    if (varsfile_given){
        fprintf(stderr, "Loading variants from %s...\n", varsfile.c_str());
        mito_chrom = load_vars_from_file(varsfile, vars); 
    }
    else{
        fprintf(stderr, "Finding variable sites on the mitochondrial genome...\n");
        find_vars_in_bam(bamfile, mito_chrom, minmapq, minbaseq, vars, 
            has_bc_whitelist, bc_whitelist, cov_filt); 
    }
    
    // If loading previously-inferred clusters, did the user provide
    // names for those clusters? If so, load them. If not, they will
    // get numeric IDs.
    vector<string> clust_ids;
    if (ids_given){
        parse_idsfile(idsfile, clust_ids);
    }

    nvars = vars.size();
    if (nvars > MAX_SITES){
        // Sort vars by decreasing frequency and keep top allowable number of them.
        vector<pair<float, int> > vs;
        for (deque<varsite>::iterator v = vars.begin(); v != vars.end(); ++v){
            vs.push_back(make_pair(-v->freq2, v->pos));
        } 
        sort(vs.begin(), vs.end());
        set<int> pos_rm;
        for (int i = MAX_SITES; i < vs.size(); ++i){
            pos_rm.insert(vs[i].second);
        }
        for (deque<varsite>::iterator v = vars.begin(); v != vars.end(); ){
            if (pos_rm.find(v->pos) != pos_rm.end()){
                v = vars.erase(v);
            }
            else{
                ++v;
            }
        }
        
        // IMPORTANT: this sets the size of arrays in varcounts
        nvars = vars.size();
    }
    
    // Now need to go back through the BAM file. This time, 
    // look at individual barcodes and count reads overlapping each variant.
    
    // Store barcodes as bit strings of length 2*nbases, interpreted as unsigned 
    // longs to save space
    robin_hood::unordered_map<unsigned long, var_counts> hap_counter;
    
    // Make a copy that will stay intact (we will iterate destructively) 
    deque<varsite> vars2 = vars;

    // Get counts of each variant allele for each cell barcode
    
    // Necessary for filtering variants (if not provided),
    // inferring cluster haplotypes (if not provided),
    // and assigning barcodes to individual IDs
    
    fprintf(stderr, "Counting alleles at variable sites in BAM file %s...\n", bamfile.c_str());

    count_vars_barcodes(bamfile, mito_chrom, minmapq, vars, 
        has_bc_whitelist, bc_whitelist, hap_counter);     
   
    // Still need to filter variant sites based on coverage across cells
    hapstr mask_global;
    // Also translate variant counts into haplotype strings
    robin_hood::unordered_map<unsigned long, hap> haplotypes;
    
    if (dump){
        // We will be dumping variants before filtering.
        write_vars(mito_chrom, output_prefix, vars2, mask_global, false);
    }
    
    map<int, int> orig_to_collapsed;
    map<int, set<int> > collapsed_to_orig;
    map<int, robin_hood::unordered_set<unsigned long> > site_minor;
    map<int, robin_hood::unordered_set<unsigned long> > site_major;
    map<int, robin_hood::unordered_set<unsigned long> > site_mask;
    vector<pair<long int, int> > clsort;
   
    // Define output file for cell haplotypes 
    string haps_out = output_prefix + ".cellhaps";
    
    if (!mixing_proportions){
        
        fprintf(stderr, "Building cell haplotypes from allele counts...\n");
        process_var_counts(hap_counter, haplotypes, varsfile_given,
            hapsfile_given || dump, mask_global, nvars, has_bc_whitelist, site_minor,
            site_major, site_mask, orig_to_collapsed, collapsed_to_orig,
            clsort, one);  
        
        if (dump){
            write_bchaps(haps_out, nvars, mask_global, haplotypes, barcode_group, 
                cellranger, seurat, underscore); 
            return 0; // finished
        }
    }
    else{
        mask_global.reset();
        for (int x = 0; x < nvars; ++x){
            mask_global.set(x);
        }
    }

    // Store haplotypes of clusters (whether inferred or loaded)
    vector<hap> clusthaps;
    
    FILE* statsfilef = NULL;
    
    int nclust_model = -1;
    double llrsum_model = 0.0;

    if (!hapsfile_given){
        
        // Need to infer clusters.

        fprintf(stderr, "Inferring clusters...\n");
        
        bool exact_matches_only = false;

        pair<int, float> results = infer_clusters(mask_global,
            haplotypes, nvars, clsort, collapsed_to_orig,
            site_minor, site_major, clusthaps, exact_matches_only,
            nclust, site_mask[clsort[0].second], hap_counter, one);

        nclust_model = results.first;
        llrsum_model = results.second;
       
        // Write haplotypes to output file
        string clusthapsfile = output_prefix + ".haps";
        write_clusthaps(clusthapsfile, mask_global, clusthaps); 
    }
    else{
        // Load clusters.
        load_clusthaps(hapsfilename, mask_global, clusthaps);
    }
    
    // Avoid writing out vars file only if it would overwrite a preexisting vars file.
    if (!varsfile_given && (output_prefix + ".vars" != varsfile)){ 
        // Write variants file
        write_vars(mito_chrom, output_prefix, vars2, mask_global, true); 
    }
    
    // Fill in numeric ID strings, if IDs file not provided
    if (clust_ids.size() == 0){
        char idbuf[50];
        for (int i = 0; i < clusthaps.size(); ++i){
            sprintf(&idbuf[0], "%d", i);
            clust_ids.push_back(idbuf);
        }
    }
    
    if (mixing_proportions){
        
        // We have everything we need to do this now.
        robin_hood::unordered_map<unsigned long, int> assignments;
        parse_assignments(assnfile, clust_ids, assignments, true);
        if (assignments.size() == 0){
            fprintf(stderr, "ERROR: inferring mixture proportions, but no doublets \
were found in input file %s\n", assnfile.c_str());
            exit(1); 
        }
        robin_hood::unordered_map<unsigned long, double> mixprops_mean;
        robin_hood::unordered_map<unsigned long, double> mixprops_sd;
        map<int, double> id_mixprop_mean;
        map<int, double> id_mixprop_sd;
        infer_mixprops(hap_counter, assignments, clusthaps, mask_global, nvars, mixprops_mean,
            mixprops_sd, id_mixprop_mean, id_mixprop_sd, clust_ids);
        
        // Spill to disk.
        string mixprops_out_name = output_prefix + ".props";
        FILE* mixprops_out = fopen(mixprops_out_name.c_str(), "w");
        write_mixprops(mixprops_out, barcode_group, cellranger, seurat, underscore, 
            mixprops_mean, mixprops_sd, assignments, clust_ids);
        fclose(mixprops_out);
        return 0;
    }

    // Write barcode haps file
    write_bchaps(haps_out, nvars, mask_global, haplotypes, barcode_group, cellranger, seurat, underscore);
    
    // Assign all cell barcodes to a haplotype ID.
    map<unsigned long, int> bc2hap;
    
    map<int, map<int, float> > sites_doublet_probs;
    
    robin_hood::unordered_map<unsigned long, int> assignments;
    robin_hood::unordered_map<unsigned long, double> assignments_llr;

    // oversight -- need to convert type of set here
    robin_hood::unordered_set<unsigned long> cell_filter;
    if (has_bc_filter_assn){
        for (set<unsigned long>::iterator cell = bc_filter_assn.begin();
            cell != bc_filter_assn.end(); ++cell){
            cell_filter.insert(*cell);
        }
    }
    assign_bcs(hap_counter, assignments, assignments_llr, clusthaps,
        mask_global, nvars, doublet_rate, has_bc_filter_assn, cell_filter, 
        one);
    
    map<int, int> id_counter;
    int tot_cells = 0;
    int doub_cells = 0;
    
    string assn_out = output_prefix + ".assignments";
    
    
    write_assignments(assn_out, assignments, assignments_llr,
        barcode_group, cellranger, seurat, underscore,
        clust_ids, clusthaps.size(), id_counter, 
        tot_cells, doub_cells);
   
    string statsfilename = output_prefix + ".summary";
    
    write_statsfile(statsfilename, output_prefix, nclust_model, 
        llrsum_model, tot_cells, doub_cells, id_counter, 
        doublet_rate, hapsfilename, varsfile, clust_ids);

    return 0;
}
