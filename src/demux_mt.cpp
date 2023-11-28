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
#include <htswrapper/bc_hash.h>
#include <mixtureDist/functions.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include "robin_hood.h"
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
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   ===== Options for both run modes =====\n");
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
    fprintf(stderr, "   --barcode_group -g A string to append to each barcode in the \n");
    fprintf(stderr, "       output files. This is equivalent to the batch ID appended to \n");
    fprintf(stderr, "       cell barcodes by scanpy when concatenating multiple data sets.\n");
    fprintf(stderr, "       Cell Ranger software automatically appends \"1\" here, so specify\n");
    fprintf(stderr, "       1 as an argument if you desire default, single-data set Cell\n");
    fprintf(stderr, "       -like output.\n");
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
    fprintf(stderr, "   |===== Inferring clusters from a BAM file =====|\n");
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
    fprintf(stderr, "   --mapq -q The minimum map quality filter (OPTIONAL; default 20)\n");
    fprintf(stderr, "   --baseq -Q The minimum base quality filter (OPTIONAL; default 20)\n");
    fprintf(stderr, "   ---------- General options ----------\n"); 
    fprintf(stderr, "   --nclustmax -n (OPTIONAL) If you know how many individuals should\n");
    fprintf(stderr, "       be in the mixture, you can specify that number here, and only\n");
    fprintf(stderr, "       up to that many haplotypes will be inferred. It's possible that\n");
    fprintf(stderr, "       fewer than that number of haplotypes will be found, however.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   |===== Assigning cells to previously-inferred haplotypes =====|\n");
    fprintf(stderr, "   ---------- I/O Options ----------\n");
    fprintf(stderr, "   --haps -H Cluster haplotypes from a previous run. Should be that\n");
    fprintf(stderr, "       run's [output_prefix].haps. (REQUIRED)\n");
    fprintf(stderr, "   --vars -v Variants used to build cluster haplotypes in previous \n");
    fprintf(stderr, "       run. Should be that run's [output_prefix].vars. (REQUIRED)\n");
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
            str2bc(bc_str.c_str(), as_bitset, 16);
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
 * Returns optimal number of clusters rather than cutoff.
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

void derivative(map<double, double>& dat, map<double, double>& datp){
    // Store half-open intervals
    bool setfirst = false;
    double prevX = -1;
    double prevY = -1;
    for (map<double, double>::iterator d = dat.begin(); d != dat.end(); ++d){
        if (setfirst){
            double slope = (d->second-prevY)/(d->first-prevX);
            // Interval has two ends
            if (datp.count(prevX) == 0){
                datp.insert(make_pair(prevX, slope));
            }
            else{
                datp[prevX] = 0.5*datp[prevX] + 0.5*slope;
            }
            if (datp.count(d->first) == 0){
                datp.insert(make_pair(d->first, slope));
            }
            else{
                datp[d->first] = 0.5*datp[d->first] + 0.5*slope;
            }
        }
        prevX = d->first;
        prevY = d->second;
        setfirst = true;
    }
}

/**
 * Returns x val with maximum curvature.
 */
double curvature(map<double, double>& x,
    map<double, double>& k){
    
    map<double, double> xprime1;
    derivative(x, xprime1);
    map<double, double> xprime2;
    derivative(xprime1, xprime2);

    double maxk = -1;
    double maxk_x = -1;

    for (map<double, double>::iterator xp = xprime1.begin(); xp != 
        xprime1.end(); ++xp){
        
        double xk = abs(xprime2[xp->first]) / 
            pow(1 + pow(xp->second, 2), 1.5);
        k.insert(make_pair(xp->first, xk));
        
        if (xk > maxk){
            maxk = xk;
            maxk_x = xp->first;
        }
    }

    return maxk_x;
}

void collapse_sites(hapstr& mask_global,
    int nvars,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_mask,
    map<int, int>& orig_to_collapsed,
    map<int, set<int> >& collapsed_to_orig,
    vector<pair<long int, int> >& clsort,
    bool keep_all_vars){
    
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
            if (dbinom(cladesize_A, cladesize_B, 0.001) > 
                dbinom(cladesize_A, cladesize_B, 0.5) &&
                dbinom(cladesize_A, cladesize_B, 0.999)){
                if (!keep_all_vars){
                    // Remove site? (risky)
                    //mask_global.reset(idx2);
                    //rm_err.insert(idx2);
                }
            }
            else{
                // Two sites are same haplotype if A&B == A == B
                if (dbinom(cladesize_A, cladesize_both, 0.999) > 
                    dbinom(cladesize_A, cladesize_both, 0.5) && 
                    dbinom(cladesize_A, cladesize_both, 0.001) &&
                    dbinom(cladesize_B, cladesize_both, 0.999) > 
                    dbinom(cladesize_B, cladesize_both, 0.5) && 
                    dbinom(cladesize_B, cladesize_both, 0.001)){
                    
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
    vector<pair<long int, int> >& clsort){

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
                    double dmaj = dbinom(maj+min, min, 0.001);
                    double dhet = dbinom(maj+min, min, 0.5);
                    double dmin = dbinom(maj+min, min, 0.999);
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
        collapse_sites(mask_global, nvars, site_minor, site_major,
            site_mask, orig_to_collapsed, collapsed_to_orig, sitesort, keep_all_vars);
    
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
    vector<hapstr>& haps_final, 
    hapstr& mask_global,
    int nvars,
    double doublet_rate,
    bool use_filter,
    robin_hood::unordered_set<unsigned long>& cell_filter){
     
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

        int tot_reads = 0;
        int missing_sites = 0;
        bool skip_bc = false;
        
        for (int site = 0; site < nvars; ++site){
            if (mask_global[site]){
                
                // Retrieve major/minor allele counts
                int count1 = hc->second.counts1[site];
                int count2 = hc->second.counts2[site];
                
                // var set in hapstr == count2 is expected
                for (int i = 0; i < all_model_idx.size()-1; ++i){
                    int idx1 = all_model_idx[i];
                    // how many minor alleles?
                    int nmin1 = 0;
                    // how many total?
                    int ntot1 = 1;
                    if (idx1 >= haps_final.size()){
                        pair<int, int> comb = idx_to_hap_comb(idx1, haps_final.size());
                        ntot1 = 2;
                        if (haps_final[comb.first][site]){
                            nmin1++;
                        }
                        if (haps_final[comb.second][site]){
                            nmin1++;
                        }
                    }
                    else{
                        if (haps_final[idx1][site]){
                            nmin1++;
                        }
                    }
                    for (int j = i + 1; j < all_model_idx.size(); ++j){
                        int idx2 = all_model_idx[j];
                        int nmin2 = 0;
                        int ntot2 = 1;
                        if (idx2 >= haps_final.size()){
                            pair<int, int> comb = idx_to_hap_comb(idx2, haps_final.size());
                            ntot2 = 2;
                            if (haps_final[comb.first][site]){
                                nmin2++;
                            }
                            if (haps_final[comb.second][site]){
                                nmin2++;
                            }
                        }
                        else{
                            if (haps_final[idx2][site]){
                                nmin2++;
                            }
                        }
                        
                        double exp_frac1 = (double)nmin1/(double)ntot1;
                        double exp_frac2 = (double)nmin2/(double)ntot2;

                        if (exp_frac1 != exp_frac2){
                            
                            if (exp_frac1 == 0){
                                exp_frac1 = 0.001;
                            }
                            else if (exp_frac1 == 1){
                                exp_frac1 = 0.999;
                            }
                            if (exp_frac2 == 0){
                                exp_frac2 = 0.001;
                            }
                            else if (exp_frac2 == 1){
                                exp_frac2 = 0.999;
                            }

                            double ll1 = dbinom(count1+count2, count2, exp_frac1);
                            double ll2 = dbinom(count1+count2, count2, exp_frac2);

                            llrs[idx1][idx2] += (ll1-ll2);

                        }
                    }
                }
            }
        }
        
        double llr;
        int best_assignment = collapse_llrs(llrs, llr);
        if (llr > 0){
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
void load_vars_from_file(string& filename, deque<varsite>& vars){
    ifstream f(filename);
    string line;
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
    set<unsigned long>& bc_whitelist){
    
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
            }
        }
    }
    
    if (n < 0){
        fprintf(stderr, "bam_mplp_auto failed for %s\n", infile.fname);
        exit(1);
    }
    bam_mplp_destroy(plp);
    
    vector<pair<double, int> > vs_sort;
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != 
        vars_unfiltered.end(); ++v){
        vs_sort.push_back(make_pair(-(double)v->second.freq2, v->first));
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
    reader.set_query_region(mito_chrom.c_str(), -1, -1);
    
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

            if (str2bc(bc_str.c_str(), bc_hashable, 16) && (
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
    vector<hapstr>& haps_final,
    bool exact_matches_only,
    int nclust_max,
    robin_hood::unordered_set<unsigned long>& cellset,
    robin_hood::unordered_map<unsigned long, var_counts>& hap_counter){
    
    hapstr mask;

    vector<hapstr> hapsites;
    vector<set<unsigned long> > hapgroups;
    vector<set<unsigned long> > hapgroups_not;
     
    vector<hapstr> hapsites_prev;
    vector<set<unsigned long> > hapgroups_prev;
    vector<set<unsigned long> > hapgroups_not_prev;

    int nsites_included = 0;
    vector<int> sites_order;
    vector<vector<hapstr> > haps_final_order;
    vector<hapstr> mask_order;
    vector<hapstr> haps_final_prev;
    
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
                haps_final.push_back(hapsites[sizepairs[i].second]);
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
                    if (haps_final_prev[i] != haps_final[i]){
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
    assign_bcs(hap_counter, assn, assn_llr, haps_final, mask, nvars, 0.0, true, cellset);
    for (robin_hood::unordered_map<unsigned long, double>::iterator al = assn_llr.begin();
        al != assn_llr.end(); ++al){
        llrsum += al->second;
    }
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
                    mask_order[site_idx], nvars, 0.0, true, cellset);
                for (robin_hood::unordered_map<unsigned long, double>::iterator al = 
                    assn_llr.begin(); al != assn_llr.end(); ++al){
                    llrsum_new += al->second;
                }
                if (llrsum_new < llrsum){
                    // Finished
                    haps_final = haps_final_order[site_idx+1];
                    mask = mask_order[site_idx+1];
                    break;
                }
                else{
                    if (haps_final_order[site_idx].size() != 
                        haps_final_order[site_idx+1].size()){
                        fprintf(stderr, "%ld haplotypes after removing %d sites\n", 
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
    string id;
    while (infile >> id ){
        ids.push_back(id);
    }
}

void write_bchaps(string& haps_out,
    int nvars,
    hapstr& mask_global,
    robin_hood::unordered_map<unsigned long, hap>& haplotypes,
    string& bc_group){
    
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
        string bcstr = bc2str(bcbits, 16);
        if (bc_group != ""){
            bcstr += "-" + bc_group;
        }

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
    vector<hapstr>& clusthaps){
    
    FILE* clusthaps_f = fopen(clusthapsfile.c_str(), "w");
    for (int i = 0; i < clusthaps.size(); ++i){
        for (int x = 0; x < nvars; ++x){
            if (mask_global[x]){
                if (clusthaps[i].test(x)){
                    fprintf(clusthaps_f, "1");
                }
                else{
                    fprintf(clusthaps_f, "0");
                }
            }
        }
        fprintf(clusthaps_f, "\n");
    }
    fclose(clusthaps_f);
}

void load_clusthaps(string& hapsfilename,
    hapstr& mask_global,
    vector<hapstr>& clusthaps){

    ifstream infile(hapsfilename.c_str());
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
        hapstr h;
        for (int i = 0; i < nvars; ++i){
            if (hapstring[i] == '0'){
            }
            else if (hapstring[i] == '1'){
                h.set(i);
            }
        }
        clusthaps.push_back(h);
    }
}

/**
 * Print barcode-identity assignments to disk.
 */
void write_assignments(string& assn_out,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    string& barcode_group,
    vector<string>& clust_ids,
    int nhaps,
    map<string, int>& id_counter,
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
        if (id_counter.count(name) == 0){
            id_counter.insert(make_pair(name, 0));
        }
        id_counter[name]++;
        bc as_bitset(assn->first);
        string bc_str = bc2str(as_bitset, 16);
        if (barcode_group != ""){
            bc_str += "-" + barcode_group;
        }
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
    map<string, int>& id_counter,
    double doublet_rate,
    string& hapsfilename,
    string& varsfilename){

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
    }
    vector<pair<int, string> > idcounts_sorted;
    for (map<string, int>::iterator idc = id_counter.begin(); idc != 
        id_counter.end(); ++idc){
        idcounts_sorted.push_back(make_pair(-idc->second, idc->first));
    }
    sort(idcounts_sorted.begin(), idcounts_sorted.end());
    for (vector<pair<int, string> >::iterator idcs = idcounts_sorted.begin(); 
        idcs != idcounts_sorted.end(); ++idcs){
        fprintf(statsfilef, "%s\t%s\t%d\n", output_prefix.c_str(), 
            idcs->second.c_str(), -idcs->first);   
    }
    fclose(statsfilef);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"barcodes_filter", required_argument, 0, 'f'},
       {"barcode_group", required_argument, 0, 'g'},
       {"mapq", required_argument, 0, 'q'},
       {"baseq", required_argument, 0, 'Q'},
       {"nclust_max", required_argument, 0, 'n'},
       {"mito", required_argument, 0, 'm'},
       {"dump", no_argument, 0, 'd'},
       {"vars", required_argument, 0, 'v'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"haps", required_argument, 0, 'H'},
       {"ids", required_argument, 0, 'i'},
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

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:B:f:g:q:Q:n:m:v:H:i:D:dh", long_options, &option_index )) != -1){
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
            case 'g':
                barcode_group = optarg;
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

    // Get an ordered list of variable sites on the mitochondrial sequence
    deque<varsite> vars;
    map<int, varsite> vars_unfiltered;
    
    if (varsfile_given){
        fprintf(stderr, "Loading variants from %s...\n", varsfile.c_str());
        load_vars_from_file(varsfile, vars); 
    }
    else{
        fprintf(stderr, "Finding variable sites on the mitochondrial genome...\n");
        find_vars_in_bam(bamfile, mito_chrom, minmapq, minbaseq, vars, 
            has_bc_whitelist, bc_whitelist); 
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
    
    fprintf(stderr, "Counting alleles at variable sites in BAM...\n");
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
    
    fprintf(stderr, "Building cell haplotypes from allele counts...\n");

    map<int, int> orig_to_collapsed;
    map<int, set<int> > collapsed_to_orig;
    map<int, robin_hood::unordered_set<unsigned long> > site_minor;
    map<int, robin_hood::unordered_set<unsigned long> > site_major;
    map<int, robin_hood::unordered_set<unsigned long> > site_mask;
    vector<pair<long int, int> > clsort;
    process_var_counts(hap_counter, haplotypes, varsfile_given,
        hapsfile_given || dump, mask_global, nvars, has_bc_whitelist, site_minor,
        site_major, site_mask, orig_to_collapsed, collapsed_to_orig,
        clsort);  
   
    
    // Define output file for cell haplotypes 
    string haps_out = output_prefix + ".cellhaps";
    if (dump){
        write_bchaps(haps_out, nvars, mask_global, haplotypes, barcode_group); 
        return 0; // finished
    }

    // Store haplotypes of clusters (whether inferred or loaded)
    vector<hapstr> clusthaps;
    
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
            nclust, site_mask[clsort[0].second], hap_counter);

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
    if (!varsfile_given || (output_prefix + ".vars" != varsfile)){ 
        // Write variants file
        write_vars(mito_chrom, output_prefix, vars2, mask_global, true); 
    }
    
    // Write barcode haps file
    write_bchaps(haps_out, nvars, mask_global, haplotypes, barcode_group);
    
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
        mask_global, nvars, doublet_rate, has_bc_filter_assn, cell_filter);
    
    map<string, int> id_counter;
    int tot_cells = 0;
    int doub_cells = 0;
    
    string assn_out = output_prefix + ".assignments";

    write_assignments(assn_out, assignments, assignments_llr,
        barcode_group, clust_ids, clusthaps.size(), id_counter, 
        tot_cells, doub_cells);
   
    string statsfilename = output_prefix + ".summary";
    
    write_statsfile(statsfilename, output_prefix, nclust_model, 
        llrsum_model, tot_cells, doub_cells, id_counter, 
        doublet_rate, hapsfilename, varsfile);

    return 0;
}
