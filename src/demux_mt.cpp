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
    fprintf(stderr, "Given one or more BAM files representing multiple individuals from \
the same species, finds variable positions in the mitochondrial DNA and determines the \
mitochondrial haplotype of each cell barcode in the BAMs. Then outputs a list of cell barcode -> \
mtDNA cluster identities.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --bam -b The BAM file of interest. If more than one BAM is desired, subset all to \
mitochondrial sequences (i.e. samtools view -bh <file> chrM) and merge (samtools merge), then index (samtools \
index).\n");
    fprintf(stderr, "    --output_prefix -o The prefix of output files to create. Creates <output_prefix>.vars, \
<output_prefix>.allhaps, <output_prefix>.clusthaps, and <output_prefix>.assignments by default. If dump data (-d) \
is specified, creates only <output_prefix>.vars, <output_prefix>.allhaps, and <output_prefix>.clusthaps. If loading \
data from a previous run, omits all files that were loaded instead of created.\n");
    fprintf(stderr, "    --barcodes -B (OPTIONAL) list of whitelisted barcodes (i.e. from CellRanger). Will exclude non-cell barcodes.\n");
    fprintf(stderr, "    --mapq -q The minimum map quality required (20)\n");
    fprintf(stderr, "    --baseq -Q The minimum base quality required (20)\n");
    fprintf(stderr, "    --freq -f The minimum frequency an allele must appear to be counted \
as genuine and not an error (0.005)\n");
    fprintf(stderr, "    --exact -e When clustering mitochondrial haplotypes, only use cells that contain \
all alleles in a growing haplotype (faster and better for good/high coverage data, worse for low-coverage)\n");
    fprintf(stderr, "    --num -n (OPTIONAL) expected number of clusters. If not provided, it will be inferred\n");
    fprintf(stderr, "    --mincount -C minimum count for edge to be considered (heuristic to speed up; default 5)\n");
    fprintf(stderr, "    --minsites -s The minimum number of sites visible to consider a cell. \
This will automatically be adjusted upward as necessary but cannot be adjusted downward. Higher \
numbers will speed up execution, but too high a number will unnecessarily discard data.\n");
    fprintf(stderr, "    --cov -c Minimum coverage to count a site as variable (OPTIONAL; default 10)\n");
    fprintf(stderr, "    --mito -m The name of the mitochondrial sequence (OPTIONAL; default: chrM)\n");
    fprintf(stderr, "    --missing_thresh -M Minimum fraction missing sites to accept a haplotype (default = 0.25)\n");
    fprintf(stderr, "    --dump -d Dump variant and haplotype (barcode and cluster) info to output files and quit. \
Skips assigning cells to barcodes. One use of this feature would be to run the program once with this option \
and -B to infer cluster haplotypes using only barcodes in a list of confidently-assigned cells (i.e. from CellRanger), \
then to run a second time, loading cluster haplotypes (-H) and variant sites (-v) from this first confident run, \
but assigning all barcodes to inferred individuals (omitting -d and -B).\n");
    fprintf(stderr, "    --dump_raw -r Like --dump/-d, but outputs raw counts, including intermediate frequency. \
Useful for visualizing doublets.\n");
    fprintf(stderr, "    --tree -t (OPTIONAL) print a Newick tree showing clustering of sites, for debugging purposes\n");
    fprintf(stderr, "    --vars -v (OPTIONAL) provide a pre-determined list of variant sites here. These will still \
be filtered for coverage, but these variants will be used instead of attempting to discover new ones.\n");
    fprintf(stderr, "    --haps -H If haplotypes were pre-computed, load them here. -vars / -v also required.\n");
    fprintf(stderr, "    --nodoublet -N Disable checking for doublets. Will increase performance if there are many individuals.\n");
    fprintf(stderr, "    --gex -g Specify that this is gene expression (single-cell RNA-seq) data\n");
    fprintf(stderr, "    --ids -i (OPTIONAL) If using pre-computed haplotypes and variants (-H and -v provided), \
provide an optional file of individual IDs, with line indices matching haplotypes in the -H file. Assignments will \
then use these IDs instead of numeric indices.\n");
    fprintf(stderr, "    --identity_thresh -I what percent different are sites allowed to be before they are \
collapsed? -1 == find automatically (SLOW), 0 == do not collapse. Default = 0.07 for GEX, 0 otherwise.\n");
    fprintf(stderr, "    --track_sites_causing_doublets -S If seeing more doublets than expected, this could be \
caused by a small number of weird sites that tend to get 50 percent coverage in many cells (maybe low mappability regions?) This option collects and outputs information about each site's contribution to each type of doublet call.\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
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
int choose_nclust(vector<int>& clust_sizes, 
    double& llr_return, 
    double& exp_size1_return, 
    double& exp_size2_return, 
    int& totcells_return, 
    int& smallest_group_return){
    
    // Sort and compute total
    int nclust_tot = 0;
    vector<int> clust_sorted;
    for (vector<int>::iterator c = clust_sizes.begin(); c != clust_sizes.end(); ++c){
        clust_sorted.push_back(-*c);
        nclust_tot += *c;
    }
    sort(clust_sorted.begin(), clust_sorted.end());
    
    int nclust_max = clust_sorted.size();
    
    // Store the highest log likelihood encountered
    double loglik_max = 0.0;
    // Store the second-highest log likelihood encountered
    // We will break iteration when we encounter a dip in LL,
    // so second best is either the last LL computed or 
    // two before the last computed.
    double loglik_second = 0.0;
    int nclust_best = -1;
    int minsize = -1;
    
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
        
        fprintf(stderr, "nclust %d LL %f means %f %f\n", nclust, loglik, meansize1, meansize2);

        if (nclust_best == -1 || loglik > loglik_max){
            // Update second likeliest to previous
            loglik_second = loglik_max;

            nclust_best = nclust;
            loglik_max = loglik;
            
            exp_size1_return = meansize1;
            exp_size2_return = meansize2;
            
            totcells_return = ncell_this;
            smallest_group_return = -clust_sorted[nclust-1]; 
        }
        else if (nclust_best != -1 && loglik < loglik_max){
            if (loglik > loglik_second){
                loglik_second = loglik;
            }
            break;
        }

    }
    
    // Return the best number of clusters; update the log likelihood ratio
    // parameter to the ratio between best and second best choice for
    // number of clusters.

    if (loglik_second == 0 ){
        llr_return = 0;
    }
    else{
        llr_return = loglik_max - loglik_second;
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
        //if (vars_blacklist.find(idx) == vars_blacklist.end()){
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

/**
 *  Go through counts of each variant in each individual and do several things:
 *  Convert read counts to bitset haplotypes.
 *  Compute mean and SD allele frequency for all variable sites.
 *  Blacklist variable sites with too high SD - create global bitset mask
 */
void process_var_counts(robin_hood::unordered_map<unsigned long, var_counts>& hap_counter,
    robin_hood::unordered_map<unsigned long, hap>& bc2hap,
    set<int>& sites_blacklist,
    bool varsfile_provided,
    hapstr& mask_global,
    int nvars,
    bool gex,
    bool barcodes_given){

    // Set up global mask
    mask_global.reset();
    for (int i = 0; i < nvars; ++i){
        mask_global.set(i);
    }
    
    int hetmiss_max = MAX_SITES;

    // Cells that have very low coverage or a lot of sites that look 
    // heterozygous are probably messed up.
    // These might interfere with clustering or at least slow it down.
    
    // While we will avoid clustering with these, we will save their var_counts
    // so they can still be identified later, after clustering.
    
    map<double, double> miss_filter_hist;
    for (int i = 1; i < MAX_SITES; ++i){
        miss_filter_hist.insert(make_pair(i, 0));
    }
    
    // Track how many sites are either missing or look heterozygous
    // in each cell
    robin_hood::unordered_map<unsigned long, int> cell2hetmiss;

    for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
        hap_counter.begin(); hc != hap_counter.end(); ++hc){
        int hetcount = 0;
        int misscount = 0;
        int totcount = 0;
        for (int i = 0; i < nvars; ++i){
            
            int count1 = hc->second.counts1[i];
            int count2 = hc->second.counts2[i];
            
            totcount += count1 + count2;
            if (count1+count2 == 0){
                ++misscount;
            }
            else{
                double dhom1 = dbinom(count1+count2, count2, 0.001);
                double dhom2 = dbinom(count1+count2, count2, 0.999);
                double dhet = dbinom(count1+count2, count2, 0.5);

                if (dhet > dhom1 && dhet > dhom2){
                    ++hetcount;
                }
            }
        }
        
        // Missing and heterozygous-looking sites are both
        // unusable for clustering.

        for (int i = misscount + hetcount; i < MAX_SITES; ++i){
            miss_filter_hist[i]++;
        }

        cell2hetmiss.emplace(hc->first, misscount + hetcount);
    }
    
    if (!gex){

        // What we should have now is a curve of how many cells are included
        // with each possible cutoff of missing sites.

        // The curve should start out mostly flat, then shoot up as we allow
        // cells with nearly all missing sites. The flat part represents
        // true cells, and the curve at the end comes from including empty
        // droplets that contain noise.
        
        // For clustering, we want to exclude these empty droplets
        // (although we can still seek to identify them later, once
        // haplotypes have been learned).

        // We will approximate the first and second derivatives of this
        // curve, and then choose as a cutoff the point with the highest
        // curvative, defined as abs(f"(x)) / (1 + (f'(x))^2)^(3/2)
        
        map<double, double> miss_filter_hist_curvature;
        double maxcurv_cutoff = curvature(miss_filter_hist, 
            miss_filter_hist_curvature);
        
        // This cutoff is the maximum allowable number of missing + het
        // sites in a cell.
        hetmiss_max = (int)round(maxcurv_cutoff);
    }

    // Now, with the cells that pass filter, do the same thing but to find
    // a minimum number of cells visible per site. If we do not wish to
    // filter variants (i.e. if a variant file was provided), we can skip
    // this step.
    
    // We also omit this step for RNA-seq data, because in that case, we
    // will attempt to collapse sites with similar sets of cells, which 
    // will rescue many sites that would otherwise be deleted here.
    
    if (!varsfile_provided && !gex){
        // Now look to see what coverage & MAF look like for each site
        // in remaining cells.
        map<double, double> site_filter_hist;
        
        vector<int> site_ncells;

        for (int i = 0; i < nvars; ++i){
            int cells_hetmiss = 0;
            int cells_good = 0;
            for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
                hap_counter.begin(); hc != hap_counter.end(); ++hc){
            
                if (cell2hetmiss[hc->first] <= hetmiss_max){
                
                    int count1 = hc->second.counts1[i];
                    int count2 = hc->second.counts2[i];
                    
                    bool hetmiss = false; 
                    if (count1+count2 == 0){
                        hetmiss = true;
                    }
                    else{
                        double dhom1 = dbinom(count1+count2, count2, 0.001);
                        double dhom2 = dbinom(count1+count2, count2, 0.999);
                        double dhet = dbinom(count1+count2, count2, 0.5);

                        if (dhet > dhom1 && dhet > dhom2){
                            hetmiss = true;
                        }
                    }
                    if (hetmiss){
                        ++cells_hetmiss;
                    }
                    else{
                        ++cells_good;
                    }
                }
            }
            site_ncells.push_back(cells_good);
            for (double i = 1; i <= cells_good; ++i){
                if (site_filter_hist.count(i) == 0){
                    site_filter_hist.insert(make_pair(i, 0));
                }
                site_filter_hist[i]++;
            }
        }
        
        map<double, double> site_filter_hist_curv;
        double sfcutoff = curvature(site_filter_hist, site_filter_hist_curv);
        int ncells_min = (int)round(sfcutoff);
        
        // Blacklist sites that don't pass this cutoff.
        for (int i = 0; i < nvars; ++i){
            if (site_ncells[i] < ncells_min){
                mask_global.reset(i);
            }
        }
        fprintf(stderr, "%ld sites passed filtering\n", mask_global.count());
    }
    
    // Now, with filters set, create a bitset haplotype string for any
    // cell that survived filtering, at the variable sites that passed
    // filtering.

    for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
        hap_counter.begin(); hc != hap_counter.end(); ++hc){
        
        if (cell2hetmiss[hc->first] <= hetmiss_max){    
            hap h;
            h.vars.reset();
            h.mask.reset();
            
            for (int i = 0; i < nvars; ++i){
                if (mask_global.test(i)){    
                    int count1 = hc->second.counts1[i];
                    int count2 = hc->second.counts2[i];
                    
                    if (count1 + count2 > 0){    
                        double dhom1 = dbinom(count1+count2, count2, 0.001);
                        double dhom2 = dbinom(count1+count2, count2, 0.999);
                        double dhet = dbinom(count1+count2, count2, 0.5);
                           
                        if (dhom1 > dhet && dhom1 > dhom2){
                            h.mask.set(i);
                        }
                        else if (dhom2 > dhet && dhom2 > dhom1){
                            h.mask.set(i);
                            h.vars.set(i);
                        }
                        else{
                            // Unusable in this cell.

                        }
                    }
                }
            }
            // Double check that the cell is informative before storing
            if (h.mask.count() > 0){
                bc2hap.emplace(hc->first, h);
            }
        }
    }
    fprintf(stderr, "%ld cell barcodes passed filtering\n", bc2hap.size());
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
    bool nodoublet){
     
    vector<int> all_model_idx;
    for (int i = 0; i < haps_final.size(); ++i){
        all_model_idx.push_back(i);
        if (!nodoublet){
            for (int j = i + 1; j < haps_final.size(); ++j){
                int k = hap_comb_to_idx(i, j, haps_final.size());
                all_model_idx.push_back(k);
            }
        }
    }
    sort(all_model_idx.begin(), all_model_idx.end());

    for (robin_hood::unordered_map<unsigned long, var_counts>::iterator hc = 
        hap_counter.begin(); hc != hap_counter.end(); ++hc){
        
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
                    llrs[idx1].insert(make_pair(idx2, 0.0));
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
        assignments.emplace(hc->first, best_assignment);
        assignments_llr.emplace(hc->first, llr);
    }
    
    hap_counter.clear();
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
    float minfreq,
    deque<varsite>& vars,
    bool gex,
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
                    //if (p->is_refskip){
                        n_skip++;
                    //}
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
                
                if (gex){
                    // If in RNA-seq mode, splicing can do weird things. Allow
                    // any read covering the position, even if it's a skipped
                    // position, to count toward the coverage total.
                    cov += n_skip;
                }
                
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
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != vars_unfiltered.end();
        ++v){
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
    return;


    // Now use a mixture model to try to separate true SNPs from fake
    // variation induced by mapping/sequencing error.
    vector<vector<double> > em_input;
    vector<int> em_input2var;
    double meancov = 0.0;
    double meanfreq = 0.0;
    double covmax = 0.0;
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != vars_unfiltered.end();
        ++v){
        vector<double> row;
        // Look at minor allele frequency
        if (v->second.freq2 > 0.01){
            row.push_back((double)v->second.cov);
            
            em_input.push_back(row);
            
            em_input2var.push_back(v->first);
        }
        if (!gex){
            // Also look at coverage, if not RNA-seq. RNA-seq coverage
            // is allowed to fluctuate.
            //row.push_back((double)v->second.cov);
        }
        meancov += (1.0/(double)vars_unfiltered.size()) * (double)v->second.cov;
        meanfreq += (1.0/(double)vars_unfiltered.size()) * (double)v->second.freq2;
        if (v->second.cov > covmax){
            covmax = v->second.cov;
        }
    }
    fprintf(stderr, "meancov %f meanfreq %f\n", meancov, meanfreq); 
    vector<mixtureDist> dists;
    //mixtureDist err("beta", vector<double>{0.01*10, (1-0.01)*10});
    mixtureDist err("normal", vector<double>{10, covmax/4.0});
    dists.push_back(err);
    //mixtureDist snp("beta", vector<double>{0.1*100, (1-0.1)*100});
    mixtureDist snp("normal", vector<double>{meancov, covmax/4.0});
    dists.push_back(snp);
    mixtureModel mod(dists);
    mod.fit(em_input);
    mod.print();
    exit(0);

    map<int, int> pos2idx;
    float cov_max = 0.0;
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != vars_unfiltered.end(); ++v){
        fprintf(stdout, "%d\t%f\n", v->second.cov, v->second.freq2);
        /*
        cov_all.push_back(v->second.cov);
        cov_em.push_back((float)v->second.cov);
        freq_em.push_back((float)v->second.freq2);
        //freq_em.push_back(v->second.freq2 * (float)v->second.cov);
        if ((float)v->second.cov > cov_max){
            cov_max = (float)v->second.cov;
        }
        pos2idx.insert(make_pair(v->first, cov_em.size()-1));
        */
    }
    exit(0);
    /*
    vector<float> cov_mu;
    cov_mu.push_back(0.0);
    cov_mu.push_back(cov_max);
    vector<float> cov_sd;
    cov_sd.push_back(cov_max/4.0);
    cov_sd.push_back(cov_max/4.0);
    vector<float> cov_weight;
    cov_weight.push_back(0.5);
    cov_weight.push_back(0.5);
    set<int> cov_bl;
    if (!gex){
        // Blacklist low-coverage sites
        gaussEM_blacklist(cov_em, cov_weight, cov_mu, cov_sd, cov_bl, true);
    }
    else{
        // Blacklist low-frequency sites
        gaussEM_blacklist(freq_em, cov_weight, cov_mu, cov_sd, cov_bl, true);
    }
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != vars_unfiltered.end(); ){
        if (cov_bl.find(pos2idx[v->first]) == cov_bl.end()){
            vars.push_back(v->second);
        }
        vars_unfiltered.erase(v++);
    }
    fprintf(stderr, "%ld vars passed coverage threshold\n", vars.size());
    return;
    
    
    
    sort(cov_all.begin(), cov_all.end());
    
    // Find median coverage
    
    int cov_med = -1;
    if (cov_all.size() % 2 == 0){
        cov_med = (int)round(((float)cov_all[cov_all.size() / 2 - 1] + (float)cov_all[cov_all.size() / 2]) / 2.0);
    }
    else{
        cov_med = cov_all[(cov_all.size() - 1) / 2];
    }
    cov_all.clear();

    fprintf(stderr, "MEDIAN COV: %d\n", cov_med);
    for (map<int, varsite>::iterator v = vars_unfiltered.begin(); v != vars_unfiltered.end(); ){
        //bool pass = true;
        //bool pass = v->second.cov >= 0.5*(float)cov_med && v->second.cov <= 2.0*(float)cov_med;
        bool pass = v->second.cov >= (int)round(0.5*(float)cov_med);
        if (pass){
            vars.push_back(v->second);
        }   
        vars_unfiltered.erase(v++);
    }
    */
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
    robin_hood::unordered_map<unsigned long, var_counts>& hap_counter,
    bool gex){

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
                !has_bc_whitelist || bc_whitelist.find(bc_hashable.to_ulong()) != bc_whitelist.end())){ 
                 
                // Remove any variants already done with
                while (vars.size() > 0 && vars.front().pos < reader.reference_start){
                    vars.pop_front();
                    vars_idx++;
                }

                if (vars.size() > 0 && vars.front().pos >= reader.reference_start && vars.front().pos <= reader.reference_end){
                    int vars_idx2 = 0;
                    for (deque<varsite>::iterator var = vars.begin(); var != vars.end(); ++var){
                        if (var->pos > reader.reference_end){
                            break;
                        }
                        else{
                            // Look for variant in read.
                            char base = reader.get_base_at(var->pos + 1);
                            
                            if (gex && base == '-'){
                                // Treat a skipped reference base as major allele.
                                //base = var->allele1;
                            }
                            
                            // Make sure an entry for this barcode exists.
                            if (hap_counter.count(bc_hashable.to_ulong()) == 0){
                                var_counts v;
                                //hap_counter.insert(make_pair(bc_hashable.to_ulong(), v));           
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

pair<int, float> infer_clusters2(hapstr& mask_global, 
    robin_hood::unordered_map<unsigned long, hap>& haplotypes, 
    int nvars,
    vector<pair<long int, int> >& clades_sort,
    map<int, set<int> >& collapsed_to_orig,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    vector<hapstr>& haps_final,
    bool exact_matches_only,
    int nclust_max){
    
    hapstr mask;

    vector<hapstr> hapsites;
    vector<set<unsigned long> > hapgroups;
    vector<set<unsigned long> > hapgroups_not;
     
    vector<hapstr> hapsites_prev;
    vector<set<unsigned long> > hapgroups_prev;
    vector<set<unsigned long> > hapgroups_not_prev;

    double llr_prev = 0.0;
    int nclust_prev = -1;
    int totcells_prev = -1;

    // What were the mean sizes of true and erroneous
    // clusters with the previous set of sites/haplotypes?
    double exp_size1_prev = -1;
    double exp_size2_prev = -1;
    
    // What were the mean sizes of true and erroneous
    // clusters with the first set of sites that give us
    // the current number of haplotypes? 
    double exp_size1_thisnclust = -1;
    double exp_size2_thisnclust = -1;
    
    int nsites_included = 0;
    
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

        // Now determine likelist number of haplotypes, given the current
        // set of sites & clades.
        
        // Compute the size of each currently-tracked haplotype    
        vector<int> sizes;
        // Store sizes of haplotypes mapped to haplotype indices
        vector<pair<int, int> > sizepairs;
        
        for (int i = 0; i < hapgroups.size(); ++i){
            sizes.push_back((int)hapgroups[i].size());
            sizepairs.push_back(make_pair((int)hapgroups[i].size(), i));
        }
        
        // Test the current group sizes to determine which seem real 
        // and which seem like errors. 
        // es1 is the mean/expected size of real clusters
        // es2 is the mean/expected size of erroneous clusters

        double llr; // log likelihood ratio of best to second best choice for # clusts
        double es1; // expected size of true clusters under best model
        double es2; // expected size of erroneous clusters under best model
        int totcells_model; // total number of cells in clusters under best model
        int smallest_group; // size of smallest true cluster under best model
        
        int nclust = choose_nclust(sizes, llr, es1, es2, totcells_model, smallest_group);
        
        // There's a danger in including too many sites -- if we add
        // many true sites in a row, we can gradually whittle down
        // the sizes of groups. This can result in a new, erroneous
        // site being treated as a real SNP and throwing off haplotypes
        // and group numbers.
        
        // For example, if we have 5 haplotypes with 200 cells each, and
        // we continue adding sites to these haplotypes until our group
        // membership narrows down to 30 cells each (but still the same 
        // haplotypes), a future site could make us jump from these 5 
        // haplotypes to 8, each containing only 15 cells.
        //
        // To guard against this, we store the expected sizes of true and
        // erroneous clusters from earlier iterations, when group sizes
        // were bigger. If the current mean size of a true cluster is more 
        // likely to be an error than a true cluster under an earlier set 
        // of sites that gave us the same number of haplotypes, then we 
        // have gone too far and should stop.
        //
        // This can pose a problem if we have not yet seen enough sites,
        // though. To do this heuristic check, we require that we've already
        // included at least 3 sites, and that we have included enough sites
        // to build at least 4*(the current number of true clusters)
        // haplotypes. This is because very early on, group sizes are huge
        // (i.e. with one site, we're including most cell barcodes, and we have
        // only 2 groups, both of which are probably supersets of multiple
        // haplotypes). We want to make sure group sizes have shrunk to more
        // realistic levels before doing the test.
        
        bool enough_sites_seen = (nsites_included > 2 && nsites_included > log2(4*nclust));
        
        if (nclust == nclust_prev){
            if (!enough_sites_seen){
                // If we haven't seen enough sites yet, disable the check to see
                // whether group sizes have become too small.
                exp_size1_thisnclust = es1;
                exp_size2_thisnclust = es2;
            }   
        }
        else{
            // This is the first time seeing a new, higher number of clusters.
            // Store the mean sizes of true / false clusters, so we can tell
            // when group sizes are becoming too small.
            exp_size1_thisnclust = es1;
            exp_size2_thisnclust = es2;
        }
            
        // Should we stop adding sites here?
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
        else if (nclust > nclust_prev){
            if (enough_sites_seen && nclust_prev != -1){
                // We just increased the number of clusters. Check to see
                // whether too many sites may have been included. 
                //double mean_clustsize = (float)totcells_prev / (float)nclust;
                float ll1 = dpois(smallest_group, exp_size1_prev);
                float ll2 = dpois(smallest_group, exp_size2_prev);
                if (ll2 > ll1){
                    terminate = true;
                }  
            } 
        }
        else if (nclust == nclust_prev){
            if (!terminate && nclust > 1 && 
                exp_size1_thisnclust != -1 &&
                exp_size2_thisnclust != -1){
                // We have stayed at the same number of clusters. Check to see
                // whether group sizes have fallen too far by testing whether
                // the current true cluster size is more likely to be an error
                // than a true cluster under the original model that produced
                // this number of clusters.
                float ll1 = dpois(smallest_group, exp_size1_thisnclust);
                float ll2 = dpois(smallest_group, exp_size2_thisnclust);
                if (ll2 > ll1){
                    terminate = true;
                }
            } 
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
            
            exp_size1_prev = es1;
            exp_size2_prev = es2;
            nclust_prev = nclust;
            totcells_prev = totcells_model;
            llr_prev = llr;

            haps_final.clear();
            sort(sizepairs.begin(), sizepairs.end());
            
            for (int i = sizepairs.size()-1; i >= sizepairs.size()-nclust; --i){
                haps_final.push_back(hapsites[sizepairs[i].second]);
                fprintf(stderr, "\tsize %d [ ", sizepairs[i].first);
                for (int x = 0; x < nvars; ++x){
                    if (mask[x] && hapsites[sizepairs[i].second][x]){
                        fprintf(stderr, "%d ", x);
                    }
                }
                fprintf(stderr, "]\n");
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
    fprintf(stderr, "%d clusts\n", nclust_prev);
    return make_pair(nclust_prev, llr_prev);
}

void get_clades(hapstr& mask_global,
    robin_hood::unordered_map<unsigned long, hap>& haplotypes,
    int nvars,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_mask,
    vector<pair<long int, int> >& clsort){
    
    // First, transform cell haplotypes & variant positions into
    // clades. A clade is a group of cells informed by a specific
    // variant.
    // Any given variant creates a set of cells that are covered/
    // visible at its position (clades_mask),
    // a set of cells that have the minor allele (clades)
    // and a set of cellst hat have the major allele (clades_not)
    
    for (int i = 0; i < nvars; ++i){
        if (mask_global[i]){
            
            robin_hood::unordered_set<unsigned long> cl;
            clades.insert(make_pair(i, cl));
            clades_not.insert(make_pair(i, cl));
            clades_mask.insert(make_pair(i, cl));
            
            for (robin_hood::unordered_map<unsigned long, hap>::iterator h = haplotypes.begin(); h != 
                haplotypes.end(); ++h){
                if (h->second.mask[i] & h->second.vars[i]){
                    clades[i].insert(h->first);
                }
                else if (h->second.mask[i] & !h->second.vars[i]){
                    clades_not[i].insert(h->first);
                }
                
                if (h->second.mask[i]){
                    clades_mask[i].insert(h->first);
                }
            }
            
            if (clades[i].size() == 0){
                // No information at this site; erase it.
                clades.erase(i); 
                clades_mask.erase(i);
                clades_not.erase(i);   
            }
            else{
                clsort.push_back(make_pair(-clades[i].size(), i));
            }
        }
    }
    
    // Now, sort sites in decreasing order of MAF (informative clade size)
    sort(clsort.begin(), clsort.end());
    

}

void collapse_sites(hapstr& mask_global,
    robin_hood::unordered_map<unsigned long, hap>& haplotypes,
    int nvars,
    map<int, robin_hood::unordered_set<unsigned long> >& clades,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_not,
    map<int, robin_hood::unordered_set<unsigned long> >& clades_mask,
    map<int, int>& orig_to_collapsed,
    map<int, set<int> >& collapsed_to_orig,
    vector<pair<long int, int> >& clsort){
    
    for (int i = 0; i < clsort.size()-1; ++i){
        int idx1 = clsort[i].second;
        
        // If it's already been collapsed, stop looking
        if (orig_to_collapsed.count(idx1) > 0){
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

        double errsize_tot_this = 0.0;
        double errsize_count_this = 0.0;

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

            // If clade B looks more like an error at this point than
            // it looks like A, it's safe to stop looking to merge
            // (clades will only get smaller).
            if (errsize_count_this > 0.0){
                double errsize_mean = errsize_tot_this/errsize_count_this;
                if (dpois(cladesize_B, errsize_mean) > 
                    dpois(cladesize_B, cladesize_A)){
                    break;
                }
            }

            // At this point, neither is an error.

            // Use smaller clade to represent "true" counts.
            vector<pair<int, int> > cladetype_sort;
            cladetype_sort.push_back(make_pair(-cladesize_both, 0));
            cladetype_sort.push_back(make_pair(-cladesize_AminusB, 1));
            cladetype_sort.push_back(make_pair(-cladesize_BminusA, 2));
            sort(cladetype_sort.begin(), cladetype_sort.end());
            
            // Ask how many things.
            double maxll = 0.0;
            int maxtruecount = -1;

            for (int truecount = 1; truecount <= 3; ++truecount){
                double mean_true = 0.0;
                for (int i = 0; i < truecount; ++i){
                    mean_true += 1.0/(double)truecount * -(double)cladetype_sort[i].first;
                }
                double ll = 0.0;
                for (int i = 0; i < truecount; ++i){
                    ll += dpois(-(double)cladetype_sort[i].first, mean_true);
                }
                if (truecount < 3){
                    double mean_err = 0.0;
                    for (int i = truecount; i < 3; ++i){
                        mean_err += 1.0/(double)(3.0-truecount) * -(double)cladetype_sort[i].first;
                    }
                    for (int i = truecount; i < 3; ++i){
                        ll += dpois(-(double)cladetype_sort[i].first, mean_err);
                    }
                }
                if (maxtruecount == -1 || ll > maxll){
                    maxll = ll;
                    maxtruecount = truecount;
                }
            }
            
            if (maxtruecount == 1 && cladetype_sort[0].second == 0){
                
                //fprintf(stdout, "Collapse %d -> %d\n", idx2, idx1);
                //fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\n", cladesize_A, cladesize_B, 
                //    cladesize_AminusB, cladesize_BminusA, cladesize_both, maxtruecount);     
                
                // There's one true clade, and it's the intersection of both.
                // Collapse these two clades into one and mark the sizes
                // of the erroneous clades for later.
                
                orig_to_collapsed.insert(make_pair(idx2, idx1));
                collapsed_to_orig[idx1].insert(idx2);
                
                errsize_tot_this += (double)cladesize_AminusB;
                errsize_tot_this += (double)cladesize_BminusA;
                errsize_count_this += 2.0;
            }

        }
    }
    fprintf(stderr, "%ld collapsed site groups remaining\n", collapsed_to_orig.size());
}

vector<hapstr> infer_clusters(hapstr& mask_global, 
    robin_hood::unordered_map<unsigned long, hap>& haplotypes, 
    int nvars,
    bool exact_matches_only,
    int num_expected,
    bool gex,
    float collapse_thresh,
    int& nclust_model,
    float& llr_model){
    
    map<int, robin_hood::unordered_set<unsigned long> > clades;
    map<int, robin_hood::unordered_set<unsigned long> > clades_not;
    map<int, robin_hood::unordered_set<unsigned long> > clades_mask;
   
    map<int, set<int> > collapsed_to_orig;
    map<int, int> orig_to_collapsed;
    
    vector<pair<long int, int> > clades_sort;
   
    get_clades(mask_global, haplotypes, nvars, clades, clades_not, clades_mask,
        clades_sort);

    collapse_sites(mask_global, haplotypes, nvars, clades, clades_not, clades_mask, 
        orig_to_collapsed, collapsed_to_orig, clades_sort);
    
    // Now that we have collapsed sets of identical sites, we need to dump the 
    // clades together and sort these new clades in decreasing order of frequency.

    clades_sort.clear();

    for (map<int, set<int> >::iterator co = collapsed_to_orig.begin(); co != collapsed_to_orig.end(); ++co){
        
        for (set<int>::iterator member = co->second.begin(); member != co->second.end(); ++member){
            for (robin_hood::unordered_set<unsigned long>::iterator cl = clades[*member].begin(); 
                cl != clades[*member].end(); ++cl){
                clades_mask[co->first].insert(*cl);
                if (clades_not[co->first].find(*cl) == clades_not[co->first].end()){
                    clades[co->first].insert(*cl);
                }
            }
            for (robin_hood::unordered_set<unsigned long>::iterator cl = clades_not[*member].begin();
                cl != clades_not[*member].end(); ++cl){
                clades_mask[co->first].insert(*cl);
                if (clades[co->first].find(*cl) == clades[co->first].end()){
                    clades_not[co->first].insert(*cl);
                }
            }
        }  
        clades_sort.push_back(make_pair(-clades[co->first].size(), co->first));
    }

    sort(clades_sort.begin(), clades_sort.end());
    
    // Store the cluster haplotypes learned by the model at the end 
    vector<hapstr> haps_final;
    
    pair<int, float> modelfit = infer_clusters2(mask_global, haplotypes, 
        nvars, clades_sort, 
        collapsed_to_orig, clades, clades_not,
        haps_final, exact_matches_only, num_expected);
    
    nclust_model = modelfit.first;
    llr_model = modelfit.second;

    mask_global.reset();
    for (int i = 0; i < haps_final.size(); ++i){
        mask_global |= haps_final[i];
    }
    return haps_final;
   
}

void parse_idsfile(string& idsfilename, vector<string>& ids){
    ifstream infile(idsfilename.c_str());
    string id;
    while (infile >> id ){
        ids.push_back(id);
    }
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"output_prefix", required_argument, 0, 'o'},
       {"barcodes", required_argument, 0, 'B'},
       {"mapq", required_argument, 0, 'q'},
       {"baseq", required_argument, 0, 'Q'},
       {"exact", no_argument, 0, 'e'},
       {"freq", required_argument, 0, 'f'},
       {"nclust", required_argument, 0, 'n'},
       {"cov", required_argument, 0, 'c'},
       {"mito", required_argument, 0, 'm'},
       {"mincount", required_argument, 0, 'C'},
       {"missing_thresh", required_argument, 0, 'M'},
       {"minsites", required_argument, 0, 's'},
       {"dump", no_argument, 0, 'd'},
       {"dump_raw", no_argument, 0, 'r'},
       {"vars", required_argument, 0, 'v'},
       {"nodoublet", no_argument, 0, 'N'},
       {"gex", no_argument, 0, 'g'},
       {"haps", required_argument, 0, 'H'},
       {"ids", required_argument, 0, 'i'},
       {"identity_thresh", required_argument, 0, 'I'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string output_prefix;
    int minmapq = 20;
    int minbaseq = 20;
    float minfreq = 0.001;
    int nclust = -1;
    bool dump = false;
    bool dump_raw = false;
    float missing_thresh = 0.25;
    string mito_chrom = "chrM";
    string varsfile;
    bool varsfile_given = false;
    string barcodesfilename;
    string hapsfilename;
    bool hapsfile_given = false;
    bool nodoublet = false;
    bool gex = false;
    bool exact_matches_only = false;
    bool ids_given = false;
    float collapse_thresh = 0;
    bool collapse_thresh_given = false;
    string idsfile;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:o:B:q:Q:f:n:m:v:H:i:I:egNrdh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
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
            case 'e':
                exact_matches_only = true;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'B':
                barcodesfilename = optarg;
                break;
            case 'q':
                minmapq = atoi(optarg);
                break;
            case 'Q':
                minbaseq = atoi(optarg);
                break;
            case 'f':
                minfreq = atof(optarg);
                break;
            case 'n':
                nclust = atoi(optarg);
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
            case 'r':
                dump_raw = true;
                break;
            case 'g':
                gex = true;
                break;
            case 'N':
                nodoublet = true;
                break;
            case 'I':
                collapse_thresh = atof(optarg);
                collapse_thresh_given = true;
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
    if (minfreq < 0 || minfreq > 0.5){
        fprintf(stderr, "ERROR: invalid value of --minfreq -f; should be 0 <= -f <= 0.5\n");
        exit(1);
    }
    if (hapsfile_given && !varsfile_given){
        fprintf(stderr, "ERROR: If --haps / -H is given, --vars / -v is also required\n");
        exit(1);
    }
    if (!collapse_thresh_given && gex){
        collapse_thresh = 0.07;
    }
    
    set<unsigned long> bc_whitelist;
    bool has_bc_whitelist = false;
    if (barcodesfilename.length() > 0){
        parse_barcode_file(barcodesfilename, bc_whitelist);
        if (bc_whitelist.size() > 0){
            has_bc_whitelist = true;
        }
    }

    // Get an ordered list of variable sites on the mitochondrial sequence
    deque<varsite> vars;
    map<int, varsite> vars_unfiltered;
    
    if (varsfile_given){
        load_vars_from_file(varsfile, vars); 
    }
    else{
        fprintf(stderr, "Finding variable sites on the mitochondrial genome...\n");
        find_vars_in_bam(bamfile, mito_chrom, minmapq, minbaseq, minfreq,
            vars, gex, has_bc_whitelist, bc_whitelist); 
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
        fprintf(stderr, "WARNING: removing %ld sites with lowest-frequency alt alleles\n", pos_rm.size());
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

        fprintf(stderr, "done\n");
    }
    
    // Now need to go back through the BAM file. This time, look at individual barcodes and count reads overlapping each variant.
    
    // Store barcodes as bit strings of length 2*nbases, interpreted as unsigned longs to save space
    //map<unsigned long, var_counts> hap_counter;
    robin_hood::unordered_map<unsigned long, var_counts> hap_counter;
        
    deque<varsite> vars2 = vars;

    // Get counts of each variant allele for each cell barcode
    
    // Necessary for filtering variants (if not provided),
    // inferring cluster haplotypes (if not provided),
    // and assigning barcodes to individual IDs
    
    count_vars_barcodes(bamfile, mito_chrom, minmapq, vars, 
        has_bc_whitelist, bc_whitelist, hap_counter, gex);     
    
    // Still need to filter variant sites based on coverage across cells
    set<int> vars_blacklist;
    //bitset<MAX_SITES> mask_global;
    hapstr mask_global;
    // Also translate variant counts into haplotype strings

    robin_hood::unordered_map<unsigned long, hap> haplotypes;
    
    if (dump_raw){
        write_vars(mito_chrom, output_prefix, vars2, mask_global, false);
    }
    
    fprintf(stderr, "turn var counts into cell haplotypes...\n");

    process_var_counts(hap_counter, haplotypes, vars_blacklist, varsfile_given,
        mask_global, nvars, gex, has_bc_whitelist);  
    
    if (dump_raw){
        return 0; // finished
    }

    // Store haplotypes of clusters (whether inferred or loaded)
    vector<hapstr> clusthaps;
    
    string statsfilename = output_prefix + ".summary";
    FILE* statsfilef = NULL;

    if (!dump){
        if (!hapsfile_given){
        
            
            int nclust_model;
            float llr_model;

            clusthaps = infer_clusters(mask_global, 
                haplotypes, nvars, exact_matches_only, 
                nclust, gex, collapse_thresh, 
                nclust_model, llr_model);
            
            if (statsfilef == NULL){
                statsfilef = fopen(statsfilename.c_str(), "w");
                
            }  
            fprintf(statsfilef, "%s\tnum_clusters\t%d\n", output_prefix.c_str(), nclust_model);
            fprintf(statsfilef, "%s\tllr_num_clusters\t%f\n", output_prefix.c_str(), llr_model);
            
            // Write haplotypes to output file
            string clusthapsfile = output_prefix + ".clusthaps";
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
        else{
            // Load clusters.
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
                //hap h;
                //h.vars.reset();
                //h.mask.reset();
                for (int i = 0; i < nvars; ++i){
                    if (hapstring[i] == '0'){
                        //h.mask.set(i);
                    }
                    else if (hapstring[i] == '1'){
                        //h.mask.set(i);
                        //h.vars.set(i);
                        h.set(i);
                    }
                }
                clusthaps.push_back(h);
            }
        }
    }
    
    // Avoid writing out vars file only if it would overwrite a preexisting vars file.
    if (!varsfile_given || (output_prefix + ".vars" != varsfile)){ 
        // Write variants file
         
        write_vars(mito_chrom, output_prefix, vars2, mask_global, true); 
    }
    
    if (true){
    //if (!hapsfile_given){
        // Write barcode haps file

        string haps_out = output_prefix + ".bchaps";
        FILE* haps_outf = fopen(haps_out.c_str(), "w");
        
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
        
        for (robin_hood::unordered_map<unsigned long, hap>::iterator h = 
            haplotypes.begin(); h != haplotypes.end(); ++h){
            
            bc bcbits(h->first);
            string bcstr = bc2str(bcbits, 16);

            //if ((float)(h->second.mask & mask_global).count() / (float)mask_global.count() >= 0.25){
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
        
        if (dump){
            return 0;
        }
    }
    
    // Assign all cell barcodes to a haplotype ID.
    map<unsigned long, int> bc2hap;
    
    map<int, map<int, float> > sites_doublet_probs;
    
    string assn_out = output_prefix + ".assignments";
    FILE* assn_outf = fopen(assn_out.c_str(), "w");

    robin_hood::unordered_map<unsigned long, int> assignments;
    robin_hood::unordered_map<unsigned long, double> assignments_llr;

    assign_bcs(hap_counter, assignments, assignments_llr, clusthaps,
        mask_global, nvars, nodoublet);
    
    map<string, int> id_counter;
    int tot_cells = 0;
    int doub_cells = 0;

    // Dump to output file
    char namebuf[100];
    for (robin_hood::unordered_map<unsigned long, int>::iterator assn = 
        assignments.begin(); assn != assignments.end(); ++assn){
        tot_cells++;
        string name;
        char s_d;
        if (clust_ids.size() > 0){
            name = idx2name(assn->second, clust_ids);
            if (assn->second >= clusthaps.size()){
                s_d = 'D';
                doub_cells++;
            }
            else{
                s_d = 'S';
            }
        }
        else{
            if (assn->second < clusthaps.size()){
                sprintf(&namebuf[0], "%d", assn->second);
                s_d = 'S';
            }
            else{
                pair<int, int> combo = idx_to_hap_comb(assn->second, 
                    clusthaps.size());
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
        fprintf(assn_outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(),
            name.c_str(), s_d, assignments_llr[assn->first]);
    } 
    fclose(assn_outf);
    
    if (statsfilef == NULL){
        statsfilef = fopen(statsfilename.c_str(), "w");
    }
    fprintf(statsfilef, "%s\ttot_cells\t%d\n", output_prefix.c_str(), tot_cells);
    if (!nodoublet){
        fprintf(statsfilef, "%s\tdoublets\t%d\n", output_prefix.c_str(), doub_cells);
        fprintf(statsfilef, "%s\tfrac_doublets\t%f\n", output_prefix.c_str(), (float)doub_cells/(float)tot_cells);
    }
    vector<pair<int, string> > idcounts_sorted;
    for (map<string, int>::iterator idc = id_counter.begin(); idc != id_counter.end(); ++idc){
        idcounts_sorted.push_back(make_pair(-idc->second, idc->first));
    }
    sort(idcounts_sorted.begin(), idcounts_sorted.end());
    for (vector<pair<int, string> >::iterator idcs = idcounts_sorted.begin(); idcs != idcounts_sorted.end(); ++idcs){
        fprintf(statsfilef, "%s\t%s\t%d\n", output_prefix.c_str(), idcs->second.c_str(), -idcs->first);   
    }

    if (statsfilef != NULL){
        fclose(statsfilef);
    }

    return 0;
}
