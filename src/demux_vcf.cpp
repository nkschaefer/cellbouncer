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
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <zlib.h>
#include <htswrapper/bc_hash.h>
#include <htswrapper/bam.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include <mixtureDist/functions.h>
#include "robin_hood.h"
#include "nnls.h"
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

// ===== Utility functions =====

/**
 * Gets the non-negative least squares solution
 * to the problem Ax = b given A and b.
 */
bool nn_lstsq(vector<vector<double> >& a,
    vector<double>& b,
    vector<double>& result_coefficients){
    
    int M = a.size(); // rows of A
    int N = a[0].size(); // cols of A
    int NRHS = 1; // cols of X
    int MDA = M; // min(1, M)
    int LDB = max(N,M); // rows of Y
    
    // Store results
    double* x = new double[N];
    
    // Allocate workspace variables    
    double* work = new double[N];
    double* zz = new double[M];
    int* index = new int[2*N];
    if (M != b.size()){
        fprintf(stderr, "dimensions do not match\n");
        exit(1);
    }
    
    // Populate A matrix
    // NOTE: entries are in order of columns, then rows
    double* a_lapack = new double[M*N];
    int k = 0;
    for (int j = 0; j < N; ++j){
        for (int i = 0; i < M; ++i){
            a_lapack[i+j*M] = a[i][j];
        }
    }
    // Populate B (y) matrix
    double* b_lapack = new double[LDB*NRHS];
    for (int i = 0; i < b.size(); ++i){
        b_lapack[i] = b[i];
    }
    int mode;
    double residual;
    int ret = nnls_c(a_lapack, &MDA, &M, &N, b_lapack, x, &residual,
        work, zz, index, &mode);
    if (mode == 3){
        fprintf(stderr, "Did not converge\n");
        return false;
    }
    else if (mode == 2){
        fprintf(stderr, "bad dimensions\n");
        exit(1);
    }
    if (result_coefficients.size() < a[0].size()){
        result_coefficients.clear();
        for (int i = 0; i < a[0].size(); ++i){
            result_coefficients.push_back(0.0);
        }
    }
    for (int i = 0; i < a[0].size(); ++i){
        result_coefficients[i] = x[i];
    }
    delete[] a_lapack;
    delete[] b_lapack; 
    delete[] x;
    delete[] work;
    delete[] zz;
    delete[] index;
    return true;
}
// ===== VCF-related functions =====

/**
 * Structure to represent population-level genotype information at a 
 * single SNP.
 */
struct var{
    char ref;
    char alt;
    
    // This limits us to 100 individuals in the input VCF
    bitset<100> haps1;
    bitset<100> haps2;
    bitset<100> haps_covered;
    vector<float> gqs;
    float vq;
    var(){
        this->ref = 'N';
        this->alt = 'N';
        this->haps1.reset();
        this->haps2.reset();
        this->haps_covered.reset();
        this->vq = 0.0;
    }
    var(const var& v){
        this->ref = v.ref;
        this->alt = v.alt;
        this->haps1 = v.haps1;
        this->haps2 = v.haps2;
        this->haps_covered = v.haps_covered;
        this->gqs = v.gqs;
        this->vq = v.vq;
    }
};

/**
 * Read variant data from VCF.
 */
int read_vcf(string& filename, 
    bam_reader& reader,
    vector<string>& samples,
    map<int, map<int, var> >& snps,
    int min_vq,
    bool hdr_only){
    
    // Map sequence names to TIDs for storage
    map<string, int> seq2tid;
    if (!hdr_only){
        seq2tid = reader.get_seq2tid();
    }

    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(filename.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", filename.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
    }
    if (hdr_only){
        return -1;
    }
    long int nvar = 0;
    
    float mingq = 30;
    
    // blacklist for duplicate variants
    set<pair<int, int> > bl;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        // Convert to BAM-compatible TID
        int tid = seq2tid[chrom];
        //tid = bcf_record->tid;

        // Get 0-based coordinate. Coordinates in BAM will also be 0-based. 
        int pos = bcf_record->pos;
        
        pair<int, int> key = make_pair(tid, pos);
        if (bl.find(key) != bl.end()){
            // If we've already seen a site twice, ignore it.
            continue;
        }

        if (snps.count(tid) > 0 && snps[tid].count(pos) > 0){
            // If we've already seen this site once, delete the previous entry 
            // for it and note that we should not store any additional variants
            // at this site.
            fprintf(stderr, "WARNING: duplicate variants at site %s:%d\n", chrom.c_str(), pos+1);
            snps[tid].erase(pos);
            bl.insert(key);
        }

        if (bcf_record->n_allele == 2){ 
            
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

            // pass = biallelic, no indels, ref/alt both A/C/G/T
            bool pass = true;
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
            else if (bcf_record->qual < min_vq){
                pass = false;
            }
            if (pass){
                if (snps.count(tid) == 0){
                    map<int, var> m;
                    snps.insert(make_pair(tid, m));
                }
                var v;
                v.ref = bcf_record->d.allele[0][0];
                v.alt = bcf_record->d.allele[1][0];
                v.vq = bcf_record->qual;
                ++nvar;

                // Get all available genotypes.
                int32_t* gts = NULL;
                int n_gts = 0;
                int nmiss = 0;
                int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                if (num_loaded <= 0){
                    fprintf(stderr, "ERROR loading genotypes at %s %lld\n", 
                        chrom.c_str(), bcf_record->pos);
                    exit(1);
                }
                
                // Assume ploidy = 2
                int ploidy = 2;
                //int ploidy = n_gts / num_samples; 
                
                // Load genotype qualities
                float* gqs = NULL;
                int n_gqs = 0;
                int num_gq_loaded = bcf_get_format_float(bcf_header, bcf_record, "GQ",
                    &gqs, &n_gqs);
                
                int nalt_alleles = 0;
                int nalt_samples = 0;

                for (int i = 0; i < num_samples; ++i){
                    int32_t* gtptr = gts + i*ploidy;
                    
                    bool gq_pass = false;
                    
                    if (num_gq_loaded < num_samples || !isnan(gqs[i]) || 
                        gqs[i] == bcf_float_missing){
                        // Missing GQ? let it slide
                        gq_pass = true;
                    }
                    else{
                        // valid GQ.
                        if (gqs[i] >= mingq){
                            v.gqs.push_back(pow(10, -(float)gqs[i] / 10.0));
                        }
                        else{
                            gq_pass = false;
                            v.gqs.push_back(-1);
                        }
                    }
                    
                    if (bcf_gt_is_missing(gtptr[0])){
                        // Missing genotype.
                        nmiss++;
                    }
                    else if (!gq_pass){
                        // Set to missing
                        nmiss++;
                    }
                    else{
                        bool alt = false;   
                        v.haps_covered.set(i);
                        if (bcf_gt_allele(gtptr[0]) == 1){
                            nalt_alleles++;
                            alt = true;
                            v.haps1.set(i);
                        }
                        if (bcf_gt_allele(gtptr[1]) == 1){
                            nalt_alleles++;
                            alt = true;
                            v.haps2.set(i);
                        }
                        if (alt){
                            nalt_samples++;
                        }
                    }    
                } 
                free(gqs);            
                snps[tid].insert(make_pair(pos, v));
            }
        }
    }
    fprintf(stderr, "Loaded %ld variants\n", nvar);
    
    bcf_hdr_destroy(bcf_header);
    bcf_destroy(bcf_record);
    bcf_close(bcf_reader);
    
    return nvar;
}

// ===== BAM-related functions =====

/**
 * Extract an allele count from the current BAM record.
 */
void process_bam_record(bam_reader& reader,
    int snppos,
    var& vardat,
    map<int, robin_hood::unordered_map<unsigned long, pair<float, float> > >& var_counts,
    bool has_bc_list,
    set<unsigned long>& bcs_valid,
    bool neg_drops,
    unsigned long negdrop_cell){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
                        
        // Get BC key
        bc bc_bits;
        str2bc(reader.cb_z, bc_bits, 16);
        unsigned long bc_key = bc_bits.to_ulong();
        
        bool neg_pass = false;
        if (neg_drops && bcs_valid.find(bc_key) == bcs_valid.end()){
            neg_pass = true;
            bc_key = negdrop_cell;
        }
        
        if (!has_bc_list || neg_pass || bcs_valid.find(bc_key) != bcs_valid.end()){
            
            // Instead of storing actual read counts, store the probability
            // that the mapping was correct.
            float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);

            if (var_counts.count(snppos) == 0){
                robin_hood::unordered_map<unsigned long, pair<float, float> > m;
                var_counts.insert(make_pair(snppos, m));
            }
            if (var_counts[snppos].count(bc_key) == 0){
                pair<float, float> p = make_pair(0,0);
                var_counts[snppos].emplace(bc_key, p);
            }

            // Note: this function expects positions to be 1-based, but 
            // BCF/BAM functions store as 0-based
            char allele = reader.get_base_at(snppos + 1);
            
            if (allele != 'N' && allele != '-'){
                if (allele == vardat.ref){
                    var_counts[snppos][bc_key].first += prob_corr;
                }
                else if (allele == vardat.alt){
                    var_counts[snppos][bc_key].second += prob_corr;
                }
            }
        }
    }
}

/**
 * Given a set of allele counts at a given site, once we are guaranteed no longer to 
 * see the site in the BAM, we can dump the information from the site-specific
 * data structure into the genome-wide data structure storing counts at different
 * types of alleles per cell.
 */
void dump_vcs_counts(robin_hood::unordered_map<unsigned long, pair<float, float> >& varcounts_site,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    var& snpdat,
    int n_samples){
    
    static int snpid = 0;
    
    vector<int> sample_idx;
    for (int i = 0; i < n_samples; ++i){
        sample_idx.push_back(i);
    }

    for (robin_hood::unordered_map<unsigned long, pair<float, float> >::iterator vcs = 
        varcounts_site.begin(); vcs != varcounts_site.end(); ++vcs){
        
        // Ensure a map exists for the current cell barcode
        if (indv_allelecounts.count(vcs->first) == 0){
            
            map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
            indv_allelecounts.emplace(vcs->first, m);
            
            // Ensure a map exists for the current cell barcode for each
            // individual + site type combination
            for (int i = 0; i < n_samples; ++i){
                map<pair<int, int>, pair<float, float> > m2;
                for (int j = 0; j < 3; ++j){
                    pair<int, int> key = make_pair(i, j);
                    indv_allelecounts[vcs->first].insert(make_pair(key, m2));
                }
            }
        } 
        
        // Only store information for SNPs with non-zero allele counts
        if (vcs->second.first + vcs->second.second > 0){
            for (int i = 0; i < n_samples; ++i){
                int n_alt_chroms = 0;
                if (snpdat.haps_covered.test(i)){
                    if (snpdat.haps1.test(i)){
                        n_alt_chroms++;
                    }
                    if (snpdat.haps2.test(i)){
                        n_alt_chroms++;
                    }
                    
                    pair<int, int> key = make_pair(i, n_alt_chroms);
                    
                    // Store total keyed to -1
                    pair<int, int> nullkey = make_pair(-1, -1);
                    if (indv_allelecounts[vcs->first][key].count(nullkey) == 0){
                        indv_allelecounts[vcs->first][key].insert(make_pair(nullkey, 
                            make_pair(0.0,0.0)));    
                    }
                    indv_allelecounts[vcs->first][key][nullkey].first += vcs->second.first;
                    indv_allelecounts[vcs->first][key][nullkey].second += vcs->second.second; 
                    
                    // Check this site's allelic state in other individuals
                    for (int j = i + 1; j < n_samples; ++j){
                        if (snpdat.haps_covered.test(j)){
                            int n_alt_chroms_j = 0;
                            if (snpdat.haps1.test(j)){
                                n_alt_chroms_j++;
                            }
                            if (snpdat.haps2.test(j)){
                                n_alt_chroms_j++;
                            }
                            pair<int, int> key_j = make_pair(j, n_alt_chroms_j);
                            
                            // This should never happen -- but ensures individual i
                            // is sort order < individual j 
                            if (key.first > key_j.first){
                                pair<int, int> tmp = key;
                                key = key_j;
                                key_j = tmp;
                            }

                            if (indv_allelecounts[vcs->first][key].count(key_j) == 0){
                                indv_allelecounts[vcs->first][key].insert(make_pair(key_j, 
                                    make_pair(0.0,0.0)));
                            }                       
                            indv_allelecounts[vcs->first][key][key_j].first += 
                                vcs->second.first;
                            indv_allelecounts[vcs->first][key][key_j].second += 
                                vcs->second.second;
                        }
                    }       
                }
            }
        }
    }
    ++snpid; 
}

// ===== Other file I/O =====

/**
 * Read file of allowed ID assignments
 */
void parse_idfile(string& idfile, 
    vector<string>& samples,
    set<int>& ids_allowed){

    map<string, int> name2idx;
    for (int i = 0; i < samples.size(); ++i){
        name2idx.insert(make_pair(samples[i], i));
    }
    ifstream infile(idfile.c_str());
    string id;
    while (infile >> id){
        size_t sep_pos = id.find("x");
        if (sep_pos == string::npos){
            sep_pos = id.find("X");
        }
        if (sep_pos == string::npos){
            sep_pos = id.find("+");
        }
        if (sep_pos != string::npos){
            string id1 = id.substr(0, sep_pos);
            string id2 = id.substr(sep_pos+1, id.length()-sep_pos-1);
            int idx1 = -1;
            int idx2 = -1;
            if (name2idx.count(id1) == 0){
                fprintf(stderr, "WARNING: indv %s from individual file not found in VCF\n", 
                    id1.c_str());
            }
            else{
                idx1 = name2idx[id1];
            }
            if (name2idx.count(id2) == 0){
                fprintf(stderr, "WARNING: indv %s from individual file not found in VCF\n", 
                    id2.c_str());
            }
            else{
                idx2 = name2idx[id2];
            }
            if (idx1 != -1 && idx2 != -1){
                int k;
                if (idx1 < idx2){
                    k = hap_comb_to_idx(idx1, idx2, samples.size());
                }
                else{
                    k = hap_comb_to_idx(idx2, idx1, samples.size());
                }
                ids_allowed.insert(k);
            }
        }
        else{
            if (name2idx.count(id) == 0){
                fprintf(stderr, "WARNING: indv %s from individual file not found in VCF\n", 
                    id.c_str());
            }
            else{
                ids_allowed.insert(name2idx[id]);
            }
        }
    }

    // One potential issue: if we allow a combination, we must allow the singlet versions of
    // both halves (both for data structures to work and logically - a doublet is the result
    // of two single cells, and a tetraploid fusion may have not worked, leaving its component
    // individuals in the pool as well)
    set<int> adds;
    for (set<int>::iterator id = ids_allowed.begin(); id != ids_allowed.end(); ++id){
        if (*id >= samples.size()){
            pair<int, int> comb = idx_to_hap_comb(*id, samples.size());
            if (ids_allowed.find(comb.first) == ids_allowed.end()){
                adds.insert(comb.first);
            }
            if (ids_allowed.find(comb.second) == ids_allowed.end()){
                adds.insert(comb.second);
            }
        }
    }
    for (set<int>::iterator a = adds.begin(); a != adds.end(); ++a){
        ids_allowed.insert(*a);
    }
}

/**
 * If a previous run was dumped to count files, load those counts instead of 
 * re-processing the BAM file.
 */
void load_counts_from_file(
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    map<int, map<int, robin_hood::unordered_map<unsigned long, pair<float, float> > > >& varcounts_site,
    vector<string>& indvs,   
    string& filename){
    
    // Use sample -> index mapping from VCF header 
    map<string, int> indv2idx;
    for (int i = 0; i < indvs.size(); ++i){
        indv2idx.insert(make_pair(indvs[i], i));
    }
    
    ifstream infile(filename.c_str());
    unsigned long cell;
    string indv;
    int type;
    string indv2;
    int type2;
    float count1;
    float count2;

    while (infile >> cell >> indv >> type >> indv2 >> type2 >> count1 >> count2){
        size_t sep_pos = indv.find("+");
        if (sep_pos == string::npos){
            sep_pos = indv.find("x");
        }
        int indv_idx = -1;
        if (sep_pos != string::npos){
            string indv1 = indv.substr(0, sep_pos);
            string indv2 = indv.substr(sep_pos+1, indv.length()-sep_pos-1);
            
            if (indv2idx.count(indv1) == 0){
                fprintf(stderr, "ERROR: could not find %s in variant data\n", indv1.c_str());
                exit(1);
            }
            if (indv2idx.count(indv2) == 0){
                fprintf(stderr, "ERROR: could not find %s in variant data\n", indv2.c_str());
                exit(1);
            }
            int idx1 = indv2idx[indv1];
            int idx2 = indv2idx[indv2];
        
            indv_idx = hap_comb_to_idx(idx1, idx2, indvs.size());
        }
        else{
            if (indv2idx.count(indv) == 0){
                fprintf(stderr, "ERROR: count not find %s in variant data\n", indv.c_str());
                exit(1);
            }
            indv_idx = indv2idx[indv];
        }
        int indv_idx2 = -1;
        if (indv2 != "NA"){
            
            size_t sep_pos2;
            sep_pos2 = indv2.find("+");
            if (sep_pos2 == string::npos){
                sep_pos2 = indv2.find("x");
            }   
            if (sep_pos2 != string::npos){
                string indv2a = indv2.substr(0, sep_pos2);
                string indv2b = indv2.substr(sep_pos2+1, indv2.length()-1-sep_pos2); 
                if (indv2idx.count(indv2a) == 0){
                    fprintf(stderr, "ERROR: could not find %s in variant data\n", indv2a.c_str());
                    exit(1);
                }
                if (indv2idx.count(indv2b) == 0){
                    fprintf(stderr, "ERROR: could not find %s in variant data\n", indv2b.c_str());
                    exit(1);
                }
                
                int idx2a = indv2idx[indv2a];
                int idx2b = indv2idx[indv2b];
        
                indv_idx2 = hap_comb_to_idx(idx2a, idx2b, indvs.size());
            }
            else{
                if (indv2idx.count(indv2) == 0){
                    fprintf(stderr, "ERROR: could not find %s in variant data\n", indv2.c_str());
                    exit(1);
                }
                indv_idx2 = indv2idx[indv2];
            }
        }
        if (indv_allelecounts.count(cell) == 0){
            map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
            indv_allelecounts.emplace(cell, m);
        }
        
        pair<int, int> key = make_pair(indv_idx, type);
        if (indv_allelecounts[cell].count(key) == 0){
            map<pair<int, int>, pair<float, float> > m;
            indv_allelecounts[cell].insert(make_pair(key, m));
        }
        pair<int, int> key2 = make_pair(indv_idx2, type2);
        if (indv_idx2 < indv_idx && indv_idx2 != -1){
            pair<int, int> tmp = key;
            key = key2;
            key2 = tmp;
        }
        indv_allelecounts[cell][key].insert(make_pair(key2, make_pair(count1, count2)));
    }
}

/**
 * Print counts to text files.
 */
void dump_cellcounts(FILE* out_cell,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts, 
    vector<string>& samples){
    
    for (robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >::iterator x = indv_allelecounts.begin();
        x != indv_allelecounts.end(); ++x){
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
                x->second.begin(); y != x->second.end(); ++y){
            string indv = idx2name(y->first.first, samples);
            int nalt = y->first.second;
            for (map<pair<int, int>, pair<float, float> >::iterator z = 
                y->second.begin(); z != y->second.end(); ++z){
                string indv2 = "NA";
                if (z->first.first != -1){
                    indv2 = idx2name(z->first.first, samples);
                }
                int nalt2 = z->first.second;
                fprintf(out_cell, "%ld\t%s\t%d\t%s\t%d\t%f\t%f\n", x->first, indv.c_str(),
                    nalt, indv2.c_str(), nalt2, z->second.first, z->second.second);
            }
        }
    } 
}



// ===== Functions related to deciding identities of cells =====



/** 
 * Helper function used by compute_kcomps().
 * 
 * Computes comparisons between two doublet combinations conditional on
 * a shared individual. 
 *
 * In other words, given that individual A is involved, can compute the
 * LLR of (A+B) vs (A+C), by considering the LLR of (A+B)/A vs (A+C)/A.
 */
void get_kcomps_cond(map<int, map<int, double> >& kcomps,
    map<int, map<int, double> >& llrs,
    int n_samples,
    vector<string>& samples,
    set<int>& allowed_assignments){

    for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end();
        ++llr){
        vector<int> js;
        vector<int> comps;
        for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != llr->second.end();
            ++llr2){
            if (llr2->first >= n_samples){
                pair<int, int> comb = idx_to_hap_comb(llr2->first, n_samples);
                int j = -1;
                if (comb.first == llr->first){
                    j = comb.second;
                }
                else{
                    j = comb.first;
                }
                js.push_back(j);
                comps.push_back(-llr2->second);
            }
        }
        
        if (comps.size() > 1){
            for (int i = 0; i < comps.size()-1; ++i){
                int idx1 = js[i];
                for (int j = i +1; j < comps.size(); ++j){
                    int idx2 = js[j];
                    double llr = comps[i] - comps[j];
                    if (idx1 < idx2){
                        if (kcomps.count(idx1) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx1, m));
                        }
                        if (kcomps[idx1].count(idx2) == 0){
                            kcomps[idx1].insert(make_pair(idx2, 0.0));
                        }
                        kcomps[idx1][idx2] += llr;
                    }
                    else{
                        if (kcomps.count(idx2) == 0){
                            map<int, double> m;
                            kcomps.insert(make_pair(idx2, m));
                        }
                        if (kcomps[idx2].count(idx1) == 0){
                            kcomps[idx2].insert(make_pair(idx1, 0.0));
                        }
                        kcomps[idx2][idx1] += llr;
                    }
                }
            }
        }
    }
}

/**
 * Used by populate_llr_table(). 
 *
 * Whereas that function can directly
 * compare the likelihood of a doublet combination with one of its
 * singlet components (i.e. individual A+B vs individual A and 
 * individual A+B vs individual B), this function computes all other
 * possible comparisons involving doublet combinations (i.e.
 * individual A+B vs individual C, individual A+B vs individual A+D,
 * and individual A+B vs individual C+D.
 * 
 * Puts these comparisons into the same log likelihood ratio table
 * so a single decision can be made across all comparisons.
 */
void compute_k_comps(map<int, map<int, double> >& llrs,
    vector<int>& ks,
    vector<string>& samples,
    set<int>& allowed_assignments){
    
    // Get a table of LLRs between different combination members, conditional on a cell
    // being half another member. In other words, keys/values here are
    // ID1 -> ID2A -> ID2B -> LLR
    // LLR is the log likelihood ratio of ID2A vs ID2B being included in a double model
    // with ID1. 
    map<int, map<int, map<int, double> > > kcomp_cond;
    
    for (map<int, map<int, double> >::iterator samp = llrs.begin(); samp != llrs.end(); ++samp){
        map<int, map<int, double> > llrstmp;
        llrstmp.insert(make_pair(samp->first, samp->second));
        map<int, map<int, double> > kcomp; 
        get_kcomps_cond(kcomp, llrstmp, samples.size(), samples, allowed_assignments);
        kcomp_cond.insert(make_pair(samp->first, kcomp));
    }
    
    // compute all single vs double comparisons
    for (int i = 0; i < samples.size(); ++i){
        if (allowed_assignments.size() != 0 && 
            allowed_assignments.find(i) == allowed_assignments.end()){
            continue;
        }
        for (int ki = 0; ki < ks.size(); ++ki){
            int k = ks[ki];
            pair<int, int> comb = idx_to_hap_comb(k, samples.size());
            if (comb.first != i && comb.second != i){
                double ll1 = llrs[comb.first][k];
                double ll2 = llrs[comb.second][k];
                if (comb.first < i){
                    ll1 += -llrs[comb.first][i];
                }
                else{
                    ll1 += llrs[i][comb.first];
                }
                if (comb.second < i){
                    ll2 += -llrs[comb.second][i];
                }
                else{
                    ll2 += llrs[i][comb.second];
                }
                llrs[i].insert(make_pair(k, 0.5*ll1 + 0.5*ll2));
            }
        }
    }
    // compute all double vs double comparisons
    for (int ki = 0; ki < ks.size()-1; ++ki){
        int k1 = ks[ki];
        map<int, double> m;
        llrs.insert(make_pair(k1, m));
        pair<int, int> comb1 = idx_to_hap_comb(k1, samples.size());
        for (int kj = ki+1; kj < ks.size(); ++kj){
            int k2 = ks[kj];
            pair<int, int> comb2 = idx_to_hap_comb(k2, samples.size());
            vector<double> llr_parts;
            int a = comb1.first;
            int b = comb1.second;
            int c = comb2.first;
            int d = comb2.second;

            if (a == c){
                if (b < d){
                    llr_parts.push_back(kcomp_cond[a][b][d]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][d][b]);
                }
            }
            else if (a == d){
                if (b < c){
                    llr_parts.push_back(kcomp_cond[a][b][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[a][c][b]);
                }
            }
            else if (b == c){
                if (a < d){
                    llr_parts.push_back(kcomp_cond[b][a][d]);
                }         
                else{
                    llr_parts.push_back(-kcomp_cond[b][d][a]);
                }
            }
            else if (b == d){
                if (a < c){
                    llr_parts.push_back(kcomp_cond[b][a][c]);
                }
                else{
                    llr_parts.push_back(-kcomp_cond[b][c][a]);
                }
            }
            else{
                double llr1 = 0;
                double llr2 = 0;
                double llr3 = 0;
                double llr4 = 0;
                double llr5 = 0;
                double llr6 = 0;
                double llr7 = 0;
                double llr8 = 0;
                // Calculate (A+B)/(C+D) 
                if (b < c){
                    llr1 += kcomp_cond[a][b][c];
                }
                else{
                    llr1 += -kcomp_cond[a][c][b];
                }
                if (a < d){
                    llr1 += kcomp_cond[c][a][d];
                }
                else{
                    llr1 += -kcomp_cond[c][d][a];
                }
                
                if (a < c){
                    llr2 += kcomp_cond[b][a][c];
                }
                else{
                    llr2 += -kcomp_cond[b][c][a];
                }
                if (b < d){
                    llr2 += kcomp_cond[c][b][d];
                }
                else{
                    llr2 += -kcomp_cond[c][d][b];
                }

                if (b < d){
                    llr3 += kcomp_cond[a][b][d];
                }
                else{
                    llr3 += -kcomp_cond[a][d][b];
                }
                if (a < c){
                    llr3 += kcomp_cond[d][a][c];
                }
                else{
                    llr3 += -kcomp_cond[d][c][a];
                }

                if (a < d){
                    llr4 += kcomp_cond[b][a][d];
                }
                else{
                    llr4 += -kcomp_cond[b][d][a];
                }
                if (b < c){
                    llr4 += kcomp_cond[d][b][c];
                }
                else{
                    llr4 += -kcomp_cond[d][c][b];
                }

                if (d < a){
                    llr5 += kcomp_cond[c][d][a];
                }
                else{
                    llr5 += -kcomp_cond[c][a][d];
                }
                if (c < b){
                    llr5 += kcomp_cond[a][c][b];
                }
                else{
                    llr5 += -kcomp_cond[a][b][c];
                }
                llr5 = -llr5;

                if (c < a){
                    llr6 += kcomp_cond[d][c][a];
                }
                else{
                    llr6 += -kcomp_cond[d][a][c];
                }
                if (d < b){
                    llr6 += kcomp_cond[a][d][b];
                }
                else{
                    llr6 += -kcomp_cond[a][b][d];
                }
                llr6 = -llr6;

                if (d < b){
                    llr7 += kcomp_cond[c][d][b];
                }
                else{
                    llr7 += -kcomp_cond[c][b][d];
                }
                if (c < a){
                    llr7 += kcomp_cond[b][c][a];
                }
                else{
                    llr7 += -kcomp_cond[b][a][c];
                }
                llr7 = -llr7;

                if (c < b){
                    llr8 += kcomp_cond[d][c][b];
                }
                else{
                    llr8 += -kcomp_cond[d][b][c];
                }
                if (d < a){
                    llr8 += kcomp_cond[b][d][a];
                }
                else{
                    llr8 += -kcomp_cond[b][a][d];
                }
                llr8 = -llr8;
                llr_parts.push_back(llr1);
                llr_parts.push_back(llr2);
                llr_parts.push_back(llr3);
                llr_parts.push_back(llr4);
                llr_parts.push_back(llr5);
                llr_parts.push_back(llr6);
                llr_parts.push_back(llr7);
                llr_parts.push_back(llr8);
            }
            
            float llrsum = 0.0;
            float llrcount = 0;
            for (int idx = 0; idx < llr_parts.size(); ++idx){
                llrsum += llr_parts[idx];
                llrcount++;
            }
            double llr = llrsum / llrcount;
            llrs[k1].insert(make_pair(k2, llr));
        }
    }
    
    // Finally, if we are filtering assignments, need to delete all disallowed identities from
    // the LLR table, now that we used them to help compute LLRs of allowed identities.
    if (allowed_assignments.size() > 0){
        for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end(); ){
            if (allowed_assignments.find(llr->first) != 
                allowed_assignments.end()){
                for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != llr->second.end(); ){
                    if (allowed_assignments.find(llr2->first) != 
                        allowed_assignments.end()){
                        ++llr2;
                    }
                    else{
                        llr->second.erase(llr2++);
                    }
                }
                ++llr;
            }
            else{
                llrs.erase(llr++);
            }
        }
    }
}

/**
 * Given a set of allele counts (at all possible SNP types) for a single cell,
 * populates a log likelihood ratio table for that cell, which gives the LLR
 * of every possible identity vs every other possible identity.
 */
void populate_llr_table(map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > >& counts,
    map<int, map<int, double> >& llrs,
    vector<string>& samples,
    set<int>& allowed_assignments,
    double doublet_rate){
    
    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
        counts.begin(); y != counts.end(); ++y){
        
        /* 
        if (allowed_assignments.size() > 0 && allowed_assignments.find(y->first.first) == 
            allowed_assignments.end()){
            continue;
        }
        */

        // Set default expectation for indv1
        // 0 = homozygous ref (~0% alt allele)
        // 1 = heterozygous (~50% alt allele)
        // 2 = homozygous alt (~100% alt allele)
        float exp1 = 0.001;
        if (y->first.second == 1){
            exp1 = 0.5;
        }
        else if (y->first.second == 2){
            exp1 = 0.999;
        }

        for (map<pair<int, int>, pair<float, float> >::iterator z = 
            y->second.begin(); z != y->second.end(); ++z){
            
            /*    
            if (allowed_assignments.size() > 0 && allowed_assignments.find(z->first.first) ==
                allowed_assignments.end()){
                continue;
            } 
            */          

            // If same site type, we can't distinguish between the
            // two individuals from this piece of information
            if (z->first.first != -1 && y->first.second != z->first.second){
                
                // Set default expectation for indv2
                float exp2 = 0.001;
                if (z->first.second == 1){
                    exp2 = 0.5;
                }
                else if (z->first.second == 2){
                    exp2 = 0.999;
                }
                
                float exp3 = 0.5*exp1 + 0.5*exp2;
               
                int i = y->first.first;
                int j = z->first.first;
                int k = hap_comb_to_idx(i, j, samples.size());

                int ref = (int)round(z->second.first);
                int alt = (int)round(z->second.second);
                
                double ll1 = dbinom(ref+alt, alt, exp1);   
                double ll2 = dbinom(ref+alt, alt, exp2);
                double ll3 = dbinom(ref+alt, alt, exp3);
                
                map<int, double> m;
                if (llrs.count(i) == 0){
                    llrs.insert(make_pair(i, m));
                }
                if (llrs.count(j) == 0){
                    llrs.insert(make_pair(j, m));
                }
                if (llrs[i].count(j) == 0){
                    llrs[i].insert(make_pair(j, 0.0));
                }
                llrs[i][j] += (ll1-ll2);
                
                if (doublet_rate > 0.0){
                    // Store comparisons between i and (i,j) combo and 
                    // between j and (i,j) combo
                    if (llrs[i].count(k) == 0){
                        llrs[i].insert(make_pair(k, 0.0));
                    }
                    if (llrs[j].count(k) == 0){
                        llrs[j].insert(make_pair(k, 0.0));
                    }
                    
                    llrs[i][k] += (ll1-ll3);
                    llrs[j][k] += (ll2-ll3);
                }
            }
        }
    }
    
    if (doublet_rate > 0.0){ 
        
        // Get a list of all possible double identities to consider.
        vector<int> ks;
        for (int i = 0; i < samples.size()-1; ++i){
            /*
            if (allowed_assignments.size() != 0 && 
                allowed_assignments.find(i) == allowed_assignments.end()){
                continue;
            }
            */
            for (int j = i + 1; j < samples.size(); ++j){
                /*
                if (allowed_assignments.size() != 0 && 
                    allowed_assignments.find(j) == allowed_assignments.end()){
                    continue;
                }
                */
                int k = hap_comb_to_idx(i, j, samples.size());
                ks.push_back(k);
            }
        }

        // put k indices in increasing order so we know what order to store comparisons
        // in the data structure
        sort(ks.begin(), ks.end());

        // With all possible values of k to consider, we already have all member component vs 
        // double model comparisons computed. We now need to compute all other possible
        // single vs double model comparisons, as well as double model vs double model comparisons.
        
        compute_k_comps(llrs, ks, samples, allowed_assignments);
    }

    // Erase impossible combinations.
    //
    // Only need to check if we've limited allowable IDs or set doublet rate to 1
    // (doublet rate set to zero already would have excluded doublet combos from
    // being added to the LLR table)
    
    if (allowed_assignments.size() > 0 || doublet_rate == 1.0){
        for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end(); ){
            bool pass = true;
            if (doublet_rate == 1 && llr->first < samples.size()){
                pass = false;
            }
            else if (allowed_assignments.size() > 0 && 
                allowed_assignments.find(llr->first) == allowed_assignments.end()){
                pass = false;
            }
            if (pass){
                for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != llr->second.end(); ){
                    
                    bool pass2 = true;
                    if (doublet_rate == 1 && llr2->first < samples.size()){
                        pass2 = false;
                    }
                    else if (allowed_assignments.size() > 0 && 
                        allowed_assignments.find(llr2->first) == allowed_assignments.end()){
                        pass2 = false;
                    }
                    if (pass2){    
                        ++llr2;
                    }
                    else{
                        llr->second.erase(llr2++);
                    }
                }
                ++llr;
            }
            else{
                llrs.erase(llr++);
            }
        }
    }
    if (doublet_rate > 0.0 && doublet_rate < 1.0){
        // Finally, need to adjust LLR table using doublet prior.
        for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end();
            ++llr){
            // Only need to adjust ratios where one is a single and the
            // other is a double. 
            // Since sorted numerically, the first index will never be 
            // a double in a single/double combination.
            if (llr->first >= samples.size()){
                break;
            }
            for (map<int, double>::iterator llr2 = llr->second.begin(); llr2 != 
                llr->second.end(); ++llr2){
                if (llr->first < samples.size() && llr2->first >= samples.size()){
                    // Adjust LLR.
                    llr2->second += log2(1.0-doublet_rate) - log2(doublet_rate);
                }  
            }
        }
    }
}

/**
 * Given genome-wide counts of different types of alleles, determines the 
 * most likely identity of each cell and stores in the data structure.
 */
void assign_ids(robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    vector<string>& samples,
    robin_hood::unordered_map<unsigned long, int>& assignments,
    robin_hood::unordered_map<unsigned long, double>& assignments_llr,
    set<int>& allowed_assignments,
    bool negdrops,
    unsigned long negdrop_cell,
    double doublet_rate){
    
    // Left in for debugging purposes only: if print_llrs is true,
    // the user will be prompted for a cell barcode of interest.
    // Once found, the log likelihood ratio table for the selected
    // cell barcode will be printed to stdout and the program will
    // exit.
    bool print_llrs = false;
    unsigned long searchbc;
    if (print_llrs){
        string searchbc_str;
        cerr << "Enter barcode: ";
        cin >> searchbc_str;
        if (searchbc_str.length() != 16){
            cerr << "ERROR: input barcode length is not 16";
            exit(1);
        }
        cerr << "Searching for " << searchbc_str;
        cerr << "\n"; 
        bc searchbc_bin;
        str2bc(searchbc_str.c_str(), searchbc_bin, 16);
        searchbc = searchbc_bin.to_ulong();
    }
    
    for (robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >::iterator x = 
            indv_allelecounts.begin(); x != indv_allelecounts.end(); ++x){
        
        // If we were able to compile negative droplets and the current
        // cell is the sum of all negative droplets, skip trying to
        // identify it.
        if (negdrops && x->first == negdrop_cell){
            continue;
        }
        if (print_llrs && x->first != searchbc){
            continue;
        }
        
        // Get a table of log likelihood ratios between every possible
        // pair of identities
        map<int, map<int, double> > llrs;
        populate_llr_table(x->second, llrs, samples, allowed_assignments, doublet_rate);

        // Debugging only: print this table and quit
        if (print_llrs){
            bc as_bitset(x->first);
            string bc_str = bc2str(as_bitset, 16);
            for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != 
                llrs.end(); ++llr){
                for (map<int, double>::iterator llr2 = llr->second.begin();
                    llr2 != llr->second.end(); ++llr2){
                    fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                        idx2name(llr->first, samples).c_str(),
                        idx2name(llr2->first, samples).c_str(),
                        llr2->second);
                    fprintf(stdout, "%s\t%s\t%s\t%f\n", bc_str.c_str(),
                        idx2name(llr2->first, samples).c_str(),
                        idx2name(llr->first, samples).c_str(),
                        -llr2->second);
                }
            }
            exit(0);
        }    
        
        // Now that we have a table of all pairwise LLRs between assignments for this cell,
        // choose the best identity.
        
        // This is done by iteratively finding the largest LLR between identities and eliminating
        // the least likely one, until only two identities are left. The final LLR is the LLR between
        // the best and second best assignment.
        
        double llr_final;
        int assn = collapse_llrs(llrs, llr_final);
        assignments.emplace(x->first, assn);
        assignments_llr.emplace(x->first, llr_final);
    }            
}

/**
 * Spill cell -> individual assignments to disk
 */
void dump_assignments(FILE* outf,
    robin_hood::unordered_map<unsigned long, int>& assn_final,
    robin_hood::unordered_map<unsigned long, double>& assn_final_llr,
    vector<string>& samples){

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn_final.begin();
        a != assn_final.end(); ++a){
        bc as_bitset(a->first);
        string bc_str = bc2str(as_bitset, 16);
        string assn = idx2name(a->second, samples);
        double llr = assn_final_llr[a->first];
        char s_d = 'S';
        if (a->second >= samples.size()){
            s_d = 'D';
        }
        fprintf(outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(), assn.c_str(), 
            s_d, llr);
    } 
}

/**
 * Given the set of read counts at each type of SNP per cell
 * and assignments of cells to individuals, sums all ref & alt
 * allele counts at each type of SNP per individual in order to
 * have expectations of those counts conditional on different
 * identities.
 */
void get_expected_allelecounts_from_ids(map<int,
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& expected_allelecounts,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    robin_hood::unordered_map<unsigned long, int>& assn_final,
    unsigned long negdrop_cell){
    
    for (robin_hood::unordered_map<unsigned long, 
        map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >::iterator x = 
        indv_allelecounts.begin(); x != indv_allelecounts.end(); ++x){
        
        if (x->first != negdrop_cell){
            int assn = assn_final[x->first];
            if (expected_allelecounts.count(assn) == 0){
                map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
                expected_allelecounts.insert(make_pair(assn, m));
            }
            for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
                x->second.begin(); y != x->second.end(); ++y){
                if (expected_allelecounts[assn].count(y->first) == 0){
                    map<pair<int, int>, pair<float, float> > m;
                    expected_allelecounts[assn].insert(make_pair(y->first, m));
                }
                for (map<pair<int, int>, pair<float, float> >::iterator z = y->second.begin();
                    z != y->second.end(); ++z){
                    if (expected_allelecounts[assn][y->first].count(z->first) == 0){
                        expected_allelecounts[assn][y->first].insert(make_pair(z->first, make_pair(0,0)));
                    }
                    expected_allelecounts[assn][y->first][z->first].first += z->second.first;
                    expected_allelecounts[assn][y->first][z->first].second += z->second.second;
                }
            }
        }
    }
}

/**
 * Given counts of ref/alt alleles at different types of SNPs in each cell,
 * along with the assignment of each cell to each identity, computes a model
 * of expected counts of ref/alt alleles of each type under each identity.
 * 
 * Then models the sum of ref+alt alleles in all empty droplets as originating
 * from a mixture of each individual in the pool. 
 * 
 * Stores the mapping of individual ID -> proportion in the indvprops data structure
 */
void model_empty_drops(robin_hood::unordered_map<unsigned long, 
    map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    unsigned long negdrop_cell,
    robin_hood::unordered_map<unsigned long, int>& assn_final,
    vector<string>& samples,
    map<int, double>& indvprops){
        
    // Get expected fractions of each type of SNP counts under each identity
    map<int, map<pair<int, int>, map<pair<int, int>, pair<float, float> > > > expected_allelecounts;
    
    get_expected_allelecounts_from_ids(expected_allelecounts,
        indv_allelecounts, assn_final, negdrop_cell);
    
    // Get a vector of all possible assignments 
    vector<int> expected_allelecounts_assn;
    for (map<int, map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >::iterator x = 
        expected_allelecounts.begin(); x != expected_allelecounts.end(); ++x){
        expected_allelecounts_assn.push_back(x->first);
    }
    
    // Use non-negative least squares to model alt allele fractions at 
    // SNPs of each type in the negative drops as originating from a
    // mixture of true identities.
        
    // In this set-up, the fraction of alt alleles at each type of SNP
    // will be a row/equation.
    
    // The right hand side (b) is the fraction of alt alleles at a type
    // of SNP observed in the combination of all negative droplets.
    
    // The left side (A) is a (number of possible identities x number of SNP
    // types) matrix, where each entry is the expected fraction of alt 
    // alleles in SNPs of that category (columns) conditional on being 
    // a specific individual (rows)

    // The result vector (x) gives the inferred fraction of each individual
    // making up the empty droplets.

    vector<vector<double> > A;
    vector<double> b;

    for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
        indv_allelecounts[negdrop_cell].begin(); y != indv_allelecounts[negdrop_cell].end();
        ++y){
        for (map<pair<int, int>, pair<float, float> >::iterator z = y->second.begin(); 
            z != y->second.end(); ++z){
            
            double frac = (double)z->second.second/(double)(z->second.first + z->second.second);
            b.push_back(frac);
            
            vector<double> A_row;
            for (map<int, map<pair<int, int>, map<pair<int, int>, pair<float, float> > > >::iterator e = 
                expected_allelecounts.begin(); e != expected_allelecounts.end(); ++e){
                double frac_this = (double)e->second[y->first][z->first].second / 
                    (double)(e->second[y->first][z->first].first + e->second[y->first][z->first].second);
                A_row.push_back(frac_this);
            } 
            A.push_back(A_row);
        }
    }

    // Make sure the individual fractions sum to 1
    vector<double> A_row_final;
    for (int i = 0; i < expected_allelecounts.size(); ++i){
        A_row_final.push_back(1.0);
    }
    A.push_back(A_row_final);
    b.push_back(1.0);
    
    vector<double> props;
    bool success = nn_lstsq(A, b, props);
    
    if (success){
        // Ensure proportions sum to 1
        // Convert proportions to more human-interpretable format
        double tot = 0.0;
        for (int i = 0; i < props.size(); ++i){
            tot += props[i];
        }
        for (int i = 0; i < props.size(); ++i){
            props[i] /= tot;
        }
        // Also, consider a doublet combination of two individuals
        // to be half of each individual
        for (int i = 0; i < props.size(); ++i){
            int id = expected_allelecounts_assn[i];
            if (id >= samples.size()){
                pair<int, int> comb = idx_to_hap_comb(id, samples.size());
                if (indvprops.count(comb.first) == 0){
                    indvprops.insert(make_pair(comb.first, 0.0));
                }
                if (indvprops.count(comb.second) == 0){
                    indvprops.insert(make_pair(comb.second, 0.0));
                }
                indvprops[comb.first] += 0.5*props[i];
                indvprops[comb.second] += 0.5*props[i];
            }
            else{
                if (indvprops.count(id) == 0){
                    indvprops.insert(make_pair(id, 0.0));
                }
                indvprops[id] += props[i];
            }
        }

    }
}

/**
 * Given a data structure storing the expected proportion of reads 
 * in empty droplets originating from each individual in the pool,
 * along with an output file handle and names of individuals,
 * writes data to the output file on disk.
 */
void dump_empty_drops(FILE* outf, map<int, double>& indvprops, vector<string>& samples){
    // Sort in decreasing order for readability
    vector<pair<double, int> > ipsort;
    for (map<int, double>::iterator ip = indvprops.begin(); ip != indvprops.end();
        ++ip){
        ipsort.push_back(make_pair(-ip->second, ip->first));
    }
    sort(ipsort.begin(), ipsort.end());
    for (int i = 0; i < ipsort.size(); ++i){
        fprintf(outf, "%s\t%f\n", idx2name(ipsort[i].second, samples).c_str(),
            -ipsort[i].first);
    }
}

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_vcf [OPTIONS]\n");
    fprintf(stderr, "Attempts to demultiplex cells in a BAM file based on genotype data \n");
    fprintf(stderr, "for source cell lines provided via VCF. Can infer both single-\n");
    fprintf(stderr, "individual cells and cells better modeled as a combination of two \n");
    fprintf(stderr, "individuals.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --bam -b The BAM file of interest\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants. Only biallelic SNPs \n");
    fprintf(stderr, "       will be considered, and phasing will be ignored.\n");
    fprintf(stderr, "    --output_prefix -o Base name for output files\n");
    fprintf(stderr, "===== RECOMMENDED =====\n");
    fprintf(stderr, "    --cell_barcodes -c To consider only barcodes that have passed \n");
    fprintf(stderr, "       filtering and represent true cells, provide a filtered list \n");
    fprintf(stderr, "       of cell barcodes here. If you mapped the data using \n");
    fprintf(stderr, "       cellranger, this file will be located in \n");
    fprintf(stderr, "       outs/filtered_feature_bc_matrix/barcodes.tsv or \n");
    fprintf(stderr, "       outs/filtered_feature_bc_matrix/barcodes.tsv.gz.\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --doublet_rate -d Prior probability of a cell being a mixture \n");
    fprintf(stderr, "       of two individuals rather than a single individual. \n");
    fprintf(stderr, "       default = 0.5 (set to 0 to disable checking for doublet \n");
    fprintf(stderr, "       combinations and to 1 to disable checking for single indvs\n");
    fprintf(stderr, "    --ids -i If the VCF file contains individuals that you do not \n");
    fprintf(stderr, "       expect to see in your sample, or if you expect to see some \n");
    fprintf(stderr, "       but not other combinations of two individuals that are \n");
    fprintf(stderr, "       present in the VCF, specify all allowed identities in a text \n");
    fprintf(stderr, "       file, one per line. Names of individuals must match those in \n");
    fprintf(stderr, "       the VCF, and combinations of two individuals can be specified \n");
    fprintf(stderr, "       by giving both individual names separated by \"+\" or \"x\", \n");
    fprintf(stderr, "       with names in either order.\n");
    fprintf(stderr, "    --qual -q Minimum variant quality to consider (default 50)\n");
    fprintf(stderr, "    --index_jump -I Instead of reading through the entire BAM file \n");
    fprintf(stderr, "       to count reads at variant positions, use the BAM index to \n");
    fprintf(stderr, "       jump to each variant position. This will be faster if you \n");
    fprintf(stderr, "       have relatively few SNPs, and much slower if you have a lot \n");
    fprintf(stderr, "       of SNPs.\n");
    fprintf(stderr, "    --dump_counts -D After loading variants and counting reads at \n");
    fprintf(stderr, "       variant sites, print write these counts to \n");
    fprintf(stderr, "       [output_prefix].counts and exit. They can then be loaded with \n");
    fprintf(stderr, "        -L on a subsequent run.\n");
    fprintf(stderr, "    --load_counts -L If a previous run was done with -d, load counts \n");
    fprintf(stderr, "       from that file and skip reading through the BAM file (output \n");
    fprintf(stderr, "       prefix must match).\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

int main(int argc, char *argv[]) {    

    static struct option long_options[] = {
       {"bam", required_argument, 0, 'b'},
       {"vcf", required_argument, 0, 'v'},
       {"output_prefix", required_argument, 0, 'o'},
       {"cell_barcode", required_argument, 0, 'c'},
       {"index_jump", no_argument, 0, 'I'},
       {"dump_counts", no_argument, 0, 'D'},
       {"load_counts", no_argument, 0, 'L'},
       {"ids", required_argument, 0, 'i'},
       {"qual", required_argument, 0, 'q'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string bamfile;
    string vcf_file;
    bool cell_barcode = false;
    string cell_barcode_file = "";
    bool stream = true; // Opposite of index-jumping
    bool dump_counts = false;
    bool load_counts = false;
    string output_prefix = "";
    // If a cell barcode list was provided, all counts for
    // "negative droplets" -- barcodes missing from the list --
    // will be summed and associated with a special barcode
    bool neg_drops = true;
    int vq = 50;
    string idfile;
    bool idfile_given = false;
    double doublet_rate = 0.5;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "b:v:o:c:i:q:d:IDLh", long_options, &option_index )) != -1){
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
            case 'c':
                cell_barcode = true;
                cell_barcode_file = optarg;
                break;
            case 'd':
                doublet_rate = atof(optarg);
                break;
            case 'i':
                idfile_given = true;
                idfile = optarg;
                break;
            case 'q':
                vq = atoi(optarg);
                break;
            case 'I':
                stream = false;
                break;
            case 'D':
                dump_counts = true;
                break;
            case 'L':
                load_counts = true;
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
    if (bamfile.length() == 0 && !load_counts){
        fprintf(stderr, "ERROR: bam file (--bam) required\n");
        exit(1);
    }
    if (vcf_file.length() == 0){
        fprintf(stderr, "ERROR: vcf/BCF file (-v) required\n");
        exit(1);
    }
    if (output_prefix.length() == 0){
        fprintf(stderr, "ERROR: output_prefix/-o required\n");
        exit(1);
    }
    if (dump_counts && load_counts){
        fprintf(stderr, "ERROR: only one of --dump_counts and --load_counts is allowed\n");
        exit(1);
    }
    if (doublet_rate < 0 || doublet_rate > 1){
        fprintf(stderr, "ERROR: doublet rate must be between 0 and 1, inclusive.\n");
        exit(1);
    }
    
    // Initialize barcode for negative drops - all A should be safe
    string negdrop_bc_str = "";
    for (int i = 0; i < 16; ++i){
        negdrop_bc_str += "A";
    } 
    bc negdrop_bc_bc;
    str2bc(negdrop_bc_str.c_str(), negdrop_bc_bc, 16);
    unsigned long negdrop_cell = negdrop_bc_bc.to_ulong();

    // Init BAM reader
    bam_reader reader = bam_reader();
    if (!load_counts){
        // We will be actually using the BAM reader, so we need to initialize properly.
        reader.set_file(bamfile);
    }
    
    // Store the names of all individuals in the VCF
    vector<string> samples;
    
    // Store data about SNPs
    map<int, map<int, var> > snpdat;
    
    // Load variant data from VCF file
    if (load_counts){
        fprintf(stderr, "Reading VCF/BCF header...\n");
    }
    else{
        fprintf(stderr, "Loading variants from VCF/BCF...\n");
    }
    int nsnps = read_vcf(vcf_file, reader, samples, snpdat, vq, load_counts);
    
    set<int> allowed_ids;
    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile.c_str());
        }
    }
    
    set<unsigned long> cell_barcodes;
    if (cell_barcode){
        parse_barcode_file(cell_barcode_file, cell_barcodes); 
    }

    // At this point, give up on neg drops if we do not have a cell barcode list.
    if (cell_barcodes.size() == 0){
        neg_drops = false;
    }
    
    // Data structure to store allele counts at SNPs of each possible type.
    // SNPs are defined by their allelic state in each pair of 2 individuals.
    
    // In other words, one category is all SNPs at which individual A is 
    // homozygous for the reference allele and individual B is heterozygous.
    // Another category is all SNPs at which individual C is homozygous for
    // the alt allele and individual D is homozygous for the reference allele.
    
    // For each cell, we sum reference and alt allele counts across all SNPs
    // falling into each category.
    
    // Here, cell barcodes are represented as unsigned long (see bc_hash.cpp in 
    // htswrapper library)
    
    // Both pair<int, int> are keys of (individual ID, nalt) where nalt is 
    // number of alt alleles for that individual in SNPs of that type
    
    // Each combination is stored only once - lower index individual always
    // comes first

    // The final two numbers are ref and alt allele counts -- where each count
    // is actually Probability(mapping of read is correct), determined by
    // map quality

    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;

    // Store counts for currently-tracked SNPs
    map<int, map<int, robin_hood::unordered_map<unsigned long, 
        pair<float, float> > > > varcounts_site;
    
    robin_hood::unordered_map<unsigned long, map<int, float> > fracs;
    
    // Print progress message every n sites
    int progress = 1000;
    // What was the last number of sites for which a message was printed?
    int last_print = 0;

    if (load_counts){

        // Figure out the appropriate file name from the previous run and
        // load the counts, instead of reading through the BAM file.
        string fn = output_prefix + ".counts";
        fprintf(stderr, "Loading counts...\n");
        load_counts_from_file(indv_allelecounts, varcounts_site, samples, fn);
        
        // If the negative droplet-representative cell index was 
        // present in the loaded counts, we can do stuff with it.
        if (indv_allelecounts.count(negdrop_cell) > 0){
            neg_drops = true;
        } 

    }
    else{
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
                            if (snpdat[curtid].count(cursnp->first) > 0){        
                                dump_vcs_counts(varcounts_site[curtid][cursnp->first], 
                                    indv_allelecounts,
                                    cursnp->second, samples.size());
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
                    
                    if (varcounts_site[reader.tid()].count(cursnp->first) > 0){
                        dump_vcs_counts(varcounts_site[reader.tid()][cursnp->first], 
                            indv_allelecounts, 
                            cursnp->second, samples.size());
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
                    process_bam_record(reader, cursnp2->first, cursnp2->second,
                        varcounts_site[reader.tid()], cell_barcode, cell_barcodes, 
                        neg_drops, negdrop_cell);
                    ++cursnp2;
                }
                if ((nsnp_processed - last_print) % progress == 0){
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
                        process_bam_record(reader, cursnp->first, cursnp->second, 
                            varcounts_site[tid], cell_barcode, cell_barcodes, 
                            neg_drops, negdrop_cell);
                    }
                    
                    dump_vcs_counts(varcounts_site[tid][cursnp->first], 
                        indv_allelecounts, cursnp->second, samples.size());
                    
                    ++nsnp_processed;        
                
                    if ((nsnp_processed - last_print) % progress == 0){
                        fprintf(stderr, "Processed %d of %d SNPs\r", nsnp_processed, 
                            nsnps); 
                        last_print = nsnp_processed;
                    }
                }
            }
        }
        fprintf(stderr, "\n");
        
        if (dump_counts){
            // Write the data just compiled to disk.

            string fname = output_prefix + ".counts";
            FILE* outf = fopen(fname.c_str(), "w");   
            fprintf(stderr, "Writing allele counts to disk...\n");
            dump_cellcounts(outf, indv_allelecounts, samples);
            fprintf(stderr, "Done\n");
            fclose(outf);
        }
    }
        
    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn_final;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_final_llr;
    
    fprintf(stderr, "Finding likeliest identities of cells...\n");
    // Get assignments of cell barcodes
    assign_ids(indv_allelecounts, samples, assn_final, assn_final_llr, 
        allowed_ids, neg_drops, negdrop_cell, doublet_rate);
    
    // Write these best assignments to disk
    {
        string fname = output_prefix + ".assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing cell-individual assignments to disk...\n");
        dump_assignments(outf, assn_final, assn_final_llr, samples);
        fclose(outf);
    }
    
    // Finally, if a cell barcode list was provided, and given the allele counts 
    // at cells assigned to each individual, we can attempt to model the 
    // sum of all empty droplets as a mixture of genotypes of different individuals
    // and write that to a file.
    if (neg_drops){
        fprintf(stderr, "Modeling empty droplets as a mixture of individuals \
present in the pool...\n");
        map<int, double> indvprops;
        model_empty_drops(indv_allelecounts, negdrop_cell, assn_final, samples, indvprops);
        string fname = output_prefix + ".emptydrops";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_empty_drops(outf, indvprops, samples);
        fclose(outf);
    }
    
    return 0;
}
