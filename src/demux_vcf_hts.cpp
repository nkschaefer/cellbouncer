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
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_hts.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * ===== Contains functions relating to processing HTSlib-format files =====
 */

// ===== VCF-related functions =====

/**
 * Read variant data from VCF.
 */
int read_vcf(string& filename, 
    bam_reader& reader,
    vector<string>& samples,
    map<int, map<int, var> >& snps,
    int min_vq,
    bool hdr_only,
    bool skip_seq2tid,
    bool allow_missing){
    
    // Map sequence names to TIDs for storage
    map<string, int> seq2tid;
    if (!hdr_only && !skip_seq2tid){
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
                    fprintf(stderr, "ERROR loading genotypes at %s %ld\n", 
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
                
                if (allow_missing || nmiss == 0){
                    snps[tid].insert(make_pair(pos, v));
                }
                else{
                    --nvar;
                }
            }
        }
    }
    fprintf(stderr, "Loaded %ld variants\n", nvar);
    
    bcf_hdr_destroy(bcf_header);
    bcf_destroy(bcf_record);
    bcf_close(bcf_reader);
    
    return nvar;
}

/**
 * Given variants from a VCF file, gets the expected number of ref, het, and alt
 * alleles per each other individual given each possible true identity.
 */
void get_conditional_match_fracs(map<int, map<int, var> >& snpdat,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs, 
    int n_samples){
    
    // (sample, sitetype) -> other_sample -> number
    map<pair<int, int>, map<int, float> > conditional_match_tots;

    float err_tol = 0.0;

    for (int i = 0; i < n_samples; ++i){
        for (int n = 0; n <= 2; ++n){
            map<int, float> m;
            pair<int, int> key = make_pair(i, n);
            conditional_match_fracs.insert(make_pair(key, m));
            conditional_match_tots.insert(make_pair(key, m));
        }
    }
    for (map<pair<int, int>, map<int, float> >::iterator cmf = conditional_match_fracs.begin();
        cmf != conditional_match_fracs.end(); ++cmf){
        for (int j = 0; j < n_samples; ++j){
            if (j != cmf->first.first){
                cmf->second.insert(make_pair(j, 0.0));
            }
        }
    }
    
    for (map<int, map<int, var> >::iterator s1 = snpdat.begin(); s1 != snpdat.end(); ++s1){
        for (map<int, var>::iterator s2 = s1->second.begin(); s2 != s1->second.end(); ++s2){
            for (int i = 0; i < n_samples-1; ++i){
                int n_i = 0;
                if (s2->second.haps_covered[i]){
                    if (s2->second.haps1[i]){
                        n_i++;
                    }
                    if (s2->second.haps2[i]){
                        n_i++;
                    }
                    for (int j = i + 1; j < n_samples; ++j){
                        int n_j = 0;
                        if (s2->second.haps_covered[j]){
                            if (s2->second.haps1[j]){
                                n_j++;
                            }
                            if (s2->second.haps2[j]){
                                n_j++;
                            }
                            pair<int, int> k1 = make_pair(i, n_i);
                            pair<int, int> k2 = make_pair(j, n_j);
                            conditional_match_tots[k1][j]++;
                            conditional_match_tots[k2][i]++;

                            conditional_match_fracs[k1][j] += (float)n_j/2.0;
                            conditional_match_fracs[k2][i] += (float)n_i/2.0;                                                 
                        }                            
                    }
                }
            }
        }
    }
    
    // normalize
    for (map<pair<int, int>, map<int, float> >::iterator cmf = conditional_match_fracs.begin();
        cmf != conditional_match_fracs.end(); ++cmf){
        
        for (map<int, float>::iterator cmf2 = cmf->second.begin(); cmf2 != cmf->second.end(); ++cmf2){
            if (conditional_match_tots[cmf->first][cmf2->first] == 0.0){
                cmf2->second = err_tol;
            }
            else{
                cmf2->second /= conditional_match_tots[cmf->first][cmf2->first];
                if (cmf2->second == 0){
                    cmf2->second = err_tol;
                }
                else if (cmf2->second == 1.0){
                    cmf2->second = 1.0-err_tol;
                }
            }
        }
    }
    
    // add in self-self
    for (int i = 0; i < n_samples; ++i){
        conditional_match_fracs[make_pair(i, 0)].insert(make_pair(i, err_tol));
        conditional_match_fracs[make_pair(i, 1)].insert(make_pair(i, 0.5));
        conditional_match_fracs[make_pair(i, 2)].insert(make_pair(i, 1.0-err_tol));
    }
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
    set<unsigned long>& bcs_valid){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
                        
        // Get BC key
        bc bc_bits;
        str2bc(reader.cb_z, bc_bits);
        unsigned long bc_key = bc_bits.to_ulong();
        
        if (!has_bc_list || bcs_valid.find(bc_key) != bcs_valid.end()){
            
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
 * Same as above, but for counting alt frequencies at individual SNPs in bulk data.
 *
 * Increments err_count if the read has a non-ref/alt allele.
 */
void process_bam_record_bulk(bam_reader& reader,
    int snppos,
    var& vardat,
    map<int, map<int, pair<float, float> > >& snp_ref_alt,
    map<int, map<int, float> >& snp_err){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
         
        // Instead of storing actual read counts, store the probability
        // that the mapping was correct.
        float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);
        
        int tid = reader.tid();
        if (snp_ref_alt.count(tid) == 0){
            map<int, pair<float, float> > m;
            snp_ref_alt.insert(make_pair(tid, m));
            map<int, float> m2;
            snp_err.insert(make_pair(tid, m2));
        }
        if (snp_ref_alt[tid].count(snppos) == 0){
            snp_ref_alt[tid].insert(make_pair(snppos, make_pair(0.0, 0.0)));
            snp_err[tid].insert(make_pair(snppos, 0.0));
        }

        // Note: this function expects positions to be 1-based, but 
        // BCF/BAM functions store as 0-based
        char allele = reader.get_base_at(snppos + 1);
        
        if (allele != 'N' && allele != '-'){
            if (allele == vardat.ref){
                snp_ref_alt[tid][snppos].first += prob_corr;
            }
            else if (allele == vardat.alt){
                snp_ref_alt[tid][snppos].second += prob_corr;
            }
            else{
                snp_err[tid][snppos] += prob_corr;
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
