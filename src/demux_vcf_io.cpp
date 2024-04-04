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
#include <zlib.h>
#include <htswrapper/bc.h>
#include <htswrapper/gzreader.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * ===== Contains functions relating to reading/writing text files ===== 
 * ===== needed by demux_vcf.                                      =====
 */

/**
 * Read file of allowed ID assignments
 */
void parse_idfile(string& idfile, 
    vector<string>& samples,
    set<int>& ids_allowed,
    set<int>& ids_allowed2,
    bool add_all_doublets){

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

    // We will also store a copy (ids_allowed2) without these singlets added in
    ids_allowed2 = ids_allowed;
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
    if (add_all_doublets){
        // We also need to add in every possible combination of two individuals provided.
        
        set<int> adds;
        vector<int> singles;
        for (set<int>::iterator id = ids_allowed.begin(); id != ids_allowed.end(); ++id){
            if (*id < samples.size()){
                singles.push_back(*id);
            }
        }
        sort(singles.begin(), singles.end());
        for (int i = 0; i < singles.size()-1; ++i){
            int idx1 = singles[i];
            for (int j = i + 1; j < singles.size(); ++j){
                int idx2 = singles[j];
                int k = hap_comb_to_idx(idx1, idx2, samples.size());
                if (ids_allowed.find(k) == ids_allowed.end()){
                    adds.insert(k);
                }
            }
        }
        for (set<int>::iterator k = adds.begin(); k != adds.end(); ++k){
            ids_allowed.insert(*k);
        }

        // If we just added all possible doublet combinations, then we couldn't have
        // started with a filtered list that was meant to include specific doublet types.
        // Therefore, ids_allowed2 should match ids_allowed.
        ids_allowed2 = ids_allowed;
    }
}

void write_samples(string& filename, vector<string>& samples){
    FILE* outf = fopen(filename.c_str(), "w");
    for (int i = 0; i < samples.size(); ++i){
        fprintf(outf, "%s\n", samples[i].c_str());
    }
    fclose(outf);
}

void load_samples(string& filename, vector<string>& samples){
    ifstream inf(filename.c_str());
    string samp;
    while (inf >> samp){
        samples.push_back(samp);
    }
}

void write_allowed(string& filename, set<int>& allowed, vector<string>& samples){
    FILE* outf = fopen(filename.c_str(), "w");
    for (set<int>::iterator a = allowed.begin(); a != allowed.end(); ++a){
        fprintf(outf, "%s\n", samples[*a].c_str());
    }
    fclose(outf);
}

void load_allowed(string& filename, set<int>& allowed, vector<string>& samples){
    map<string, int> samp2idx;
    for (int i = 0; i < samples.size(); ++i){
        samp2idx.insert(make_pair(samples[i], i));
    }
    ifstream inf(filename.c_str());
    string samp;
    while (inf >> samp){
        size_t sep_pos = samp.find("+");
        if (sep_pos == string::npos){
            sep_pos = samp.find("x");
        }
        if (sep_pos != string::npos){
            string samp1 = samp.substr(0, sep_pos);
            string samp2 = samp.substr(sep_pos+1, samp.length()-sep_pos-1);
            if (samp2idx.count(samp1) > 0){
                if (samp2idx.count(samp2) > 0){
                    int idx1 = samp2idx[samp1];
                    int idx2 = samp2idx[samp2];
                    if (idx1 == idx2){
                        allowed.insert(idx1);
                    }
                    else if (idx2 < idx1){
                        int tmp = idx1;
                        idx1 = idx2;
                        idx2 = tmp;
                        int idx3 = hap_comb_to_idx(idx1, idx2, samples.size());
                        allowed.insert(idx3);
                        // Also allow either component individual
                        allowed.insert(idx1);
                        allowed.insert(idx2);
                    }        
                }
                else{
                    fprintf(stderr, "WARNING: sample %s not found in data\n", samp2.c_str());
                }
            }
            else{
                fprintf(stderr, "WARNING: sample %s not found in data\n", samp1.c_str());
            }
        }
        else{
            if (samp2idx.count(samp) > 0){
                allowed.insert(samp2idx[samp]);
            }
            else{
                fprintf(stderr, "WARNING: sample %s not found in data\n", samp.c_str());
            }
        }
    }
}

/**
 * If a previous run was dumped to count files, load those counts instead of 
 * re-processing the BAM file.
 */
void load_counts_from_file(
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    vector<string>& indvs,   
    string& filename){
   
    gzreader reader(filename);
    while(reader.next()){
        istringstream splitter(reader.line);
        string field;
        int idx = 0;

        unsigned long cell;
        int indv1;
        int indv2;
        int type1;
        int type2;
        float ref;
        float alt;
    
        pair<int, int> indv1key;
        pair<int, int> indv2key;

        while(getline(splitter, field, '\t')){
            if (idx == 0){
                // cell
                cell = atol(field.c_str());
                if (indv_allelecounts.count(cell) == 0){
                    map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
                    indv_allelecounts.emplace(cell, m);
                }
            }
            else if (idx == 1){
                // indv 1
                indv1 = atoi(field.c_str());
            }
            else if (idx == 2){
                // type 1
                type1 = atoi(field.c_str());
                indv1key = make_pair(indv1, type1);
                    if (indv_allelecounts[cell].count(indv1key) == 0){
                        map<pair<int, int>, pair<float, float> > m;
                        indv_allelecounts[cell].insert(make_pair(indv1key, m));
                    }
            }
            else if (idx == 3){
                // indv 2
                indv2 = atoi(field.c_str());
            }
            else if (idx == 4){
                // type 2
                type2 = atoi(field.c_str());
                indv2key = make_pair(indv2, type2);
            }
            else if (idx == 5){
                // ref count
                ref = atof(field.c_str());
            }
            else if (idx == 6){
                // alt count
                alt = atof(field.c_str());
                indv_allelecounts[cell][indv1key].insert(make_pair(indv2key, make_pair(ref, alt)));
            }
            ++idx;
        }
    }
    
    /*
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
    */
}

/**
 * Print counts to text files.
 */
void dump_cellcounts(gzFile& out_cell,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts, 
    vector<string>& samples){
    
    char linebuf[1024];

    for (robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >::iterator x = indv_allelecounts.begin();
        x != indv_allelecounts.end(); ++x){
        
        for (map<pair<int, int>, map<pair<int, int>, pair<float, float> > >::iterator y = 
                x->second.begin(); y != x->second.end(); ++y){
            //string indv = idx2name(y->first.first, samples);
            int nalt = y->first.second;
            for (map<pair<int, int>, pair<float, float> >::iterator z = 
                y->second.begin(); z != y->second.end(); ++z){
                //string indv2 = "NA";
                //if (z->first.first != -1){
                //    indv2 = idx2name(z->first.first, samples);
                //}
                int nalt2 = z->first.second;
                
                sprintf(&linebuf[0], "%ld\t%d\t%d\t%d\t%d\t%f\t%f\n", x->first, y->first.first,
                    nalt, z->first.first, nalt2, z->second.first, z->second.second);
                
                gzwrite(out_cell, &linebuf[0], strlen(linebuf));
                   
                //fprintf(out_cell, "%ld\t%d\t%d\t%d\t%d\t%f\t%f\n", x->first, y->first.first,
                //   nalt, z->first.first, nalt2, z->second.first, z->second.second);

                //fprintf(out_cell, "%ld\t%s\t%d\t%s\t%d\t%f\t%f\n", x->first, indv.c_str(),
                //    nalt, indv2.c_str(), nalt2, z->second.first, z->second.second);
            }
        }
    } 
}

/**
 * Read in previously computed expected matching fractions conditional
 * on cell identities.
 */
void load_exp_fracs(string& filename,   
    map<pair<int, int>, map<int, float> >& conditional_match_frac,
    vector<string>& samples){
    
    map<string, int> sample2idx;
    for (int i = 0; i < samples.size(); ++i){
        sample2idx.insert(make_pair(samples[i], i));
    }

    ifstream infile(filename.c_str());
    string sample1;
    int type;
    string sample2;
    float frac; 
    while (infile >> sample1 >> type >> sample2 >> frac){
        if (sample2idx.count(sample1) == 0){
            fprintf(stderr, "ERROR: indv %s from expected fraction file not found in VCF\n", sample1.c_str());
            exit(1);
        }
        if (sample2idx.count(sample2) == 0){
            fprintf(stderr, "ERROR: indv %s from expected fraction file not found in VCF\n", sample2.c_str());
            exit(1);
        }       
        else{
            int idx1 = sample2idx[sample1];
            int idx2 = sample2idx[sample2];

            pair<int, int> key1 = make_pair(idx1, type);
            if (conditional_match_frac.count(key1) == 0){
                map<int, float> m;
                conditional_match_frac.insert(make_pair(key1, m));
            }
            conditional_match_frac[key1].insert(make_pair(idx2, frac));
        }
    }
}

/**
 * Write conditional match fracs to output file.
 */
void dump_exp_fracs(FILE* exp_frac_out,
    map<pair<int, int>, map<int, float> >& conditional_match_frac,
    vector<string>& samples){ 
    for (map<pair<int, int>, map<int, float> >::iterator cmf = conditional_match_frac.begin();
        cmf != conditional_match_frac.end(); ++cmf){
        for (map<int, float>::iterator cmf2 = cmf->second.begin(); cmf2 != cmf->second.end(); ++cmf2){
            fprintf(exp_frac_out, "%s\t%d\t%s\t%f\n", idx2name(cmf->first.first, samples).c_str(),
                cmf->first.second, idx2name(cmf2->first, samples).c_str(), cmf2->second);
        }
    }
}

/**
 * Spill cell -> individual assignments to disk
 */
void dump_assignments(FILE* outf,
    robin_hood::unordered_map<unsigned long, int>& assn_final,
    robin_hood::unordered_map<unsigned long, double>& assn_final_llr,
    vector<string>& samples,
    string& barcode_group){

    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn_final.begin();
        a != assn_final.end(); ++a){
        bc as_bitset(a->first);
        string bc_str = bc2str(as_bitset);
        if (barcode_group != ""){
            bc_str += "-" + barcode_group;
        }
        string assn = idx2name(a->second, samples);
        double llr = assn_final_llr[a->first];
        // Only write assignments that are better than an even guess
        if (llr > 0){
            char s_d = 'S';
            if (a->second >= samples.size()){
                s_d = 'D';
            }
            fprintf(outf, "%s\t%s\t%c\t%f\n", bc_str.c_str(), assn.c_str(), 
                s_d, llr);
        }
    } 
}

/**
 * Load cell -> individual assignments from text file
 */
void load_assignments_from_file(string& filename,
    robin_hood::unordered_map<unsigned long, int>& assn,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    vector<string>& samples,
    string& barcode_group){
    
    map<string, int> sample2idx;
    for (int i = 0; i < samples.size(); ++i){
        sample2idx.insert(make_pair(samples[i], i));
    }

    ifstream inf(filename.c_str());
    string bc_str;
    string idstr;
    char s_d;
    double llr;
    
    barcode_group = "";

    while (inf >> bc_str >> idstr >> s_d >> llr){
        size_t dashpos = bc_str.rfind("-");
        if (dashpos != string::npos){
            // Assumes all "barcode_group" signifiers are the same. We shouldn't be 
            // trying to quanitfy contamination across more than one lane at a time.
            barcode_group = bc_str.substr(dashpos + 1, bc_str.length()-dashpos-1);
            bc_str = bc_str.substr(0, dashpos);
        }
        bc as_bitset;
        str2bc(bc_str.c_str(), as_bitset);
        unsigned long as_ul = as_bitset.to_ulong();

        int indv_idx = -1;
        size_t splitpos = idstr.find("+");
        if (splitpos == string::npos){
            splitpos = idstr.find("x");
        }
        if (splitpos != string::npos){
            string id1 = idstr.substr(0, splitpos);
            string id2 = idstr.substr(splitpos+1, idstr.length()-splitpos-1);
            if (sample2idx.count(id1) == 0){
                fprintf(stderr, "ERROR: individual %s not found in VCF\n", id1.c_str());
                exit(1);
            }
            if (sample2idx.count(id2) == 0){
                fprintf(stderr, "ERROR: individual %s not found in VCF\n", id2.c_str());
            }
            int idx1 = sample2idx[id1];
            int idx2 = sample2idx[id2];
            if (idx1 > idx2){
                indv_idx = hap_comb_to_idx(idx2, idx1, samples.size());
            }
            else{
                indv_idx = hap_comb_to_idx(idx1, idx2, samples.size());
            }
        }
        else{
            if (s_d == 'D'){
                fprintf(stderr, "ERROR: ID %s has no delimiter but labeled as double\n", idstr.c_str());
                exit(1);
            }
            if (sample2idx.count(idstr) == 0){
                fprintf(stderr, "ERROR: individual %s not found in VCF\n", idstr.c_str());
                exit(1);
            }
            indv_idx = sample2idx[idstr];
        }
        assn.emplace(as_ul, indv_idx);
        assn_llr.emplace(as_ul, llr);
    }
}

/**
 * Write summary information about a completed run to output file.
 */
void write_summary(FILE* outf, 
    string& outpre,
    robin_hood::unordered_map<unsigned long, int>& assn,
    vector<string>& samples,
    double error_ref,
    double error_alt,
    double error_sigma,
    double ref_mm_rate_posterior,
    double alt_mm_rate_posterior,
    string& vcf_file,
    int vq_filter,
    double doublet_rate,
    map<int, double>& p_ncell,
    map<int, double>& p_llr){
    
    if (vcf_file != ""){
        fprintf(outf, "%s\tparam\tvcf_file\t%s\n", outpre.c_str(),
            vcf_file.c_str());
        fprintf(outf, "%s\tparam\tvariant_qual\t%d\n", outpre.c_str(),
            vq_filter);
    }
    fprintf(outf, "%s\tparam\tdoublet_rate\t%f\n", outpre.c_str(),
        doublet_rate);
    fprintf(outf, "%s\tparam\terror_ref\t%f\n", outpre.c_str(),
        error_ref);
    fprintf(outf, "%s\tparam\terror_alt\t%f\n", outpre.c_str(),
        error_alt);
    fprintf(outf, "%s\tparam\terror_sigma\t%f\n", outpre.c_str(),
        error_sigma);
    
    fprintf(outf, "%s\tresult\terror_ref_posterior\t%f\n", outpre.c_str(),
        ref_mm_rate_posterior);
    fprintf(outf, "%s\tresult\terror_alt_posterior\t%f\n", outpre.c_str(),
        alt_mm_rate_posterior);
    
    int count_doublets = 0;
    int tot_singlets = 0;
    map<int, int> idcounts;
    map<int, int> idcounts_in_doublet;
    for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
        a != assn.end(); ++a){
        if (a->second >= samples.size()){
            count_doublets++;
            pair<int, int> combo = idx_to_hap_comb(a->second, samples.size());
            if (idcounts_in_doublet.count(combo.first) == 0){
                idcounts_in_doublet.insert(make_pair(combo.first, 0));
            }
            if (idcounts_in_doublet.count(combo.second) == 0){
                idcounts_in_doublet.insert(make_pair(combo.second, 0));
            }
            idcounts_in_doublet[combo.first]++;
            idcounts_in_doublet[combo.second]++;
        }
        else{
            tot_singlets++;
        }
        if (idcounts.count(a->second) == 0){
            idcounts.insert(make_pair(a->second, 0));
        }
        idcounts[a->second]++;
    }
    
    fprintf(outf, "%s\tresult\ttot_cells\t%d\n", outpre.c_str(), tot_singlets+count_doublets);
    if (count_doublets > 0){
        fprintf(outf, "%s\tresult\tfrac_doublets\t%f\n", outpre.c_str(), 
            (double)count_doublets / (double)assn.size());
        
        double doub_chisq = doublet_chisq(idcounts, samples.size());
        if (doub_chisq >= 0.0){
            fprintf(outf, "%s\tresult\tdoublet_chisq.p\t%f\n", outpre.c_str(),
                doub_chisq);
        }

        vector<pair<double, int> > fdisort;
        for (map<int, int>::iterator icd = idcounts_in_doublet.begin();
            icd != idcounts_in_doublet.end(); ++icd){
            fdisort.push_back(make_pair(-(double)icd->second/(double)count_doublets, icd->first));
        } 
        sort(fdisort.begin(), fdisort.end());
        for (int i = 0; i < fdisort.size(); ++i){
            fprintf(outf, "%s\tperc_doublets_include\t%s\t%f\n", outpre.c_str(),
                idx2name(fdisort[i].second, samples).c_str(),
                -fdisort[i].first);
        }
    }
    
    for (map<int, double>::iterator pn = p_ncell.begin(); pn != p_ncell.end();
        ++pn){
        fprintf(outf, "%s\tp_fewer_cells\t%s\t%f\n", outpre.c_str(),
            idx2name(pn->first, samples).c_str(), pn->second);
        fprintf(outf, "%s\tp_lower_LLRs\t%s\t%f\n", outpre.c_str(),
            idx2name(pn->first, samples).c_str(), p_llr[pn->first]);
    }

    vector<pair<int, int> > idcsort;
    for (map<int, int>::iterator ic = idcounts.begin(); ic != idcounts.end(); ++ic){
        idcsort.push_back(make_pair(-ic->second, ic->first));
    }
    sort(idcsort.begin(), idcsort.end());
    for (int i = 0; i < idcsort.size(); ++i){
        fprintf(outf, "%s\tnum_cells\t%s\t%d\n", outpre.c_str(), 
            idx2name(idcsort[i].second, samples).c_str(), -idcsort[i].first);
    }
}

/**
 * Write contamination profile (likeliest makeup of samples from VCF) to disk
 */
void dump_contam_prof(FILE* outf,
    map<int, double>& contam_prof,
    vector<string>& samples){
    
    vector<pair<double, string> > cpsort;
    for (map<int, double>::iterator cp = contam_prof.begin(); cp != contam_prof.end(); ++cp){
        string name;
        if (cp->first < 0){
            name = "other_species";
        }
        else{
            name = idx2name(cp->first, samples);
        }
        cpsort.push_back(make_pair(-cp->second, name));
    }
    sort(cpsort.begin(), cpsort.end());
    for (int i = 0; i < cpsort.size(); ++i){
        fprintf(outf, "%s\t%f\n", cpsort[i].second.c_str(), -cpsort[i].first);
    }
}

/**
 * Write contamination rate & SE per cell to disk
 */
void dump_contam_rates(FILE* outf,
    robin_hood::unordered_map<unsigned long, double>& contam_rate,
    robin_hood::unordered_map<unsigned long, double>& contam_rate_se,
    vector<string>& samples,
    string& barcode_group){
    
    for (robin_hood::unordered_map<unsigned long, double>::iterator cr = contam_rate.begin();
       cr != contam_rate.end(); ++cr){
        bc as_bitset(cr->first);
        string bc_str = bc2str(as_bitset);
        if (barcode_group != ""){
            bc_str += "-" + barcode_group;
        }
        fprintf(outf, "%s\t%f\t%f\n", bc_str.c_str(), cr->second, contam_rate_se[cr->first]);
    }
}

/**
 * Write to disk the set of alt allele frequencies inferred to exist in ambient RNA
 * in a given data set.
 */
void dump_amb_fracs(FILE* outf, 
    map<pair<int, int>, map<pair<int, int>, double> >& amb_mu){

    for (map<pair<int, int>, map<pair<int, int>, double> >::iterator x = amb_mu.begin();
        x != amb_mu.end(); ++x){
        for (map<pair<int, int>, double>::iterator y = x->second.begin(); y != x->second.end();
            ++y){
            fprintf(outf, "%d\t%d\t%d\t%d\t%f\n", x->first.first, x->first.second,
                y->first.first, y->first.second, y->second);
        }
    } 
}

