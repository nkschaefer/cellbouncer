#define STRINGIZE_(x) #x
#define STRINGIZE(x) STRINGIZE_(x)
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
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/seq_fuzzy_match.h>
#include <htswrapper/bc_scanner.h>
#include <mixtureDist/functions.h>
#include <mixtureDist/mixtureDist.h>
#include <mixtureDist/mixtureModel.h>
#include "common.h"


using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "demux_multiseq [OPTIONS]\n");
    fprintf(stderr, "Given reads containing MULTIseq data, counts the occurrences of each\n");
    fprintf(stderr, "  MULTIseq barcode per cell barcode and assigns each cell barcode\n");
    fprintf(stderr, "  an identity based on MULTIseq barcode meaning.\n");
    fprintf(stderr, "A single run will produce a .counts file and an .assignments file. If\n");
    fprintf(stderr, "  you run again with the same output_prefix, the counts from the prior\n");
    fprintf(stderr, "  run will be loaded.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "\n   ===== REQUIRED =====\n");
    fprintf(stderr, "   --output_prefix -o The prefix for output file names. Will create\n");
    fprintf(stderr, "       a file ending in .counts and another ending in .assignments.\n");
    fprintf(stderr, "\n   ===== REQUIRED UNLESS LOADING COUNTS FROM A PREVIOUS RUN =====\n");
    fprintf(stderr, "   --read1 -1 Forward read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --read2 -2 Reverse read file for MULTIseq data. Can specify multiple times.\n");
    fprintf(stderr, "   --whitelist -w Cell barcode whitelist file (i.e. included with 10X Genomics\n");
    fprintf(stderr, "       Cellranger software in cellranger-x.y.z/lib/python/cellranger/barcodes/\n");
    fprintf(stderr, "\n   ===== STRONGLY RECOMMENDED =====\n");
    fprintf(stderr, "   --mapping -m 2-column tab separated file mapping string interpretation of\n");
    fprintf(stderr, "       MULTIseq label (i.e. unique identifier, or treatment name) to MULTIseq\n");
    fprintf(stderr, "       well identifier.\n");
    fprintf(stderr, "   --cell_barcodes -B The filtered list of valid cell barcodes, i.e. from cellranger or\n");
    fprintf(stderr, "       STARsolo. Only cells in this list will be processed.\n");
    fprintf(stderr, "\n   ===== OPTIONAL =====\n");
    fprintf(stderr, "   --help -h Display this message and exit.\n");
    
    fprintf(stderr, "   --barcodes -b A path to a file listing MULTIseq well IDs and corresponding\n");
    fprintf(stderr, "       barcodes, tab separated. If you do not provide one, the default file (in\n");
    fprintf(stderr, "       the data directory) will be used.\n");
    fprintf(stderr, "   --mismatches -M The number of allowable mismatches to a MULTIseq barcode to accept\n");
    fprintf(stderr, "       it. Default = -1 (take best match overall, if there are no ties)\n");
    fprintf(stderr, "   --doublet_rate -D What is the prior expected doublet rate?\n");
    fprintf(stderr, "       (OPTIONAL; default = 0.1). Must be a decimal between 0 and 1,\n");
    fprintf(stderr, "       exclusive.\n");
    fprintf(stderr, "   --output_unsassigned -u Print entries for unassigned barcodes in assignments file\n");
    exit(code);
}

/**
 * Returns barcode length.
 */
int load_well2bc(string& filename, 
    map<unsigned long, string>& bc2well, 
    vector<string>& bclist){
    
    ifstream inf(filename);
    string bc_str;
    string well;
    bc bc_bitset;
    int bc_len = -1;
    while (inf >> well >> bc_str){
        if (bc_len == -1){
            bc_len = bc_str.length();
        }
        else if (bc_len != bc_str.length()){
            fprintf(stderr, "ERROR: mismatching MULTIseq barcode lengths in file %s:\n",
                filename.c_str());
            fprintf(stderr, "%d vs %ld\n", bc_len, bc_str.length());
            exit(1);
        }
        str2bc(bc_str.c_str(), bc_bitset, bc_len);
        bc2well.insert(make_pair(bc_bitset.to_ulong(), well));
        bclist.push_back(bc_str);
    }
    return bc_len;
}

void load_well_mapping(string& filename, map<string, string>& well2id){
    ifstream inf(filename);
    string well;
    string id;
    while (inf >> well >> id){
        well2id.insert(make_pair(well, id));
    }
}

void dump_counts(string& filename, 
    robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_ms_umis,
    vector<string>& ms_names,
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts){

    FILE* f = fopen(filename.c_str(), "w");
    
    // Print header line
    fprintf(f, "bc");
    for (vector<string>::iterator n = ms_names.begin(); n != ms_names.end(); ++n){
        if (n->length() > 0){
            fprintf(f, "\t%s", n->c_str());
        }
    }
    fprintf(f, "\n");

    for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x =
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){

        vector<pair<int, string> > v;
        bc_ms_counts.emplace(x->first, v);

        string bc_str = bc2str(x->first);
        fprintf(f, "%s", bc_str.c_str());
        for (int i = 0; i < x->second.size(); ++i){
            if (ms_names[i].length() > 0){
                // This MULTIseq barcode was intended to be present in the data set
                if (x->second[i] == NULL){
                    fprintf(f, "\t0");
                    bc_ms_counts[x->first].push_back(make_pair(0,
                        ms_names[i]));
                }
                else{
                    bc_ms_counts[x->first].push_back(make_pair(-x->second[i]->count(),
                        ms_names[i]));
                    fprintf(f, "\t%d", x->second[i]->count());
                }
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
     
int load_counts(string& filename,
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    vector<string>& labels){
    
    string line;
    ifstream inf(filename);
    bool first = true;
    bc cur_bc;

    vector<string> header;

    while (getline(inf, line)){
        istringstream splitter(line);
        string field;
        int idx = 0;
        while(getline(splitter, field, '\t')){
            if (first){
                header.push_back(field);
            }
            else{
                if (idx == 0){
                    str2bc(field.c_str(), cur_bc);
                    vector<pair<int, string> > v;
                    bc_ms_counts.emplace(cur_bc.to_ulong(), v);      
                }
                else{
                    int val = atoi(field.c_str());
                    if (val > 0){
                        bc_ms_counts[cur_bc.to_ulong()].push_back(make_pair(-val, header[idx]));
                    }
                }
            }
            ++idx;
        }
        first = false;
    }

    for (int i = 0; i < header.size()-1; ++i){
        labels.push_back(header[i+1]);
    }
    return header.size()-1;
}

void assign_ids(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    vector<string> labels,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned){

    int n_labels = labels.size();

    vector<vector<double> > obs;
    vector<unsigned long> bcs;

    vector<mixtureDist> dists;
    
    double off_target = 0.1;
    double covhi = 100;
    double covlow = 1;

    vector<string> names = {"multinomial", "poisson"};
    for (int i = 0; i < n_labels; ++i){
        vector<double> params;
        for (int j = 0; j < n_labels; ++j){
            if (j == i){
                params.push_back(1.0 - off_target);
            }
            else{
                params.push_back(off_target/(double)(n_labels-1));
            }
        }
        mixtureDist dist(names, vector<vector<double> >{ params, vector<double>{ covhi }});
        dist.name = labels[i];
        dist.set_num_inputs(0, n_labels);
        dists.push_back(dist);
        for (int j = i + 1; j < n_labels; ++j){
            vector<double> params2;
            for (int k = 0; k < n_labels; ++k){
                if (k == i || k == j){
                    params2.push_back((1.0-off_target)/2.0);
                }
                else{
                    params2.push_back(off_target/(n_labels-2));
                }
            }
            mixtureDist dist2(names, vector<vector<double> >{ params2, vector<double>{ covhi }});
            if (labels[i] < labels[j]){
                dist2.name = labels[i] + "+" + labels[j];
            }
            else{
                dist2.name = labels[j] + "+" + labels[i];
            }
            dist2.set_num_inputs(0, n_labels);
            dists.push_back(dist2);
        }
    }

    // Create a dist for garbage
    vector<double> params_noise;
    for (int i = 0; i < i < n_labels; ++i){
        params_noise.push_back(1.0/(double)n_labels);
    }
    mixtureDist dist_noise(names, vector<vector<double> >{ params_noise, vector<double>{ covlow }});
    dist_noise.name = "noise";
    dist_noise.set_num_inputs(0, n_labels);
    dists.push_back(dist_noise);

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        int tot = 0;
        vector<double> row;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
            row.push_back((double)-x->second[i].first);
        }
        if (tot > 0){
            row.push_back((double)tot);
            obs.push_back(row);
            bcs.push_back(x->first);
        }
    }

    mixtureModel mod(dists);
    mod.fit(obs);
    mod.print();
    exit(0);

}

void assign_ids3(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    int n_labels,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned){

    vector<vector<double> > obs;
    vector<unsigned long> bcs;

    // Define 3 distributions:
    // Single-label, high count
    // Double-label, high count
    // Noise (3+ labels), low count
    vector<mixtureDist> dists;
    
    double off_target = 0.1;
    vector<double> params_single = { 1.0-off_target, off_target/2, off_target/2 };
    vector<double> params_double = { (1.0-off_target)/2, (1.0-off_target)/2, off_target };
    vector<double> params_noise = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
    
    /* 
    vector<string> names = {"multinomial", "poisson"};
    mixtureDist singlet(names, vector<vector<double> >{ params_single, vector<double>{ 100 } });
    singlet.name = "singlet";
    mixtureDist doublet(names, vector<vector<double> >{ params_double, vector<double>{ 100 } });
    doublet.name = "doublet";
    mixtureDist noise(names, vector<vector<double> >{ params_noise, vector<double>{ 1 } });
    noise.name = "noise";
    singlet.set_num_inputs(0, 3);
    doublet.set_num_inputs(0, 3);
    noise.set_num_inputs(0, 3); 
    dists.push_back(singlet);
    dists.push_back(doublet);
    dists.push_back(noise);
    */
    mixtureDist singlet("multinomial", params_single);
    mixtureDist doublet("multinomial", params_double);
    mixtureDist noise("multinomial", params_noise);
    singlet.name = "singlet";
    doublet.name = "doublet";
    noise.name = "noise";
    singlet.set_num_inputs(0, 3);
    doublet.set_num_inputs(0, 3);
    noise.set_num_inputs(0, 3);
    dists.push_back(singlet);
    dists.push_back(doublet);
    dists.push_back(noise);

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
        
        sort(x->second.begin(), x->second.end());
        if (x->second[0].first == 0){
            // No counts
            continue;
        }    

        int tot = 0;
        int past2 = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
            if (i > 1){
                past2 += -x->second[i].first;
            }
        }
        vector<double> row = { (double)-x->second[0].first, (double)-x->second[1].first, (double)past2 };
        fprintf(stderr, "ROW %d %d %d\n", -x->second[0].first, -x->second[1].first, past2);
        obs.push_back(row);
        bcs.push_back(x->first);
    }
    
    mixtureModel mod(dists);
    mod.fit(obs);
    mod.print();
    exit(0);

}

double median(vector<int>& vals){
    sort(vals.begin(), vals.end());
    if (vals.size() == 0){
        return 0;
    }
    else if (vals.size() == 1){
        return vals[0];
    }
    if (vals.size() % 2 == 0){
        double v1 = (double)vals[vals.size()/2];
        double v2 = (double)vals[vals.size()/2-1];
        return (0.5*v1 + 0.5*v2);
    }
    else{
        return vals[(vals.size()-1)/2];
    }
}

void assign_ids2(robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >& bc_ms_counts,
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr,
    bool output_unassigned){

    vector<int> valid;
    vector<int> doub;
    vector<int> noise;
    
    double sum1 = 0.0;
    double tot1 = 0.0;
    double sum2 = 0.0;
    double tot2 = 0.0;
    double sum3 = 0.0;
    double tot3 = 0.0;

    for (robin_hood::unordered_map<unsigned long, vector<pair<int, string> > >::iterator x = 
        bc_ms_counts.begin(); x != bc_ms_counts.end(); ++x){
         
        // If we don't have a positive count for the second BC, add a pseudocount
        bool pseudo2 = false;
        // If we don't have a positive count for the third BC, add a pseudocount
        bool pseudo3 = false;

        sort(x->second.begin(), x->second.end());
        if (x->second.size() == 0){
            // No counts
            continue;
        }
        else if (x->second.size() == 1){
            pseudo2 = true;
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
            x->second.push_back(make_pair(-1, ""));
        }
        else if (x->second.size() == 2){
            pseudo3 = true;
            x->second.push_back(make_pair(-1, ""));
        }
        
        int tot = 0;
        for (int i = 0; i < x->second.size(); ++i){
            tot += -x->second[i].first;
        }
        double ll1 = 0;
        double ll2 = 0;
        double ll3 = 0;

        for (int ngroups = 1; ngroups <= 3; ++ngroups){
            double mean_true = 0;
            double tot_true = 0;
            for (int i = 0; i < ngroups; ++i){
                mean_true += -(double)x->second[i].first;
                tot_true++;
            }
            if (tot_true == 0){
                tot_true = 1;
            }
            mean_true /= tot_true;

            double mean_err = 0;
            double tot_err = 0;
            for (int i = ngroups; i < x->second.size(); ++i){
                mean_err += -(double)x->second[i].first;
                tot_err++;
            }
            if (tot_err == 0){
                tot_err = 1;
            }
            mean_err /= tot_err;
            if (mean_err == 0){
                // Prevent issues using Poisson dist
                mean_err = 1;
            }
            double ll = 0;
            for (int i = 0; i < ngroups; ++i){
                ll += dpois(-x->second[i].first, mean_true);
            }
            for (int i = ngroups; i < x->second.size(); ++i){
                ll += dpois(-x->second[i].first, mean_err);
            }
            
            if (ngroups == 1){
                ll1 = ll;
            }
            else if (ngroups == 2){
                ll2 = ll;
            }
            else if (ngroups == 3){
                ll3 = ll;
            }
        }
        
        if (pseudo2){
            tot--;
        }
        if (pseudo3){
            tot--;
        }
        
        double ll1b = dpois(tot, 17);
        double ll2b = dpois(tot, 20);
        double ll3b = dpois(tot, 13);

        double p_contam = 1.0-ppois(tot, 13);
        double frac1 = (double)-x->second[0].first / (double)tot;
        
        string name;
        double llr;
        bool store = false;
        
        double llmax = ll1;
        if (ll2 > llmax){
            llmax = ll2;
        }
        if (ll3 > llmax){
            llmax = ll3;
        }
        double probsum = pow(2, ll1-llmax) + pow(2, ll2-llmax) + pow(2, ll3-llmax);
        probsum = log2(probsum) + llmax;
        double p1 = pow(2, ll1-probsum);
        double p2 = pow(2, ll2-probsum);
        double p3 = pow(2, ll3-probsum);
        
        sum1 += p1*(double)tot;
        tot1 += p1;
        sum2 += p2*(double)tot;
        tot2 += p2;
        sum3 += p3*(double)tot;
        tot3 += p3;

        if (ll3 > ll1 && ll3 > ll2){
            // 3 or more == noise. Don't make an assignment.        
            //fprintf(stdout, "N\t%f\t%f\t%f\t%f\t%f\n", llr, p_contam, p1, p2, p3);
            fprintf(stdout, "N\t%f\t%f\t%f\t%f\t%f\t%f\n", ll1, ll2, ll3, ll1b, ll2b, ll3b);
            noise.push_back(tot);
        }
        else if (ll2 > ll1 && ll2 > ll3){
            if (!pseudo2){
                // Doublet.
                store = true;
                doub.push_back(tot);
                if (ll1 > ll3){
                    llr = ll2 - ll1;
                }
                else{
                    llr = ll2 - ll3;
                }
                //fprintf(stdout, "D\t%f\t%f\t%f\t%f\t%f\n", llr, p_contam, p1, p2, p3);
                fprintf(stdout, "D\t%f\t%f\t%f\t%f\t%f\t%f\n", ll1, ll2, ll3, ll1b, ll2b, ll3b);
                if (x->second[0].second < x->second[1].second){
                    name = x->second[0].second + "+" + x->second[1].second;
                }
                else{
                    name = x->second[1].second + "+" + x->second[0].second;
                }
            }
            else{
                // Noise
                noise.push_back(tot);
            }
        }
        else if (ll1 > ll2 && ll1 > ll3){
            store = true;
            llr = ll1 - ll2;
            name = x->second[0].second;
            valid.push_back(tot);
            //fprintf(stdout, "S\t%f\t%f\t%f\t%f\t%f\n", llr, p_contam, p1, p2, p3);
            fprintf(stdout, "S\t%f\t%f\t%f\t%f\t%f\t%f\n", ll1, ll2, ll3, ll1b, ll2b, ll3b);
        }
        else{
            //noise.push_back(tot);
        }
    }
    
    sum1 /= tot1;
    sum2 /= tot2;
    sum3 /= tot3;

    fprintf(stderr, "1 %f 2 %f 3 %f\n", sum1, sum2, sum3);

    double med_noise = median(noise);
    double l_noise = (med_noise - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "noise %f\n", l_noise);
    double med_valid = median(valid);
    double l_valid = (med_valid - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "valid %f\n", l_valid); 
    double med_doub = median(doub);
    double l_doub = (med_doub - 1.0/3.0)*(50.0/49.0);
    fprintf(stderr, "doub %f\n", l_doub); 
    exit(0);
    // Optional second pass through data: find mean successful count and mean unsuccessful
    // count; throw out identifications where successful count is more likely under the
    // data set-wide unsuccessful count distribution.

    // This seems unnecessary - now, LLRs reflect what we can say about the observed data; 
    // we don't need to care what other cells look like.
    
    bool second_pass = false;

    if (second_pass){
        double mean_valid = 0;
        double mean_noise = 0;

        double valid_weight = 1.0/(double)valid.size();
        for (int i = 0; i < valid.size(); ++i){
            mean_valid += valid_weight*(double)valid[i];    
        }
        double noise_weight = 1.0/(double)noise.size();
        for (int i = 0; i < noise.size(); ++i){
            mean_noise += noise_weight*(double)noise[i];
        }
        
        for (robin_hood::unordered_map<unsigned long, double>::iterator a = assn_llr.begin();
            a != assn_llr.end();){
            bool rm = false;
            double llr2 = 0.0;
            if (assn2.count(a->first) > 0){
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) -
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
                llr2 += dpois(-bc_ms_counts[a->first][1].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][1].first, mean_noise);
            }
            else{
                int count1 = -bc_ms_counts[a->first][0].first;
                int count2 = -bc_ms_counts[a->first][1].first;
                llr2 += dpois(-bc_ms_counts[a->first][0].first, mean_valid) - 
                    dpois(-bc_ms_counts[a->first][0].first, mean_noise);
            }
            
            double llr_new = a->second + llr2;
            if (llr_new < 0){
                // Eliminate
                assn1.erase(a->first);
                if (assn2.count(a->first) > 0){
                    assn2.erase(a->first);
                }
                rm = true;
            }
            else{
                a->second = llr_new;
            }

            if (rm){
                assn_llr.erase(a++);
            }
            else{
                a++;
            }
        }
    }
}

void write_assignments(string filename, 
    robin_hood::unordered_map<unsigned long, string>& assn1,
    robin_hood::unordered_map<unsigned long, string>& assn2,
    robin_hood::unordered_map<unsigned long, double>& assn_llr){
    
    FILE* f_out = fopen(filename.c_str(), "w");
    for (robin_hood::unordered_map<unsigned long, string>::iterator a = assn1.begin(); 
        a != assn1.end(); ++a){
        char type = 'S';
        string name = a->second;
        if (assn2.count(a->first) > 0){
            type = 'D';
            name = a->second + "+" + assn2[a->first];
        }
        string bc_str = bc2str(a->first);
        fprintf(f_out, "%s\t%s\t%c\t%f\n", bc_str.c_str(), name.c_str(), type, assn_llr[a->first]); 
    }
    fclose(f_out);
}

void check_missing_ms_wells(robin_hood::unordered_map<unsigned long, vector<umi_set*> >& bc_ms_umis,
    vector<string>& ms_wells,
    map<string, string>& well2name,
    string& outfilename){
    
    fprintf(stderr, "check missing MS wells\n");    
    // Count each well across all cells
    vector<pair<int, string> > counts;
    for (int i = 0; i < ms_wells.size(); ++i){
        counts.push_back(make_pair(0, ms_wells[i]));
    }
    for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x = 
        bc_ms_umis.begin(); x != bc_ms_umis.end(); ++x){
        for (int i = 0; i < x->second.size(); ++i){
            if (x->second[i] != NULL){
                counts[i].first += -x->second[i]->count();
            }
        }
    }
    sort(counts.begin(), counts.end());
    
    int minrank_missing = -1;
    int toprank_present = -1;
    FILE* outf = fopen(outfilename.c_str(), "w");
    for (int i = 0; i < counts.size(); ++i){
        if (well2name.count(counts[i].second) == 0 && well2name[counts[i].second] != ""){
            if (minrank_missing == -1){
                minrank_missing = i;
            }
        }
        else{
            toprank_present = i;
        }
        fprintf(outf, "%d\t%s\t%s\n", -counts[i].first, counts[i].second.c_str(),
            well2name[counts[i].second].c_str());
    }
    fclose(outf);
    if (minrank_missing < toprank_present){
        fprintf(stderr, "WARNING: one or more un-named wells has higher counts than one or \
more named/expected wells. You may have mis-specified well-to-ID mappings in provided \
metadata. Please see the output well counts file to evaluate.\n");
    }

}

int main(int argc, char *argv[]) {    
    
    // Define long-form program options 
    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"read1", required_argument, 0, '1'},
       {"read2", required_argument, 0, '2'},
       {"mapping", required_argument, 0, 'm'},
       {"whitelist", required_argument, 0, 'w'},
       {"barcodes", required_argument, 0, 'b'},
       {"cell_barcodes", required_argument, 0, 'B'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"mismatches", required_argument, 0, 'M'},
       {"output_unassigned", no_argument, 0, 'u'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string output_prefix = "";
    vector<string> read1fn;
    vector<string> read2fn;
    string mapfn = "";
    string wlfn = "";
    string cell_barcodesfn;
    int mismatches = -1;
    bool output_unassigned = false;

    fprintf(stderr, "%s\n", STRINGIZE(PROJ_ROOT));
    string pr = STRINGIZE(PROJ_ROOT);
    if (pr[pr.length()-1] != '/'){
        pr += "/";
    }
    string barcodesfn = pr + "data/multiseq_indices.txt";
    double doublet_rate = 0.1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:1:2:m:M:w:b:B:D:uh", 
        long_options, &option_index )) != -1){
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
            case '1':
                read1fn.push_back(optarg);
               break;
            case '2':
                read2fn.push_back(optarg);
                break;
            case 'm':
                mapfn = optarg;
                break;
            case 'M':
                mismatches = atoi(optarg);
                break;
            case 'w':
                wlfn = optarg;
                break;
            case 'b':
                barcodesfn = optarg;
                break;
            case 'B':
                cell_barcodesfn = optarg;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'u':
                output_unassigned = true;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: --output_prefix/-o is required\n");
        exit(1);
    }
    if (mismatches < -1){
        fprintf(stderr, "ERROR: mismatches must be >= 0, or -1 to allow all matches\n");
        exit(1);
    }
    
    // Data structure to store counts
    robin_hood::unordered_map<unsigned long, vector<pair<int, string> > > bc_ms_counts;
    
    // How many possible MULTI-seq labels are there?
    int n_labels;
    vector<string> labels;

    string countsfilename = output_prefix + ".counts";
    if (file_exists(countsfilename)){
        fprintf(stderr, "Loading previously-computed counts (delete file to avoid this \
next time)...\n");
        n_labels = load_counts(countsfilename, bc_ms_counts, labels);
    }
    else{
        if (read1fn.size() == 0 || read2fn.size() == 0){    
            fprintf(stderr, "ERROR: input FASTQs for reads (-1 and -2) are required.\n");
            exit(1);
        }
        else if (read1fn.size() != read2fn.size()){
            fprintf(stderr, "ERROR: unequal number of forward (-1) and reverse (-2) reads provided\n");
            exit(1);
        }
        else if (wlfn == ""){
            fprintf(stderr, "ERROR: cell barcode whitelist / -w is required\n");
            exit(1);
        }
        
        // Load MULTIseq barcode -> well ID mapping
        map<unsigned long, string> bc2well;
        vector<string> bclist;
        int ms_bc_len = load_well2bc(barcodesfn, bc2well, bclist);
        
        // Load cell barcodes, if provided
        set<unsigned long> cell_barcodes;
        bool has_cell_barcodes;
        if (cell_barcodesfn != ""){
            has_cell_barcodes = true;
            parse_barcode_file(cell_barcodesfn, cell_barcodes);
        }    

        // Load MULTIseq well -> unique identifier mapping
        map<string, string> well2id;
        if (mapfn == ""){
            fprintf(stderr, "WARNING: no MULTIseq well to identifier mapping provided.\n");
            fprintf(stderr, "Using well names as identifiers.\n");    
            for (map<unsigned long, string>::iterator bw = bc2well.begin(); 
                bw != bc2well.end(); ++bw){
                if (well2id.count(bw->second) == 0){
                    well2id.insert(make_pair(bw->second, bw->second));
                }
            }
        }
        else{
            load_well_mapping(mapfn, well2id);
            // Make sure wells specified here have associated valid barcodes.
            set<string> wellnames;
            for (map<unsigned long, string>::iterator bw = bc2well.begin(); bw != 
                bc2well.end(); ++bw){
                wellnames.insert(bw->second);
            }
            for (map<string, string>::iterator wi = well2id.begin(); wi != well2id.end();
                ++wi){
                if (wellnames.find(wi->first) == wellnames.end()){
                    fprintf(stderr, "ERROR: well -> ID mapping file contains a well named %s\n",
                        wi->first.c_str());
                    fprintf(stderr, "This well ID is not present in the MULTIseq barcode -> \
well mappings\n");
                    exit(1);
                }
            }
            // If the user has omitted one or more wells that turn up often in the data, we will
            // alert the user about this later.
        }
        n_labels = well2id.size();

        // Initiate MULTIseq barcode whitelist
        seq_fuzzy_match msmatch(bclist, mismatches, true, true);
        
        vector<unsigned long> ms_bc;
        bc cur_bc;
        for (int i = 0; i < bclist.size(); ++i){
            str2bc(bclist[i].c_str(), cur_bc);
            ms_bc.push_back(cur_bc.to_ulong());
        }
        /*
        // Determine appropriate k-mer length for fuzzy matching
        // k must be < (bc_length + 1)/2
        int ms_k;
        if (ms_bc_len % 2 == 0){
            // Even barcode length.
            ms_k = ms_bc_len/2;
        }
        else{
            ms_k = (ms_bc_len-1)/2;
        }
        */
        //bc_whitelist wl_ms(bclist, ms_bc_len, ms_k);
        
        map<unsigned long, int> ms_bc2idx;
        vector<string> ms_names;
        vector<string> ms_wells;
        int ix = 0;
        for (map<unsigned long, string>::iterator bw = bc2well.begin(); 
            bw != bc2well.end(); ++bw){
            
            ms_bc2idx.insert(make_pair(bw->first, ix));
            ++ix;
            ms_wells.push_back(bw->second);
            if (well2id.count(bw->second) > 0){
                ms_names.push_back(well2id[bw->second]);
                labels.push_back(well2id[bw->second]);
            }
            else{
                ms_names.push_back("");
            }
        }
        
        // Data structure to store MULTIseq barcode counts per cell barcode
        // counts come from UMIs
        robin_hood::unordered_map<unsigned long, vector<umi_set*> > bc_ms_umis;
        
        // Set up object that will scan each pair of read files
        bc_scanner scanner;
        scanner.init_multiseq_v3(wlfn);
        
        char ms_bc_buf[ms_bc_len+1];

        for (int i = 0; i < read1fn.size(); ++i){
            
            // Initiate object to read through FASTQs and find cell barcodes.
            scanner.add_reads(read1fn[i], read2fn[i]);
            
            while (scanner.next()){
                
                // Skip reads if we have a cell barcode whitelist and the current
                // barcode isn't in it
                if (has_cell_barcodes && 
                    cell_barcodes.find(scanner.barcode) == cell_barcodes.end()){
                    continue;
                }
                
                // At this point, there's a valid cell barcode.
                // scanner.barcode_read holds R1
                // scanner.read_f holds R2
                // scanner.umi holds UMI 
                if (scanner.has_umi){
                    // MULTIseq barcodes are the first 8 bp (or however long supplied barcodes are)
                    // of R2, forward orientation
                    unsigned long ms_ul;
                    
                    strncpy(&ms_bc_buf[0], scanner.read_f, ms_bc_len);
                    ms_bc_buf[ms_bc_len] = '\0';            
                    int idx = msmatch.match(&ms_bc_buf[0]);
                    if (idx != -1){
                        ms_ul = ms_bc[idx];
                    //if (wl_ms.lookup(scanner.read_f, false, ms_ul)){
                        // Store UMI
                        umi this_umi(scanner.umi, scanner.umi_len);    
                        if (bc_ms_umis.count(scanner.barcode) == 0){
                            vector<umi_set*> v(bclist.size(), NULL);
                            bc_ms_umis.emplace(scanner.barcode, v);
                            //map<unsigned long, umi_set> m;
                            //bc_ms_umis.emplace(scanner.barcode, m);
                        }
                        umi_set* ptr = bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]];
                        if (ptr == NULL){
                            // Initialize.
                            bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]] = new umi_set(scanner.umi_len);
                            ptr = bc_ms_umis[scanner.barcode][ms_bc2idx[ms_ul]];
                        }
                        ptr->add(this_umi); 
                    }
                }
            }
        }

        // Write counts to disk
        string countsfn = output_prefix + ".counts";
        dump_counts(countsfn, bc_ms_umis, ms_names, bc_ms_counts);
       
        // Check to see that the desired MULTIseq barcodes are the most common ones.
        // Warn the user if unexpected ones are more common. 
        string wellcountsfn = output_prefix + ".wells";
        check_missing_ms_wells(bc_ms_umis, ms_wells, well2id, wellcountsfn);        

        // Free stuff
        for (robin_hood::unordered_map<unsigned long, vector<umi_set*> >::iterator x = 
            bc_ms_umis.begin(); x != bc_ms_umis.end(); ){
            for (int i = 0; i < x->second.size(); ++i){
                if (x->second[i] != NULL){
                    delete x->second[i];
                }
            }
            bc_ms_umis.erase(x++);
        }
    }

    // Now we can do the actual demultiplexing
    robin_hood::unordered_map<unsigned long, string> assn1;
    robin_hood::unordered_map<unsigned long, string> assn2;
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    assign_ids2(bc_ms_counts, assn1, assn2, assn_llr, output_unassigned);
    
    string assnfilename = output_prefix + ".assignments";
    write_assignments(assnfilename, assn1, assn2, assn_llr);

    return 0;  
}
