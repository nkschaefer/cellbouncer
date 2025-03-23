#include <getopt.h>
#include <argp.h>
#include <string>
#include <algorithm>
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
#include <optimML/functions.h>
#include <optimML/multivar_ml.h>
#include <mixtureDist/functions.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;


/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "doublet_dragon [OUTPUT_PREFIX] [FILES]\n");
    fprintf(stderr, "[OUTPUT_PREFIX] is the base file name for output files.\n");
    fprintf(stderr, "[OUTPUT_PREFIX].dd.all and [OUTPUT_PREFIX].dd.indv will be created.\n");
    fprintf(stderr, "[FILES] is a list of one or more .assignments files corresponding to the same library.\n");
    fprintf(stderr, "-----\n");
    fprintf(stderr, "Given one or more .assignments files from the same library, estimates\n");
    fprintf(stderr, "the likeliest underlying doublet rate (via maximum likelihood).\n");
    fprintf(stderr, "Also infers the frequency of each individual type in the data set.\n");
    fprintf(stderr, "If the same names are encountered in multiple assignments files, it is\n");
    fprintf(stderr, "assumed they correspond to the same individual/ID. If this is not the case\n");
    fprintf(stderr, "(for example, if you have clustered mitochondria from multiple species in the\n");
    fprintf(stderr, "same library, and they both have numeric IDs), be sure to change the IDs in one\n");
    fprintf(stderr, "of the files (in this example, create .ids files for both).\n");
    fprintf(stderr, "\nEXAMPLE 1:\n");
    fprintf(stderr, "If you pooled data from three species, separated these data using\n");
    fprintf(stderr, "demux_species, inferred mitochondrial haplotypes via demux_mt,\n");
    fprintf(stderr, "demultiplexed individuals from a VCF using demux_vcf, and processed\n");
    fprintf(stderr, "MULTIseq data on the same library using demux_multiseq, then you can\n");
    fprintf(stderr, "run doublet_dragon on all these output files to estimate the underlying\n");
    fprintf(stderr, "global doublet rate.\n");
    fprintf(stderr, "\nEXAMPLE 2:\n");
    fprintf(stderr, "If you identified individuals in the same data set multiple ways (i.e. via\n");
    fprintf(stderr, "cell-hashing, SNP-based demultiplexing, and mitochondrial clustering), then\n");
    fprintf(stderr, "make sure all of these .assignments files have the same names for the individuals\n");
    fprintf(stderr, "and run doublet_dragon to get a global doublet rate and find over/underrepresentation\n");
    fprintf(stderr, "of individual IDs in all of the .assignments files.\n");
    fprintf(stderr, "\nEXAMPLE 3:\n");
    fprintf(stderr, "You have just run demux_vcf on a single data set and want to find the most likely\n");
    fprintf(stderr, "doublet rate, including self-self doublets.\n");
    exit(code);
}

/**
 * Returns total number of cells
 */
int parse_assignments(string filename, 
    map<string, int>& counts_singlet, 
    map<string, int>& counts_doublet,
    map<string, double>& weights){
    
    ifstream inf(filename.c_str());
    string bc_str;
    string id;
    string type;
    double llr;
    
    map<string, int> id_count;
    
    int ntot = 0;
    
    string line;
    while(getline(inf, line)){
        istringstream splitter(line);
        string field;
        int idx = 0;
        while(getline(splitter, field, '\t')){
            if (idx == 0){
                bc_str = field;
            }
            else if (idx == 1){
                id = field;
            }
            else if (idx == 2){
                type = field;
            }
            else if (idx == 3){
                llr = atof(field.c_str());
            }
            else{
                fprintf(stderr, "ERROR: unexpected number of fields in %s.\n", filename.c_str());
                exit(1);
            }
            ++idx;
        }

        ntot++;
        if (id_count.count(id) == 0){
            id_count.insert(make_pair(id, 0));
            weights.insert(make_pair(id, 0.0));
        }    
        
        id_count[id]++;
        weights[id] += llr;
        
        if (type == "S"){
            if (counts_singlet.count(id) == 0){
                counts_singlet.insert(make_pair(id, 0));
            }
            counts_singlet[id]++;
        }
        else if (type == "D"){
            if (counts_doublet.count(id) == 0){
                counts_doublet.insert(make_pair(id, 0));
            }
            counts_doublet[id]++;
        }
        else if (type == "M"){
            // Skip these for now.
        }
        else{
            fprintf(stderr, "ERROR: unexpected value encountered in type field: %s\n", type.c_str());
            fprintf(stderr, "doublet_dragon cannot accept assignments files for sgRNA capture data.\n");
            fprintf(stderr, "Allowed values for this field are S/D/M.\n");
            exit(1);
        }
    }

    double wsum = 0.0;
    for (map<string, double>::iterator w = weights.begin(); w != weights.end(); ++w){
        w->second = w->second / (double)id_count[w->first];
        wsum += w->second;
    }
    for (map<string, double>::iterator w = weights.begin(); w != weights.end(); ++w){
        w->second /= wsum;
    }

    return ntot;
}

double loglik(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double d = params[0];
    double n = (double)data_i.at("n");
    double k = (double)data_i.at("k");
    int p_idx = data_i.at("p_idx");
    int p_idx2 = data_i.at("p_idx2");
    //int grp_start = data_i.at("grp_start");
    //int grp_end = data_i.at("grp_end");
    
    /*
    double grptot = 0;
    for (int i = grp_start; i <= grp_end; ++i){
        grptot += expit(params[i]);
    }
    double p1 = expit(params[p_idx]) / grptot;
    */
    double p1 = params[p_idx];

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d)*p1 + d*p1 * p1;
    }
    else{
        // Doublet
        //double p2 = expit(params[p_idx2]) / grptot;
        double p2 = params[p_idx2];
        // Account for two ways of sampling this doublet
        p = 2*d*p1*p2;
    }
    return logbinom(n, k, p);
}

void dloglik_dx(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double d = params[0];
    double n = (double)data_i.at("n");
    double k = (double)data_i.at("k");
    int p_idx = data_i.at("p_idx");
    int p_idx2 = data_i.at("p_idx2");
    //int grp_start = data_i.at("grp_start");
    //int grp_end = data_i.at("grp_end");
    /*
    double grptot = 0;
    for (int i = grp_start; i <= grp_end; ++i){
        grptot += expit(params[i]);
    }
    double p1 = expit(params[p_idx]) / grptot;
    */
    double p1 = params[p_idx];
    double p2;

    double p;
    if (p_idx2 == -1){
        // Singlet
        p = (1.0 - d)*p1 + d*p1 * p1;
    }
    else{
        // Doublet
        //double p2 = expit(params[p_idx2])/grptot;
        p2 = params[p_idx2];
        p = 2.0*d*p1*p2;
    }

    double dy_dp = (k - n*p)/(p - p*p);
    /*
    double e_negx1 = exp(-params[p_idx]);
    double e_negx1_p1_2 = pow(e_negx1 + 1, 2);
    double e_negx1_p1_3 = e_negx1_p1_2 * (e_negx1 + 1);
    double der_comp_1 = e_negx1 / (e_negx1_p1_2 * grptot) - e_negx1 / (e_negx1_p1_3 *grptot *grptot);
    */
    if (p_idx2 == -1){
        //double a = expit(params[p_idx])/grptot;
        //double a2 = a*a;
        //results[0] += dy_dp * (-a + a2);
        //results[p_idx] += dy_dp * (1.0 - d + 2*d*a) * der_comp_1;
        results[0] += dy_dp * (-p1 + p1*p1);
        results[p_idx] += dy_dp * (2*d*p1 + 1 - d);
    }
    else{
        /*
        double a = expit(params[p_idx])/grptot;
        double b = expit(params[p_idx2])/grptot;
        results[0] += 2.0 * dy_dp * (a*b);
        results[p_idx] += 2.0 * dy_dp * (d*b) * der_comp_1;

        double e_negx2 = exp(-params[p_idx2]);
        double e_negx2_p1_2 = pow(e_negx2 + 1, 2);
        double e_negx2_p1_3 = e_negx2_p1_2 * (e_negx2 + 1);
        double der_comp_2 = e_negx2 / (e_negx2_p1_2 * grptot) - e_negx2/ ( e_negx2_p1_3 * grptot * grptot);

        results[p_idx2] += 2.0 * dy_dp * (d*a) * der_comp_2;
        */
        results[0] += dy_dp*(2*p1*p2);
        results[p_idx] += dy_dp*(2*d*p2);
        results[p_idx2] += dy_dp*(2*d*p1);
    }
}

int main(int argc, char *argv[]) {    
    
    if (argc < 3){
        help(1);
    } 
    
    string outbase = argv[1];
    if (is_dir(outbase)){
        fprintf(stderr, "ERROR: first argument, %s, is a directory but should be \
a file name prefix.\n", outbase.c_str());
        exit(1);
    } 
    // Make sure the user didn't accidentally pass an output file as first argument
    if (outbase.length() > 12 && outbase.substr(outbase.length()-12,12) == ".assignments"){
        fprintf(stderr, "ERROR: first argument, [output_prefix], is an .assignments file.\n");
        exit(1);
    }
    else if (file_exists(outbase)){
        fprintf(stderr, "ERROR: first argument, [output_prefix], is a file that exists.\n");
        exit(1);
    }
    else{
        string out_test = outbase + ".assignments";
        if (file_exists(out_test)){
            fprintf(stderr, "ERROR: first argument, [output_prefix], is a preexisting output_prefix\n");
            fprintf(stderr, "  (to be loaded?)\n");
            exit(1);
        }
    }

    vector<int> n;
    vector<int> k;
    vector<int> p_idx;
    vector<int> p_idx2;
    vector<int> grp_idx;

    vector<double> weights;
    
    vector<double> params;
    params.push_back(0.01); // doublet rate
    
    map<pair<int, int>, string> idx2name;
    map<string, pair<int, int> > name2idx;
    
    vector<vector<string> > grp_names;
    vector<vector<double> > grp_props;
    // To set initial values: above will be sums and below will be counts - 
    // take means
    vector<vector<double> > grp_props_denom;

    map<string, map<string, int> > file2singlets;
    map<string, map<string, int> > file2doublets;
    map<string, int> file2grp;

    map<string, int> file_tot;
    
    for (int i = 2; i < argc; ++i){
        map<string, int> counts_singlet;
        map<string, int> counts_doublet;
        map<string, double> weights_this;
        
        // Ensure this is an .assignments file.
        string fn = argv[i];
        if (fn.length() <= 12 || fn.substr(fn.length()-12, 12) != ".assignments"){
            // Assume it's a base file name instead
            string fn2 = fn + ".assignments";
            if (!file_exists(fn2)){
                fprintf(stderr, "ERROR: file %s does not appear to be an .assignments file or\n", fn.c_str());
                fprintf(stderr, "  an output_prefix with an assignments file.\n");
                exit(1);
            }
            else{
                fn = fn2;
            }
        }
        int totcells = parse_assignments(fn, counts_singlet, counts_doublet, weights_this);
        if (counts_doublet.size() == 0){
            fprintf(stderr, "ERROR: no doublet identifications in %s\n", argv[i]);
            exit(1);
        }

        file2singlets.insert(make_pair(fn, counts_singlet));
        file2doublets.insert(make_pair(fn, counts_doublet));
        file_tot.insert(make_pair(fn, totcells));
         
        int tot_singlet = 0;
        for (map<string, int>::iterator cs = counts_singlet.begin();
            cs != counts_singlet.end(); ++cs){
            tot_singlet += cs->second;
        }
        
        int grp_idx_this = -1;

        for (map<string, int>::iterator cs = counts_singlet.begin(); 
            cs != counts_singlet.end(); ++cs){
            int p_idx_this;
            if (name2idx.count(cs->first) > 0){
                
                // This individual already exists as part of a set.

                pair<int, int> grp_and_idx = name2idx[cs->first];

                if (grp_idx_this != -1 && grp_and_idx.first != grp_idx_this){
                    fprintf(stderr, "ERROR: %s contains individuals from two different sets.\n",
                        fn.c_str());
                    fprintf(stderr, "Set 1:\n");
                    for (vector<string>::iterator x = grp_names[grp_idx_this].begin(); x !=
                        grp_names[grp_idx_this].end(); ++x){
                        fprintf(stderr, "%s\n", x->c_str());
                    }
                    fprintf(stderr, "Set 2:\n");
                    for (vector<string>::iterator x = grp_names[grp_and_idx.first].begin(); 
                        x != grp_names[grp_and_idx.first].end(); ++x){
                        fprintf(stderr, "%s\n", x->c_str());
                    }
                    fprintf(stderr, "Overlapping name: %s\n", cs->first.c_str());
                    exit(1);
                }

                // If name matches an already-existing name, treat them as the
                // same variable.
                grp_idx_this = grp_and_idx.first;
                p_idx_this = grp_and_idx.second;

                // Make the initial value the mean initial value across all data sets
                grp_props[grp_idx_this][p_idx_this] += ((double)cs->second / (double)tot_singlet);
                grp_props_denom[grp_idx_this][p_idx_this]++;

            }
            else if (grp_idx_this != -1){
                // This file already corresponds to a group, but this individual has not yet been
                // added to it.
                p_idx_this = grp_props[grp_idx_this].size();
                grp_props[grp_idx_this].push_back((double)cs->second / (double)tot_singlet);
                grp_props_denom[grp_idx_this].push_back(1.0);
                grp_names[grp_idx_this].push_back(cs->first);
                
                name2idx.insert(make_pair(cs->first, make_pair(grp_idx_this, p_idx_this)));
                idx2name.insert(make_pair(make_pair(grp_idx_this, p_idx_this), cs->first));
            }
            else{
                // A new group needs to be created.
                grp_idx_this = grp_props.size();
                p_idx_this = 0;
                vector<double> grp_prop_new = { (double)cs->second / (double)tot_singlet };
                vector<double> grp_denom_new = { 1.0 };
                vector<string> grp_name_new = { cs->first };
                
                grp_props.push_back(grp_prop_new);
                grp_props_denom.push_back(grp_denom_new);
                grp_names.push_back(grp_name_new);
                
                name2idx.insert(make_pair(cs->first, make_pair(grp_idx_this, p_idx_this)));
                idx2name.insert(make_pair(make_pair(grp_idx_this, p_idx_this), cs->first));
            }
                    
            n.push_back(totcells);
            k.push_back(cs->second);
            p_idx.push_back(p_idx_this);
            p_idx2.push_back(-1);
            grp_idx.push_back(grp_idx_this);
            weights.push_back(weights_this[cs->first]);
        }
        for (map<string, int>::iterator cd = counts_doublet.begin();
            cd != counts_doublet.end(); ++cd){
            size_t pos = cd->first.find("+");
            if (pos == string::npos){
                fprintf(stderr, "ERROR: unable to identify components of doublet ID %s\n", cd->first.c_str());
                exit(1);
            }
            string id1 = cd->first.substr(0, pos);
            string id2 = cd->first.substr(pos+1, cd->first.length()-pos-1);
            if (name2idx.count(id1) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s (%s)\n", id1.c_str(), argv[i]);
                for (map<string, pair<int, int> >::iterator ni = name2idx.begin(); ni != name2idx.end(); ++ni){
                    fprintf(stderr, "%s\t%d %d\n", ni->first.c_str(), ni->second.first, ni->second.second);
                }
                exit(1);
            }
            if (name2idx.count(id2) == 0){
                fprintf(stderr, "ERROR: no singlets for ID %s (%s)\n", id2.c_str(), argv[i]);
                for (map<string, pair<int, int> >::iterator ni = name2idx.begin(); ni != name2idx.end(); ++ni){
                    fprintf(stderr, "%s\t%d %d\n", ni->first.c_str(), ni->second.first, ni->second.second);
                }
                exit(1);
            }
            
            pair<int, int> grp_and_idx1 = name2idx[id1];
            pair<int, int> grp_and_idx2 = name2idx[id2];

            if (grp_and_idx1.first != grp_and_idx2.first){
                fprintf(stderr, "ERROR: doublet %s corresponds to singlets from two different sets\n",
                    cd->first.c_str());
                exit(1);
            }
            else if (grp_and_idx1.first != grp_idx_this || grp_and_idx2.first != grp_idx_this){
                fprintf(stderr, "ERROR: file %s has individuals from multiple sets\n", fn.c_str());
                exit(1);
            }

            n.push_back(totcells);
            k.push_back(cd->second);
            p_idx.push_back(grp_and_idx1.second);
            p_idx2.push_back(grp_and_idx2.second);
            grp_idx.push_back(grp_and_idx1.first);
            weights.push_back(weights_this[cd->first]);
        }
        file2grp.insert(make_pair(fn, grp_idx_this));
    }
    
    // Transform starting props into means, and then normalize per group.
    for (int i = 0; i < grp_props.size(); ++i){
        double grptot = 0.0;
        for (int j = 0; j < grp_props[i].size(); ++j){
            grp_props[i][j] /= grp_props_denom[i][j];
            grptot += grp_props[i][j];
        }
        for (int j = 0; j < grp_props[i].size(); ++j){
            grp_props[i][j] /= grptot;
        }
    }
    
    // Transform group indices into overall indices in the final parameter vector
    vector<int> grpstart;
    int idx_global = 1; // account for first being doublet rate
    for (int i = 0; i < grp_props.size(); ++i){
        grpstart.push_back(idx_global);
        idx_global += grp_props[i].size();
    }
    for (int i = 0; i < p_idx.size(); ++i){
        int gi = grp_idx[i];
        // Current index is index into group. Convert to global index
        p_idx[i] += grpstart[gi];
        if (p_idx2[i] != -1){
            p_idx2[i] += grpstart[gi];
        }
    }
   
    optimML::multivar_ml_solver solver(params, loglik, dloglik_dx);
    // Constrain doublet rate
    solver.constrain_01(0);
    // Add in parameter groups
    for (int i = 0; i < grp_props.size(); ++i){
        solver.add_param_grp(grp_props[i]);
    } 
    
    solver.add_data("n", n);
    solver.add_data("k", k);
    solver.add_data("p_idx", p_idx);
    solver.add_data("p_idx2", p_idx2);

    solver.add_weights(weights);
    
    solver.solve();
    vector<double> results_max = solver.results;
    double llmax = solver.log_likelihood;
    for (double p = 0.05; p < 1; p += 0.05){
        solver.set_param(0, p);
        solver.solve();
        if (solver.log_likelihood > llmax){
            results_max = solver.results;
            llmax = solver.log_likelihood;
        } 
    }

    if (outbase[outbase.length()-1] == '.'){
        outbase = outbase.substr(0, outbase.length()-1);
    }
    string fn_all = outbase + ".dd.all";
    string fn_indv = outbase + ".dd.indv";
    FILE* out_all = fopen(fn_all.c_str(), "w");
    FILE* out_indv = fopen(fn_indv.c_str(), "w");

    double doublet_rate = results_max[0];
    map<string, double> indv_freq;
    double ifsum = 0.0;
    for (map<string, pair<int, int> >::iterator ni = name2idx.begin(); ni != name2idx.end(); ++ni){
        int solver_idx = grpstart[ni->second.first] + ni->second.second;
        indv_freq.insert(make_pair(ni->first, results_max[solver_idx]));
        fprintf(stderr, "%s\t%f\n", ni->first.c_str(), indv_freq[ni->first]);
    }
    for (int gi = 0; gi < grp_names.size(); ++gi){
        // Sort alphabetically
        vector<pair<string, int> > grpnamesort;
        for (int ii = 0; ii < grp_names[gi].size(); ++ii){
            grpnamesort.push_back(make_pair(grp_names[gi][ii], ii));
        }
        sort(grpnamesort.begin(), grpnamesort.end());
        for (int j = 0; j < grpnamesort.size(); ++j){
            fprintf(out_all, "group_%d\t%s\t%f\n", gi, grpnamesort[j].first.c_str(),
                indv_freq[grpnamesort[j].first]);
        }
    }
    
    fprintf(stderr, "Global doublet rate: %f\n", doublet_rate);
    fprintf(out_all, "all\tdoublet_rate\t%f\n", doublet_rate);
    fprintf(stderr, "Multinomial log likelihoods:\n");
    //fprintf(stdout, "data_set\ttype\tindv\tcount\texpected\n");    
    
    for (map<string, int>::iterator ft = file_tot.begin(); ft != file_tot.end(); ++ft){

        string fntrunc = ft->first.substr(0, ft->first.length()-12);
        double tot = ft->second;
        vector<double> x_all;
        vector<double> p_all;
        double ptot = 0.0;
        double chisq = 0.0;
        int chisq_df = 0;

        int grp_idx = file2grp[ft->first];
        vector<string> singlets_vec;

        // Go through singlets
        for (vector<string>::iterator s = grp_names[grp_idx].begin(); 
            s != grp_names[grp_idx].end(); ++s){
            
            singlets_vec.push_back(*s);
            double p1 = indv_freq[*s];
            double p = (1.0-doublet_rate)*p1 + doublet_rate*p1*p1;
            double expec = p*tot;
            p_all.push_back(p);
            int num = 0;
            if (file2singlets[ft->first].count(*s) > 0){
                num = file2singlets[ft->first][*s];
            }
            x_all.push_back((double)num);
            double p_indv = pbinom(tot, (double)num, p);
            chisq += pow((double)num - expec, 2)/expec;
            chisq_df++;
            fprintf(out_indv, "%s\tS\t%s\t%d\t%f\t%f\n", fntrunc.c_str(), s->c_str(),
                num, expec, p_indv);
        }
        
        // Go through doublets
        for (int i = 0; i < singlets_vec.size()-1; ++i){
            double p1 = indv_freq[singlets_vec[i]];
            for (int j = i +1; j < singlets_vec.size(); ++j){
                double p2 = indv_freq[singlets_vec[j]];
                double p = 2.0*doublet_rate*p1*p2;
                double expected = p*tot;
                p_all.push_back(p);
                string dname = singlets_vec[i] + "+" + singlets_vec[j];
                int num = 0;
                if (file2doublets[ft->first].count(dname) > 0){
                    num = file2doublets[ft->first][dname];
                }
                x_all.push_back((double)num);
                double expec = p*tot;
                double p_indv = pbinom(tot, (double)num, p);
                chisq += pow((double)num - expec, 2)/expec;
                chisq_df++;
                fprintf(out_indv, "%s\tD\t%s\t%d\t%f\t%f\n", fntrunc.c_str(), dname.c_str(),
                    num, expec, p_indv);
            }
        }
        double loglik = dmultinom(x_all, p_all);
        chisq_df--;
        // With this many observations, chi-squared p values are often all significant - ignore
        double chisq_p = pchisq(chisq, chisq_df);
        fprintf(stderr, "%s\t%f\n", fntrunc.c_str(), loglik);
        fprintf(out_all, "data_set\t%s\t%f\n", fntrunc.c_str(), loglik);
    }
    
    fclose(out_all);
    fclose(out_indv);
    return 0;
}
