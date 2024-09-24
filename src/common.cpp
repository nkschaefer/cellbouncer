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
#include <sys/stat.h>
#include <mixtureDist/functions.h>
#include <htswrapper/bc.h>
#include <htswrapper/gzreader.h>
#include <optimML/multivar_ml.h>

/**
 * Contains functions used by more than one program in this
 * repository.
 */

using std::cout;
using std::endl;
using namespace std;

/**
 * Print data to help screens describing the --libname argument.
 */
void print_libname_help(){
    fprintf(stderr, "   --libname -n By default, cell barcodes will be printed in output files\n");
    fprintf(stderr, "       as DNA sequences only. If you plan to combine multiple data sets for\n");
    fprintf(stderr, "       analysis, this can lead to barcode collisions, where the same barcode sequence\n");
    fprintf(stderr, "       appears by chance in multiple data sets, but representing different cells.\n");
    fprintf(stderr, "       To account for this, a unique name must be appended to barcode sequences.\n");
    fprintf(stderr, "       If you provide a unique name here, it will be appended to cell barcode\n");
    fprintf(stderr, "       sequences, separated by a hyphen (-) (unless specified otherwise; see below).\n");
    fprintf(stderr, "       By default, this mimics the behavior of the batch_categories argument to\n");
    fprintf(stderr, "       anndata.concatenate() in scanpy.\n");
    fprintf(stderr, "       If you combine data using anndata.concatenate() without batch_categories, then\n");
    fprintf(stderr, "       this argument should be set to the 1-based numeric index of this data set in the\n");
    fprintf(stderr, "       list of all data sets you plan to combine.\n");
    fprintf(stderr, "       If you are combining multiple runs using cellranger aggr, you should set this\n");
    fprintf(stderr, "       to the 1-based numeric index of this data set among all data sets you will\n");
    fprintf(stderr, "       combine and omit the --cellranger/-C argument (below).\n");
    fprintf(stderr, "       If you are combining data sets using Seurat merge(), this should be set to what you\n");
    fprintf(stderr, "       will pass to the add_cell_ids argument, and you should specify Seurat format (below).\n");
    fprintf(stderr, "       If you want to use an underscore instead of a hyphen, see below.\n");
    fprintf(stderr, "   --cellranger -C Set this flag to append a -1 to the end of all output barcodes, like\n");
    fprintf(stderr ,"       10X Genomics CellRanger.\n");
    fprintf(stderr, "   --seurat -S Set this flag to append --libname to the beginning of barcode sequences, with an\n");
    fprintf(stderr, "       underscore separator.\n");
    fprintf(stderr, "   --underscore -U Set this flag to use an underscore instead of a hyphen when appending --libnames.\n");
}

/**
 * Figures out how to modify the text of a barcode string before printing, based on user
 * options.
 */
void mod_bc_libname(string& bc_str, const string& libname, bool cellranger, bool seurat, bool underscore){
    if (libname != ""){
        if (cellranger){
            bc_str += "-1";
        }
        if (seurat){
            bc_str = libname + "_" + bc_str;
        }
        else{
            if (underscore){
                bc_str += "_" + libname;
            }
            else{
                bc_str += "-" + libname;
            }
        }
    }
}

/**
 * Parse a file output by demux_mt or demux_vcf, and store barcodes mapped
 * to individual IDs.
 */
void parse_barcode_map(string& fn, 
    map<unsigned long, string>& bc2hap,
    set<string>& barcode_groups,
    double llr_cutoff,
    bool keep_doublets){

    ifstream infile(fn.c_str());
    
    string bc_str;
    string hap_str;
    char singdoub;
    double llr;

    while (infile >> bc_str >> hap_str >> singdoub >> llr){
        if ((keep_doublets || singdoub == 'S') && llr >= llr_cutoff){
            // Trim any trailing stuff off of the barcode, so it can
            // be interpreted as a bitset
            unsigned long bc_hashed = bc_ul(bc_str); 
            bc2hap.insert(make_pair(bc_hashed, hap_str));
            barcode_groups.insert(hap_str);
        }
    }
}

/**
 * For doublet identification (in a data structure), convert
 * a combination of haplotype indices i and j into a single index.
 */
short hap_comb_to_idx(short i, short j, short nhaps){
    if (i >= nhaps || j >= nhaps){
        return -1;
    }
    if (i > j){
        // Swap.
        short tmp = j;
        j = i;
        i = tmp;
    }
    // The overall index in the data structure (unique idx for combination)
    short combined_idx = nhaps;
    // The "i" haplotype in the combination
    short first_idx = 0;
    // The "j" haplotype in the combination
    short second_idx = 1;
    short k = 1;
    while (first_idx < i){
        combined_idx += (nhaps - k);
        first_idx++;
        second_idx++;
        k++;
        if (k >= nhaps){
            return -1;
        }
    }
    // first_idx should now == i. Find j
    while (second_idx < j){
        second_idx++;
        combined_idx++;
    }
    return combined_idx;
}

/**
 * Undo the above combination - convert a combined haplotype idx (for doublets)
 * into a <i, j> pair denoting the two indices that went into creating it.
 */
pair<short, short> idx_to_hap_comb(short idx, short nhaps){
    if (idx < nhaps){
        return make_pair(-1, -1);
    }
    short combined_idx = nhaps;
    short i = 0;
    short j = 1;
    short k = 1;
    while (combined_idx < idx){
        if (idx > combined_idx && idx < combined_idx + (nhaps - k)){
            // Increment j till we find it.
            j += (idx - combined_idx);
            return make_pair(i, j);
        }
        combined_idx += (nhaps - k);
        j++;
        i++;
        k++;
        if (k >= nhaps){
            return make_pair(-1, -1);
        }       
    }
    if (combined_idx == idx){
        return make_pair(i, j);
    }
    return make_pair(-1, -1);
}

/**
 * Given a numeric index (single or doublet combination) and a vector
 * of (string) sample names, returns the name of the given sample
 * or sample combination.
 *
 * Ensures that names of doublet combinations are always given in 
 * alphabetic order.
 */
string idx2name(int x, vector<string>& samples){
    string indv_name;
    if (x < samples.size()){
        indv_name = samples[x];
    }
    else{
        pair<short, short> hc = idx_to_hap_comb(x, samples.size());
        if (samples[hc.first] < samples[hc.second]){
            indv_name = samples[hc.first] + "+" + samples[hc.second];
        }
        else{
            indv_name = samples[hc.second] + "+" + samples[hc.first];
        }
    }
    return indv_name;
}

/**
 * Initialize a pre-declared distance matrix by setitng all elements
 * that will be used and accessed to 0 and all others to -1.
 */
void init_distmat(vector<vector<float> >& dist_mat, int dim){
    for (int i = 0; i < dim; ++i){
        vector<float> row;
        for (int j = 0; j <= i; ++j){
            row.push_back(-1);
        }
        for (int j = i +1; j < dim; ++j){
            row.push_back(0);
        }
        dist_mat.push_back(row);
    }
}

/**
 * Print distance matrix to stderr.
 */
void print_distmat(vector<vector<float> >& dist_mat){
    for (int i = 0; i < dist_mat.size(); ++i){
        for (int j = i + 1; j < dist_mat.size(); ++j){
            if (j > i + 1){
                fprintf(stdout, " ");
            }
            fprintf(stdout, "%0.2f", dist_mat[i][j]);
        }
        fprintf(stdout, "\n");
    }
}

void print_distmat_square(vector<vector<float> >& dist_mat){
    for (int i = 0; i < dist_mat.size(); ++i){
        for (int j = 0; j < i; ++j){
            if (j != 0){
                fprintf(stdout, "\t");
            }
            fprintf(stdout, "%0.2f", dist_mat[j][i]);
        }
        fprintf(stdout, "\t0.00\t");
        for (int j = i + 1; j < dist_mat.size(); ++j){
            fprintf(stdout, "\t%0.2f", dist_mat[i][j]);
        }
        fprintf(stdout, "\n");
    }
}

/**
 * Given a table of log likelihood ratios between each possible pair of 
 * identities for a cell, iteratively eliminates the least likely possibility
 * (the individual deemed less likely by the highest-magnitude LLR in the table)
 * until only two individuals are left. The more likely of those individuals
 * is chosen as the true individual, and the LLR between the two final individuals
 * is chosen as the LLR of assignment.
 */
int collapse_llrs(map<int, map<int, double> >& llrs, double& llr_final){
    bool done = false;
    int idx_elim = -1;
    while (!done){
        int min_idx = -1;
        int max_idx = -1;
        double min_llr;
        int ncomps = 0;
        for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end(); ){
            if (idx_elim != -1 && llr->first == idx_elim){
                llrs.erase(llr++);
            }
            else{
                for (map<int, double>::iterator llr2 = llr->second.begin(); 
                    llr2 != llr->second.end(); ){
                    if (idx_elim != -1 && llr2->first == idx_elim){
                        llr->second.erase(llr2++);   
                    } 
                    else{
                        if (min_idx == -1){
                            if (llr2->second < 0){
                                min_idx = llr->first;
                                max_idx = llr2->first;
                                min_llr = llr2->second;
                            }
                            else{
                                min_idx = llr2->first;
                                max_idx = llr->first;
                                min_llr = -llr2->second;
                            }
                        }
                        else if (abs(llr2->second) > abs(min_llr)){
                            if (llr2->second < 0){
                                min_idx = llr->first;
                                max_idx = llr2->first;
                                min_llr = llr2->second;
                            }
                            else{
                                min_idx = llr2->first;
                                max_idx = llr->first;
                                min_llr = -llr2->second;
                            }
                        }
                        ++ncomps;
                        ++llr2;
                    }
                }
                ++llr;
            }
        }

        if (ncomps == 0 || min_idx == -1){
            llr_final = 0;
            return -1;
        }
        else if (ncomps == 1){
            llr_final = -min_llr;
            return max_idx;
        }
        else{
            idx_elim = min_idx;   
        }
    }
    // Nothing found
    llr_final = 0.0;
    return -1;
}

/**
 * Given counts of identifications in a sample, compares the 
 * counts of different doublet combinations to expectations,
 * baed on the frequencies of each individual in the sample.
 *
 * Returns a chi-squared goodness of fit p-value.
 *
 * Low p-values indicate possibly incorrect doublet 
 * identifications, which suggests demultiplexing was
 * inaccurate.
 */
double doublet_chisq(map<int, int>& idcounts, int n_samples){
    
    if (n_samples <= 2){
        // Can't do test with up to only one doublet type
        return -1.0;
    }

    int tot_single = 0;
    int tot_double = 0;
    map<int, int> singles;
    map<int, int> doubles;
    for (map<int, int>::iterator ic = idcounts.begin(); ic != idcounts.end();
        ++ic){
        if (ic->first < n_samples){
            tot_single += ic->second;
            singles.insert(make_pair(ic->first, ic->second));
        }
        else{
            tot_double += ic->second;
            doubles.insert(make_pair(ic->first, ic->second));
        }
    }
     
    // If no doublets, can't do anything
    if (tot_double == 0){
        return 1.0;
    }

    // Get frequency of each single combination
    map<int, double> singfreq;
    for (map<int, int>::iterator s = singles.begin(); s != singles.end(); ++s){
        singfreq.insert(make_pair(s->first, (double)s->second/(double)tot_single));
    }
    
    // Check for missing doublet combinations, and store a count of 0 for each 
    for (int i = 0; i < n_samples-1; ++i){
        for (int j = i + 1; j < n_samples; ++j){
            int k = hap_comb_to_idx(i, j, n_samples);
            if (doubles.count(k) == 0){
                doubles.insert(make_pair(k, 0));
            }
        }
    }

    // Get expectation of each doublet combination
    map<int, double> doubfreq;
    double doubfreq_tot = 0.0; // Need to re-scale probs since not including self+self doublets
    for (map<int, int>::iterator d = doubles.begin(); d != doubles.end(); ++d){
        pair<int, int> combo = idx_to_hap_comb(d->first, n_samples);
        double expected = singfreq[combo.first] * singfreq[combo.second];
        doubfreq_tot += expected;
        doubfreq.insert(make_pair(d->first, expected));
    }
    for (map<int, double>::iterator df = doubfreq.begin(); df != doubfreq.end(); ++df){
        df->second /= doubfreq_tot;
    }
    
    // Now do Chi-squared test
    double chisq = 0.0;
    int df = 0;
    for (map<int, int>::iterator d = doubles.begin(); d != doubles.end(); ++d){
        double expected = (double)tot_double*doubfreq[d->first];
        if (expected > 0){
            chisq += pow((double)d->second - expected, 2) / expected;
        }
        ++df;
    }
    // Degrees of freedom = number of categories minus 1
    df -= 1;
    return pchisq(chisq, (double)df);
}

/**
 * Trim the directory off of a full filename path
 */
string filename_nopath(string& filename){
    size_t trim_idx = filename.find_last_of("\\/");
    if (trim_idx != string::npos){
        return filename.substr(trim_idx + 1, filename.length() - trim_idx - 1);
    }
    else{
        return filename;
    }
}

/**
 * Log PDF of binomial distribution wrt n, k, p
 */
double logbinom(double n, double k, double p){
    double ll = k * log(p) + (n-k)*log(1.0-p);
    // Compute log binomial coefficient
    if (k < n && k != 0){
        // Use Stirling's approximation
        double logn = log(n);
        double logk = log(k);
        double logn_k = log(n-k);
        ll += n*logn - k*logk - (n-k)*logn_k + 0.5*(logn - logk - logn_k - log(2*M_PI));
    }
    return ll;
}

bool file_exists(string filename){
    struct stat buf;
    if (stat(filename.c_str(), &buf) != 0){
        return false;
    }
    return true;
}

/**
 * Given a histogram, approximates the first derivative of
 * the histogram and puts it in datp. Does not apply
 * kernel smoothing.
 */
void derivative(map<double, double>& dat, map<double, double>& datp, int smooth){
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
    
    if (smooth > 0){
        // Apply smoothing
        vector<double> keys;
        vector<double> sums;
        for (map<double, double>::iterator d = datp.begin(); d != datp.end(); ++d){
            keys.push_back(d->first);
            sums.push_back(d->second);
        }
        for (int i = 0; i < keys.size(); ++i){
            int count = 0;
            int jlim = i-smooth+1;
            if (jlim < 0){
                jlim = 0;
            }
            for (int j = i-1; j >= jlim; --j){
                count++;
                sums[i] += datp[keys[j]];
            }
            jlim = i+smooth-1;
            if (jlim > keys.size()-1){
                jlim = keys.size()-1;
            }
            for (int j = i + 1; j <= jlim; ++j){
                count++;
                sums[i] += datp[keys[j]];
            }
            sums[i] /= (double)count;
        }
        datp.clear();
        for (int i = 0; i < keys.size(); ++i){
            datp.insert(make_pair(keys[i], sums[i]));
        } 
    }
}

/**
 * Given a histogram, finds the "knee" point, defined as the
 * point of maximum curvature. Returns the x-coordinate of
 * the maximum curvature point.
 */
double find_knee(map<double, double>& x, double min_frac_to_allow){
    int smooth = 3;

    map<double, double> xprime1;
    derivative(x, xprime1, smooth);
    map<double, double> xprime2;
    derivative(xprime1, xprime2, smooth);

    double maxk = -1;
    double maxk_x = -1;
    
    vector<double> keys;
    vector<double> curv;

    for (map<double, double>::iterator xp = xprime1.begin(); xp != 
        xprime1.end(); ++xp){
        
        // Compute curvature at this point 
        double xk = abs(xprime2[xp->first]) / 
            pow(1 + pow(xp->second, 2), 1.5);
        keys.push_back(xp->first);
        curv.push_back(xk);    
    }

    // How many total cells?
    double tot_dat = x.begin()->second;

    // Do smoothing
    vector<double> smoothed;
    for (int i = 0; i < keys.size(); ++i){
        smoothed.push_back(curv[i]);
        int jlim_low = i-smooth+1;
        if (jlim_low < 0){
            jlim_low = 0;
        }
        int jlim_high = i+smooth-1;
        if (jlim_high > keys.size()-1){
            jlim_high = keys.size()-1;
        }
        int count = 0;
        for (int j = i -1 ; j >= jlim_low; --j){
            smoothed[i] += curv[j];
            count++;
        }
        for (int j = i + 1; j <= jlim_high; ++j){
            smoothed[i] += curv[j];
            count++;
        }
        smoothed[i] /= (double)count;
        if (x[keys[i]] >= min_frac_to_allow*tot_dat){
            if (maxk_x == -1 || smoothed[i] > maxk){
                maxk = smoothed[i];
                maxk_x = keys[i];
            }
        }
    }

    return maxk_x;
}

/**
 * Parses input files in Market exchange format (MEX),
 * used to represent sparse matrices. Takes as input the three MEX files
 * (gzipped or not), along with two data structures to store results.
 * counts will store (hashed) cell barcodes mapped to counts of label indices,
 * and labels will store the names of the labels (i.e. genes).
 */
void parse_mex(const string& barcodesfile,
    const string& featuresfile,
    const string& matrixfile, 
    robin_hood::unordered_map<unsigned long, map<int, long int> >& counts,
    vector<string>& labels,
    const string& featuretype){
     
    // First, parse barcodes
    vector<unsigned long> barcodes;
    parse_barcode_file(barcodesfile, barcodes);
    
    set<string> unique_featuretype; 
    
    // Map feature index in whole file to feature index in subset 
    map<int, int> feature_inds;
    
    // Then, parse genes
    // Feature index in file
    int feature_idx = 0;
    // Feature index in subset
    int feature_idx_subset = 0;

    gzreader labelsreader(featuresfile);
    int nfeatures_match = 0;
    while(labelsreader.next()){
        istringstream splitter(labelsreader.line);
        string field;
        int idx = 0;
        string id = "";
        string name = "";
        string type = "";
        while(getline(splitter, field, '\t')){
            if (idx == 0){
                id = field;
            }
            else if (idx == 1){
                name = field;
            }
            else if (idx == 2){
                type = field;
            }
            idx++;
        }
        if (featuretype == "" || type == "" || featuretype == type){
            if (name != ""){
                labels.push_back(name);
            }
            else{
                labels.push_back(id);
            }
            feature_inds.insert(make_pair(feature_idx, feature_idx_subset));
            ++feature_idx_subset;
            ++nfeatures_match;
        }
        if (featuretype == "" && featuretype != ""){
            unique_featuretype.insert(featuretype);
        }
        ++feature_idx;
    }
    
    if (nfeatures_match == 0){
        fprintf(stderr, "ERROR: 0 features match provided type %s\n", featuretype.c_str());
        fprintf(stderr, "Allowed feature types:\n");
        for (set<string>::iterator ut = unique_featuretype.begin(); ut != unique_featuretype.end();
            ++ut){
            fprintf(stderr, "%s\n", ut->c_str());
        }
        exit(1);
    }
    
    set<int> bciuniq;

    bool barcodes_first = false;
    
    // Then, populate data structure
    gzreader mtxreader(matrixfile);
    int mexline = 0;
    while(mtxreader.next()){
        if (mtxreader.line[0] == '%'){
            continue;
        }
        else{
            if (mexline == 0){
                // Read header.
                istringstream splitter(mtxreader.line);
                long int n1;
                long int n2;
                int idx = 0;
                string token;
                while(getline(splitter, token, ' ' )){
                    if (idx == 0){
                        n1 = atol(token.c_str());
                    }
                    else if (idx == 1){
                        n2 = atol(token.c_str());
                    }
                    ++idx;
                }
                if (n1 == barcodes.size()){
                    barcodes_first = true;
                }
                else if (n2 == barcodes.size()){
                    barcodes_first = false;
                }
                else{
                    fprintf(stderr, "ERROR: %ld barcodes; does not match MTX file: %ld or %ld\n",
                        barcodes.size(), n1, n2);
                    exit(1);
                }
            }
            else if (mexline > 0){
                istringstream splitter(mtxreader.line);
                int idx = 0;
                string token;
                int bc_idx;
                int feature_idx;
                long int count;
                while (getline(splitter, token, ' ')){
                    if ((!barcodes_first && idx == 0) || (barcodes_first && idx == 1)){
                        // feature index
                        feature_idx = atoi(token.c_str());
                        feature_idx--; // Make 0-based
                    }
                    else if ((!barcodes_first && idx == 1) || (barcodes_first && idx == 0)){
                        // barcode index
                        bc_idx = atoi(token.c_str());
                        bc_idx--; // Make 0-based
                    }
                    else if (idx == 2){
                        // UMI count
                        double count_d = atof(token.c_str());
                        count = (long int)round(count_d);
                        
                        if (feature_inds.count(feature_idx) > 0){ 
                            // This feature is in the list
                            bciuniq.insert(bc_idx);
                            unsigned long barcode = barcodes[bc_idx];
                            if (counts.count(barcode) == 0){
                                map<int, long int> m;
                                counts.emplace(barcode, m);
                            }
                            counts[barcode].emplace(feature_inds[feature_idx], count);
                        }
                    }
                    ++idx;
                }
            }
            ++mexline;
        }
    }
    if (featuretype == "" && unique_featuretype.size() > 1){
        fprintf(stderr, "ERROR: no feature type filter provided, but %ld feature types\n", 
            unique_featuretype.size());
        fprintf(stderr, "encountered in MEX data\n");
        exit(1);
    }
    fprintf(stderr, "Loaded %ld barcodes and %d features\n", barcodes.size(), nfeatures_match);

}

/**
 * Log likelihood of dirichlet distribution
 *
 * Format expected by optimML::multivar_ml_solver for finding MLE parameters
 */
double ll_dirichlet(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
   
    char buf[30];
    string bufstr;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    for (int i = 0; i < params.size(); ++i){
        sprintf(&buf[0], "f_%d", i);
        bufstr = buf;
        double f = data_d.at(bufstr);
        term1 += params[i];
        term2 += lgammaf(params[i]);
        term3 += (params[i] - 1.0)*log(f);
    }
    
    return lgammaf(term1) - term2 + term3;
}

/**
 * Derivative of log likelihood of dirichlet distribution wrt each concentration
 * parameter.
 *
 * Format expected by optimML::multivar_ml_solver for finding MLE parameters
 */
void dll_dirichlet(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    char buf[30];
    string bufstr;
    double term1 = 0.0;
    for (int i = 0; i < params.size(); ++i){
        sprintf(&buf[0], "f_%d", i);
        bufstr = buf;
        double f = data_d.at(bufstr);
        term1 += params[i];
        results[i] -= digamma(params[i]) + log(f);
    }
    for (int i = 0; i < params.size(); ++i){
        results[i] += digamma(term1);
    }
}

/**
 * Given the maximum likelihood parameters of a mixture component
 * problem - where we are inferring the proportion of a mixture
 * made up of each of a set of individuals that must each be 0 < x < 1
 * and sum to 1, which has been bootstrapped, fits a Dirichlet
 * distribution to the bootstrap samples and provides the MLE
 * Dirichlet concentration parameters.
 *
 * Parameters:
 *   mle_fracs: maximum likelihood estimates of mixture components
 *   dirprops: Vector of length == mle_fracs.size()
 *      each sub-vector contains every bootstrap sample of that specific
 *      mixture component.
 *   dirichlet_mle: vector to store result concentration parameters.
 *      Will be modified by reference.
 */
void fit_dirichlet(vector<double>& mle_fracs,
    vector<vector<double> >& dirprops,
    vector<double>& dirichlet_mle){
    
    // How many mixture components are there?
    int n_samples = mle_fracs.size();

    // Ensure result vector is empty
    dirichlet_mle.clear();

    // Come up with initial guesses for each concentration parameter.
    vector<double> dir_init;
    // Initialize to max likelihood params with n=1 observation
    for (int j = 0; j < n_samples; ++j){
        dir_init.push_back(mle_fracs[j]);
    }
    
    // Create solver
    optimML::multivar_ml_solver dirsolver(dir_init, ll_dirichlet, dll_dirichlet);
    char buf[30];
    string bufstr;
    for (int j = 0; j < n_samples; ++j){
        // Concentration parameters must be positive
        dirsolver.constrain_pos(j);
        sprintf(&buf[0], "f_%d", j);
        bufstr = buf;
        dirsolver.add_data(bufstr, dirprops[j]);
    }
    
    dirsolver.solve();
    for (int j = 0; j < n_samples; ++j){
        dirichlet_mle.push_back(dirsolver.results[j]);
    }
}
