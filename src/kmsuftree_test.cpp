#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include <map>
#include <string>
#include <algorithm>
#include "kmsuftree.h"

using namespace std;

string revcomp(const string& seq){
    char seqrc[seq.length() + 1];
    seqrc[seq.length()] = '\0';
    int ind_rc = 0;
    for (int i = seq.length()-1; i >= 0; --i){
        switch(seq[i]){
            case 'A':
                seqrc[ind_rc] = 'T';
            break;
            case 'C':
                seqrc[ind_rc] = 'G';
            break;
            case 'G':
                seqrc[ind_rc] = 'C';
            break;
            case 'T':
                seqrc[ind_rc] = 'A';
            break;
            default:
                seqrc[ind_rc] = 'N';
            break;
        }
        ++ind_rc;
    }
    return string(seqrc);
}

int main(int argc, char *argv[]) {   
    
    // Decide value of k to use
    int k = 17; 
    
    // Initialize data structure (returns pointer to it)
    kmer_tree_p kt = init_kmer_tree(k);
    
    // Also store counts in a map for comparison.
    map<string, int> kmer_map;
    
    fprintf(stderr, "===== Testing insertion and lookup: =====\n\n");
    string seq = "ACGTGGGGTAGAGTAGGACCGCGATTTGAGATGCGACCTACTCTACCCCACGTGGTCGCATCTCAAATCGCGG";
    // Create null-terminated buffer to store k-mers as we read from the sequence
    char kmer[k + 1];
    kmer[k] = '\0';
    
    fprintf(stderr, "counts of forward and reverse k-mers from forward seq:\n");
    for (int i = 0; i < seq.length() - k + 1; ++i){
        strncpy(kmer, &seq[i], k);
        
        // Insert into standard map
        string kmer_f = string(kmer);
        string kmer_r = revcomp(kmer_f);
        if (kmer_map.count(kmer_f) > 0){
            kmer_map[kmer_f]++;
        }
        else if (kmer_map.count(kmer_r) > 0){
            kmer_map[kmer_r]++;
        }
        else{
            kmer_map.insert(make_pair(kmer_f, 1));
        }
        
        // Add a count for the current k-mer to the data structure
        kmer_tree_increment(kmer, kt);
        // Check whether the result was stored in forward orientation
        void* result = kmer_tree_lookup(kmer, kt, 0);
        if (result != NULL){
            fprintf(stderr, "%s F %ld\n", kmer, *(long int*)result);
        }
        // Check whether the result was stored in reverse-complement orientation
        result = kmer_tree_lookup(kmer, kt, 1);
        if (result != NULL){
            fprintf(stderr, "%s R %ld\n", kmer, *(long int*)result);
        }
    }
    // Do the same thing for the whole sequence reverse-complemented
    fprintf(stderr, "\ncounts of forward and reverse k-mers from reverse-complemented seq:\n");
    seq = "CCGCGATTTGAGATGCGACCACGTGGGGTAGAGTAGGTCGCATCTCAAATCGCGGTCCTACTCTACCCCACGT";
    for (int i = 0; i < seq.length() - k + 1; ++i){
        strncpy(kmer, &seq[i], k);
        
        // Insert into standard map
        string kmer_f = string(kmer);
        string kmer_r = revcomp(kmer_f); 
        if (kmer_map.count(kmer_f) > 0){
            kmer_map[kmer_f]++;
        }
        else if (kmer_map.count(kmer_r) > 0){
            kmer_map[kmer_r]++;
        }
        else{
            kmer_map.insert(make_pair(kmer_f, 1));
        }
        
        kmer_tree_increment(kmer, kt);
        void* result = kmer_tree_lookup(kmer, kt, 0);
        if (result != NULL){
            fprintf(stderr, "%s F %ld\n", kmer, *(long int*)result);
        }
        result = kmer_tree_lookup(kmer, kt, 1);
        if (result != NULL){
            fprintf(stderr, "%s R %ld\n", kmer, *(long int*)result);
        }
    }
    
    fprintf(stderr, "\n");
    
    // Traverse the whole data structure and print counts (to stdout)
    fprintf(stderr, "===== Testing print_kmer_counts(): =====\n");
    fprintf(stderr, "\noverall (all) k-mer counts:\n");
    print_kmer_counts(kt);
    
    fprintf(stderr, "\n===== Testing accuracy of stored values: =====\n");
    fprintf(stderr, "\ncomparison of standard map and kmsuftree counts:\n");
    //fprintf(stderr, "\nstandard map counts:\n");
    for (map<string, int>::iterator kc = kmer_map.begin(); kc != kmer_map.end(); ++kc){
        void* ktree_entry = kmer_tree_lookup(kc->first.c_str(), kt, 0);
        long int ktree_count = -1;
        if (ktree_entry == NULL){
            fprintf(stderr, "ERROR: null value encountered\n");
        }
        else{
            ktree_count = *(long int*)ktree_entry;
        }
        fprintf(stderr, "%s\t%d\t%ld\n", kc->first.c_str(), kc->second, ktree_count);
        if (ktree_count != kc->second){
            fprintf(stderr, "ERROR: counts do not match!\n");
        }
    }
    
    // Test removing individual k-mers from the data structure
    fprintf(stderr, "\n===== Testing k-mer removal: =====\n\n");
    int to_remove = 42;
    int removed = 0;
    fprintf(stderr, "removing %d k-mers from lookup table (leaving %ld)\n", to_remove,
        kmer_map.size()-to_remove);
    for (map<string, int>::iterator kc = kmer_map.begin(); kc != kmer_map.end(); ){
        strncpy(kmer, kc->first.c_str(), k);
        if (kmer_tree_lookup(kmer, kt, 0) != NULL){
            fprintf(stderr, "%s\n", kmer);
            bool success = kmer_tree_remove(kmer, kt, 0, 0);
            if (!success){
                fprintf(stderr, "ERROR removing k-mer from tree\n");
            }
            else{
                fprintf(stderr, "success\n");
                void* retrieved = kmer_tree_lookup(kmer, kt, 0);
                if (retrieved == NULL){
                    fprintf(stderr, "\tsuccess: retrieved == NULL after removal\n");
                }
            }
            removed++;
            if (removed >= to_remove){
                break;
            }
        }
        kmer_map.erase(kc++);
    }
    fprintf(stderr, "\nprinting %d remaining k-mers\n", (int)kmer_map.size());
    print_kmer_counts(kt);
    
    fprintf(stderr, "\nremoving remaining k-mers in reverse-complement orientation\n");
    for (map<string, int>::iterator kc = kmer_map.begin(); kc != kmer_map.end(); ){
        string rc = revcomp(kc->first);
        strncpy(kmer, rc.c_str(), k);
        if (kmer_tree_lookup(kmer, kt, 1) != NULL){
            fprintf(stderr, "%s\n", kc->first.c_str());
            bool success = kmer_tree_remove(kmer, kt, 1, 0);
            if (!success){
                fprintf(stderr, "ERROR removing k-mer from tree\n");
            }
            else{
                fprintf(stderr, "success\n");
                void* retrieved = kmer_tree_lookup(kmer, kt, 1);
                if (retrieved == NULL){
                    fprintf(stderr, "\tsuccess: retrieved == NULL after removal\n");
                }
            }
        }
        kmer_map.erase(kc++);
    }
    
    fprintf(stderr, "\nprinting remaining k-mers (should be none):\n");
    print_kmer_counts(kt);
    
    fprintf(stderr, "\n===== Testing destructor method =====\n");
    fprintf(stderr, "\nfreeing remaining k-mer lookup table entries\n");
    kmsuftree_destruct(kt, 0);
    
    fprintf(stderr, "success\n");
    
    
}
