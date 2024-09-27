/**
 * 
 * Suffix tree to store data mapped to k-mers
 *
 * This data structure was originally conceived of and written by 
 * Ed Green, ed@soe.ucsc.edu
 *
 */
#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include "kmsuftree.h"

/**
 * Initialize kmer suffix tree.
 *
 * Args:
 *  int k the (complete) value of k for all k-mers to be stored
 * Returns:
 *  A pointer to the newly allocated data structure
 */
kmer_tree_p init_kmer_tree( int k ) {
    kmer_tree_p ks;
    size_t i, len;
    len = 1<<(K_AR_SIZE*2); // Length of the array for the array part
    ks = (kmer_tree_p)malloc(sizeof(Kmers));

    ks->k_ar_size = K_AR_SIZE;
    ks->k = k;
    ks->lookup_array = (ktnP*)malloc(sizeof(ktnP) * len);
    for( i = 0; i < len; i++ ) {
        ks->lookup_array[i] = NULL;
    }
    return ks;
}

/**
 * Initialize a k-mer suffix tree node.
 *
 * Sets all pointers to NULL.
 */
ktnP init_ktn( void ) {
    ktnP new_ktn = (ktnP)malloc(sizeof(ktn));
    new_ktn->Ap = NULL;
    new_ktn->Cp = NULL;
    new_ktn->Gp = NULL;
    new_ktn->Tp = NULL;
    new_ktn->dat = NULL;
    return new_ktn;
}

/**
 * Stores data in the kmer tree. If the k-mer already exists in the tree,
 * its data will be overwritten.
 *
 * Args:
 *  const char* kmer the k-mer to use as the key
 *  kmer_tree_p ks The k-mer tree data structure
 *  void* node The item being stored
 *  int rc Should this k-mer be reverse-complemented when keying into the
 *      data structure?
 *  Returns:
 *      0 = failure
 *      1 = success
 */
int kmer_tree_add( const char* kmer, kmer_tree_p ks, void* node, int rc ) {
    size_t inx;
    ktnP curr_node, next_node;
//    size_t kmer_pos;
    char next_base;

    /* Find what the index position is for the first K_AR_SIZE
     * bases in this kmer 
     */

    if ( kmer2inx( kmer, ks->k_ar_size, &inx, rc, ks->k ) ) {

        /* Start at the index position of the first ks->k_ar_size bases */
        curr_node = ks->lookup_array[inx];


        /* If we've never seen that before, then initialize it */
        if ( curr_node == NULL ) {
            curr_node = init_ktn();
            ks->lookup_array[inx] = curr_node;
        }
        
        // Go backwards if reverse-complemented
        if (rc){
            for (int kp = (ks->k - ks->k_ar_size-1); kp >= 0; --kp){  
                next_base = toupper(kmer[kp]);
                switch (next_base){
                    case 'T':
                        if (curr_node->Ap == NULL){
                            curr_node->Ap = init_ktn();
                        }
                        next_node = curr_node->Ap;
                        break;
                    case 'G':
                        if (curr_node->Cp == NULL){
                            curr_node->Cp = init_ktn();
                        }
                        next_node = curr_node->Cp;
                        break;
                    case 'C':
                        if (curr_node->Gp == NULL){
                            curr_node->Gp = init_ktn();
                        }
                        next_node = curr_node->Gp;
                        break;
                    case 'A':
                        if (curr_node->Tp == NULL){
                            curr_node->Tp = init_ktn();
                        }
                        next_node = curr_node->Tp;
                        break;
                    default:
                        return 0;
                        break;
                }
            }
        }
        else{
            for (int kmer_pos = ks->k_ar_size; kmer_pos < ks->k; ++kmer_pos){
                next_base = toupper(kmer[kmer_pos]);
                switch( next_base ) {
                    case 'A' :
                        if ( curr_node->Ap == NULL ) {
                          curr_node->Ap =  init_ktn();
                        }
                        next_node = curr_node->Ap;
                        break;
                    case 'C' :
                        if ( curr_node->Cp == NULL ) {
                            curr_node->Cp = init_ktn();
                        }
                        next_node = curr_node->Cp;
                        break;
                    case 'G' :
                        if ( curr_node->Gp == NULL ) {
                            curr_node->Gp = init_ktn();
                        }
                        next_node = curr_node->Gp;
                        break;
                    case 'T' :
                        if ( curr_node->Tp == NULL ) {
                            curr_node->Tp = init_ktn();
                        }
                        next_node = curr_node->Tp;
                        break;
                    default :
                        // not a good base, not a good kmer, we're done
                        return 0;
                        break;
                }
                curr_node = next_node;
            }
        }
        // Store data in data structure
        curr_node->dat = node;
        return 1;
    }
    else {
        return 0;
    }
}

/**
 * Looks up the data mapped to a specific k-mer. Returns NULL on failure.
 *
 * Args:
 *  const char* kmer the k-mer to look up
 *  const kmer_tree_p ks A pointer to the k-mer tree data structure
 *  int rc - should the k-mer be reverse complemented?
 *
 * Returns:
 *  A pointer to the data mapped to the k-mer. NULL if not found.
 *
 */
void* kmer_tree_lookup( const char* kmer, const kmer_tree_p ks, int rc) {
    size_t inx, kmer_pos;
    ktnP curr_node;
    char curr_base;
    
    /* First, find the inx of the array part, i.e., the beginning */
    if ( kmer2inx( kmer, ks->k_ar_size, &inx, rc, ks->k ) ) {
        curr_node = ks->lookup_array[inx];
    }
    else { 
        return NULL;
    }

    /* If curr_node == NULL, we've never seen this kmer in a 
     cluster before.  */
    if ( curr_node == NULL ) {
        return NULL;
    }
    
    /* Now, follow the pointers until you get to the end of the kmer
     or until we never saw this kmer before */
    
    // Go backwards if reverse complemented
    if (rc){
        for (int kp = (ks->k - ks->k_ar_size-1); kp >= 0; --kp){
            curr_base = toupper(kmer[kp]);
            switch (curr_base){
                case 'T':
                    curr_node = curr_node->Ap;
                    break;
                case 'G':
                    curr_node = curr_node->Cp;
                    break;
                case 'C':
                    curr_node = curr_node->Gp;
                    break;
                case 'A':
                    curr_node = curr_node->Tp;
                    break;
                default:
                    curr_node = NULL;
                    break;
            }
            if (curr_node == NULL){
                return NULL;
            }
        }
    }
    else{
        for( kmer_pos = ks->k_ar_size; kmer_pos < ks->k; kmer_pos++ ) {
            curr_base = toupper(kmer[kmer_pos]);
            switch( curr_base ) {
                case 'A' :
                    curr_node = curr_node->Ap;
                    break;
                case 'C' :
                    curr_node = curr_node->Cp;
                    break;
                case 'G' :
                    curr_node = curr_node->Gp;
                    break;
                case 'T' :
                    curr_node = curr_node->Tp;
                    break;
                default :
                    curr_node = NULL;
                    break;
            }
            if ( curr_node == NULL ) {
                return NULL;
            }
        }
    }
    
    return curr_node->dat;
}

int kmer_tree_remove_aux(char* kmerchr, 
    ktnP curr_node, 
    int bases_left, 
    int rc, 
    int careful){
    
    // Base case: we are at the destination.
    if (bases_left == 0){
        if (!careful){
            free(curr_node->dat);
        }
        free(curr_node);
        return 1;
    }
    
    // Otherwise, go through the tree to the destination. If the destination
    // is found, check to see if the next level of nodes is free-able.
    char curr_base = toupper(kmerchr[0]);
    
    ktnP nextP = NULL;
    switch(curr_base){
        case 'A':
            if (!rc){
                if (kmer_tree_remove_aux(kmerchr+1, curr_node->Ap, bases_left-1, rc, careful)){
                    curr_node->Ap = NULL;
                }
            }
            else{
                if (kmer_tree_remove_aux(kmerchr-1, curr_node->Tp, bases_left-1, rc, careful)){
                    curr_node->Tp = NULL;
                }
            }
        break;
        case 'C':
            if (!rc){
                if (kmer_tree_remove_aux(kmerchr+1, curr_node->Cp, bases_left-1, rc, careful)){
                    curr_node->Cp = NULL;
                }
            }
            else{
                if (kmer_tree_remove_aux(kmerchr-1, curr_node->Gp, bases_left-1, rc, careful)){
                    curr_node->Gp = NULL;
                }
            }
        break;
        case 'G':
            if (!rc){
                if (kmer_tree_remove_aux(kmerchr+1, curr_node->Gp, bases_left-1, rc, careful)){
                    curr_node->Gp = NULL;
                }
            }
            else{
                if (kmer_tree_remove_aux(kmerchr-1, curr_node->Cp, bases_left-1, rc, careful)){
                    curr_node->Cp = NULL;
                }
            }
        break;
        case 'T':
            if (!rc){
                if (kmer_tree_remove_aux(kmerchr+1, curr_node->Tp, bases_left-1, rc, careful)){
                    curr_node->Tp = NULL;
                }
            }
            else{
                if (kmer_tree_remove_aux(kmerchr-1, curr_node->Ap, bases_left-1, rc, careful)){
                    curr_node->Ap = NULL;
                }
            }
        break;
        default:
            // Can't be in this tree because it contains an invalid character.
            fprintf(stderr, "C\n");
            return 0;
        break;
    }
    
    if (curr_node->Ap == NULL && curr_node->Cp == NULL && curr_node->Gp == NULL &&
        curr_node->Tp == NULL && curr_node->dat == NULL){
        // This node can be removed too.
        //free(curr_node);
        return 1;
        //return 0;
    }
    else{
        // This node still is connected to something; keep it around.
        return 0;
    }
}

/**
 * Remove a stored item from the data structure. 
 */
int kmer_tree_remove(char* kmer, const kmer_tree_p ks, int rc, int careful){
    
    // This is the same as looking an item up, but when the item is found,
    // we need to recursively remove pointers to empty stuff.
    
    size_t inx, kmer_pos;
    ktnP curr_node;
    char curr_base;
    
    /* First, find the inx of the array part, i.e., the beginning */
    if ( kmer2inx( kmer, ks->k_ar_size, &inx, rc, ks->k ) ) {
        curr_node = ks->lookup_array[inx];
    }
    else { 
        return 0;
    }

    /* If curr_node == NULL, we've never seen this kmer in a 
     cluster before.  */
    if ( curr_node == NULL ) {
        return 0;
    }
    
    char* kmerptr = kmer + ks->k_ar_size;
    int bases_left = (ks->k - ks->k_ar_size);
    if (rc){
        kmerptr = kmer + (ks->k - ks->k_ar_size-1);
    }
    
    return kmer_tree_remove_aux(kmerptr, curr_node, bases_left, rc, careful);
}

/**
 * BE CAREFUL; this function assumes that the data (void*) stored in the
 * kmer suffix tree is integers. If not, this will have weird, undefined
 * behavior.
 *
 * Given a kmer, kmer suffix tree, and integer, looks to see if that k-mer
 *  already has an entry in the tree. If not, it stores the value 1 (as a void*)
 *  in the tree. If so, it increments the stored value by 1.
 *
 * This function will automatically check both the forward and reverse 
 *  complement versions of the k-mer.
 *
 * Returns: 
 *  0 = failure
 *  1 = success
 */
int kmer_tree_increment(const char* kmer, const kmer_tree_p ks){
    
    void* stored = NULL;
    
    // Try to look up in forward orientation.
    stored = kmer_tree_lookup(kmer, ks, 0);
    
    if (stored != NULL){
        // Already exists in forward orientation.
        long int* intptr = (long int*) stored;
        (*intptr)++;
    }
    else{
        // Try to look up in reverse complement orientation.
        stored = kmer_tree_lookup(kmer, ks, 1);
        if (stored != NULL){
            long int* intptr = (long int*) stored;
            (*intptr)++;
        }
        else{
            // Need to create (do so in forward orientation).
            long int* count = (long int*) malloc(sizeof(long int));
            *count = (long int)1;
            kmer_tree_add(kmer, ks, (void*) count, 0);
        }
    }
    return 1;
}

/**
 * Auxiliary recursive helper function for print_kmer_counts().
 * Prints data if we've hit an end point in the tree where data are stored,
 * or keeps iterating through the tree if not.
 *
 * Arguments:
 *  kmer_tree_node* ptr The pointer to the current node in the suffix tree
 *  char* prefix The string buffer to store the k-mer (should be populated
 *      up until the current base)
 *  char curBase the current base in the k-mer (based on the pointer that led
 *      us here)
 *  int kmer_index The base index in the current k-mer
 *  int k The total length of k-mers in the suffix tree
 */
void print_kmer_counts_aux(struct kmer_tree_node* ptr, char* prefix, char curBase, int kmer_index, int k){
    char* baseptr = prefix + kmer_index;
    *baseptr = curBase;
    if (kmer_index == k - 1 && ptr->dat != NULL){
        // Done.
        // Get count and print.
        long int* intptr = (long int*) ptr->dat;
        fprintf(stdout, "%s\t%ld\n", prefix, *intptr);
        return;
    }
    else{
        if (ptr->Ap != NULL){
            print_kmer_counts_aux(ptr->Ap, prefix, 'A', kmer_index+1, k);
        }
        if (ptr->Cp != NULL){
            print_kmer_counts_aux(ptr->Cp, prefix, 'C', kmer_index+1, k);
        }
        if (ptr->Gp != NULL){
            print_kmer_counts_aux(ptr->Gp, prefix, 'G', kmer_index+1, k);
        }
        if (ptr->Tp != NULL){
            print_kmer_counts_aux(ptr->Tp, prefix, 'T', kmer_index+1, k);
        }
    }
}

/**
 * BE CAREFUL; this function assumes that the data (void*) stored in the
 * kmer suffix tree is integers. If not, this will have weird, undefined
 * behavior.
 *
 * Given kmer suffix tree, traverses it and prints out all k-mers and values
 * stored within as if they're integers. This should only be used if the
 * suffix tree is being used as a counter.
 *
 * Arguments:
 *  const kmer_tree_p ks The pointer to the data structure
 */
void print_kmer_counts(const kmer_tree_p ks){
    // First, traverse the lookup array and note any position that
    // doesn't contain a null pointer
    for (size_t index = 0; index < 1<<(K_AR_SIZE*2); ++index){
        if (ks->lookup_array[index] != NULL){
            
            // Populate what we can of the k-mer
            char kmer[ks->k + 1];
            kmer[ks->k] = '\0';
            
            // Convert index to kmer prefix
            inx2kmer(index, ks->k_ar_size, kmer);
            
            //for (int i = ks->k_ar_size; i < ks->k; ++i){
            //    kmer[i] = 'N';
            //}
            
            if (ks->lookup_array[index]->Ap != NULL){
                print_kmer_counts_aux(ks->lookup_array[index]->Ap, 
                    kmer, 'A', ks->k_ar_size, ks->k);
            }
            if (ks->lookup_array[index]->Cp != NULL){
                print_kmer_counts_aux(ks->lookup_array[index]->Cp, 
                    kmer, 'C', ks->k_ar_size, ks->k);
            }
            if (ks->lookup_array[index]->Gp != NULL){
                print_kmer_counts_aux(ks->lookup_array[index]->Gp, 
                    kmer, 'G', ks->k_ar_size, ks->k);
            }
            if (ks->lookup_array[index]->Tp != NULL){
                print_kmer_counts_aux(ks->lookup_array[index]->Tp, 
                    kmer, 'T', ks->k_ar_size, ks->k);
            }
        }
    }
}

/**
 * Translates a string k-mer into a size_t (unsigned long) representation,
 * which can be used as an index in the lookup_array for the first stretch
 * of bases in a k-mer.
 *
 * Arguments:
 *  const char* kmer The string containing the k-mer. Might not be 
 *      null-terminated
 *  const size_t kmer_len The length of the k-mer
 *  size_t* inx The pointer to the size_t that will store the index
 *
 * Returns:
 *  1 / true if the index was set
 *  0 / false if the index was not set
 *
 * Notes:
 *  Uses the formula A => 00, C => 01, G => 10, T => 11 to make a bitstring
 *  for the k-mer. Any other character (i.e. N) is not allowed and will
 *  cause the function to return 0. 
 *  The bitstring is constructed by reading the k-mer from left to right.
 *  This bitstring is then interpreted as a variable of type size_t and
 *  is appropriate as an array index.
 */
int kmer2inx( const char* kmer,
    const size_t kmer_len,
    size_t* inx,
    int rc,
    int k ) {
    
    size_t l_inx  = 0;
    char curr_char;

    if (rc){
        for (int i = k-1; i >= (k-kmer_len); --i){
            l_inx = l_inx << 2;
            curr_char = toupper(kmer[i]);
            switch( curr_char ){
                case 'T':
                    l_inx += 0;
                    break;
                case 'G':
                    l_inx += 1;
                    break;
                case 'C':
                    l_inx += 2;
                    break;
                case 'A':
                    l_inx += 3;
                    break;
                default:
                    return 0;
                    break;
            }
        }
    }
    else{
        for (int i = 0; i < kmer_len; ++i){
            l_inx = l_inx << 2;
            curr_char = toupper(kmer[i]); // Upper case it in case it is not
            switch( curr_char ) {
                case 'A' :
                    l_inx += 0;
                    break;
                case 'C' :
                    l_inx += 1;
                    break;
                case 'G' :
                    l_inx += 2;
                    break;
                case 'T' :
                    l_inx += 3;
                    break;
                default :
                    return 0; // not valid!
                    break;
            }
            
        }
    }
    
    *inx = l_inx;
    return 1; // valid!
}

/**
 * Translates an unsigned long index value into the lookup array for
 * a k-mer suffix tree into the k-mer it represents.
 *
 * Arguments:
 *  const size_t inx The index into the lookup array
 *  const int k The total length of k-mers in the suffix tree
 *  char* kmer The string buffer to populate with the k-mer
 */
void inx2kmer( const size_t inx, const int k, char* kmer ) {
    size_t i;
    size_t code = 0;
    size_t mask = 3; // just the two lowest bits
    for( i = 1; i <= k; i++ ) {
        code = inx >> ((k*2) - (i*2));
        code = (code & mask);

        switch(code) {
        case 0 :
            kmer[i-1] = 'A';
            break;
        case 1 :
            kmer[i-1] = 'C';
            break;
        case 2 :
            kmer[i-1] = 'G';
            break;
        case 3 :
            kmer[i-1] = 'T';
            break;
        }
    }
    kmer[k] = '\0';
}

/**
 * Destructor for a node in the kmsufftree.
 *
 * Arguments:
 *  kmer_tree_node* ptr     The node
 *  careful:    1 = don't free stored data
 *              0 = free stored data
*/
void knode_destruct(struct kmer_tree_node* ptr, int careful) {
    // Recursively destruct the other nodes
    if (ptr->Ap != NULL) {
        knode_destruct(ptr->Ap, careful);
    }
    if (ptr->Cp != NULL) {
        knode_destruct(ptr->Cp, careful);
    }
    if (ptr->Gp != NULL) {
        knode_destruct(ptr->Gp, careful);
    }
    if (ptr->Tp != NULL) {
        knode_destruct(ptr->Tp, careful);
    }
    if (!careful){
        free(ptr->dat);
    }
    free(ptr);  
}

/**
 * Destructor for the kmsuftree.
 *
 * Arguments:
 *  kmer_tree_p ks  The suffix tree
*/
void kmsuftree_destruct(kmer_tree_p ks, int careful) {

    for (size_t index = 0; index < 1<<(K_AR_SIZE*2); ++index){
        ktnP curr_node = ks->lookup_array[index];
        if (curr_node!= NULL){
            // Recursively destroy nodes
            knode_destruct(curr_node, careful);
        }
    }

    free(ks->lookup_array);
    free(ks);
}

