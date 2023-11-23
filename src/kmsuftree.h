#ifndef KMER_SUF_TREE_H
#define KMER_SUF_TREE_H
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
#define MAX_K (35) // biggest K we can deal with
#define K_AR_SIZE (13) // default length of the array size of the kmer structure

typedef struct kmer_tree_node {
    struct kmer_tree_node* Ap;
    struct kmer_tree_node* Cp;
    struct kmer_tree_node* Gp;
    struct kmer_tree_node* Tp;
    void* dat;
} ktn;

typedef struct kmer_tree_node* ktnP;

typedef struct kmers {
    // How many bases per k-mer?
    size_t k;
    // The first x bases per k-mer that will be handled by the array part
    // and not the tree (if 9, then the first 9 bases of each k-mer will be
    // interpreted as a bit string, which in turn will be interpreted as
    // an unsigned long, which is an index into this array)
    size_t k_ar_size ;
    // A pointer to said array
    ktnP* lookup_array;
} Kmers;

typedef struct kmers* kmer_tree_p;

/* Function prototypes */
kmer_tree_p init_kmer_tree( int );
int kmer_tree_add( const char*, kmer_tree_p, void*, int);
void* kmer_tree_lookup( const char*, const kmer_tree_p, int);
int kmer_tree_remove_aux(char*, ktnP, int, int, int);
int kmer_tree_remove(char*, const kmer_tree_p, int, int);
int kmer_tree_increment( const char*, const kmer_tree_p);
void print_kmer_counts_aux(struct kmer_tree_node*, char*, char, int, int);
void print_kmer_counts(kmer_tree_p);
int kmer2inx( const char*, const size_t, size_t*, int, int);
void inx2kmer( const size_t, const int, char* );
ktnP init_ktn( void );
void knode_destruct(struct kmer_tree_node*, int);
void kmsuftree_destruct(kmer_tree_p, int);


#endif // KMER_SUF_TREE_H
