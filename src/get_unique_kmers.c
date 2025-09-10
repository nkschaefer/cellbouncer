#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "FASTK/libfastk.h"

void help(int code){
    fprintf(stderr, "USAGE: get_unique_kmers\n");
    fprintf(stderr, "First, use a utility like gffread (https://github.com/gpertea/gffread)\n");
    fprintf(stderr, "   to extract FASTA sequences of all transcripts in each transcriptome\n");
    fprintf(stderr, "   of interest, if you do not already have transcript FASTAs. To do this,\n");
    fprintf(stderr, "   run gffread -F -w {output_tx.fa} -g {genome_fasta.fa} {unzipped_gtf.gtf}\n");
    fprintf(stderr, "   where items in braces should be changed to your files.\n");
    fprintf(stderr, "   Now, {output_tx.fa} contains FASTA sequences of transcripts.\n");  
    fprintf(stderr, "Next, run the FastK program in the FASTK package (https://github.com/thegenemyers/FASTK)\n");
    fprintf(stderr, "   to count k-mers in each transcriptome. To do this, run:\n");
    fprintf(stderr, "   FastK -N{output_prefix} -k{kmer_size} -t1 {output_tx.fa}\n");
    fprintf(stderr, "   Note that there should be no space between argument names and argument values.\n");
    fprintf(stderr, "   {kmer_size} should be a value long enough to be unique, but short enough to allow\n");
    fprintf(stderr, "   for some errors in reads and to fit into memory. A value between 25-45 should be\n");
    fprintf(stderr, "   reasonable.\n");
    fprintf(stderr, "   After running this program, there should be a file called {output_prefix}.ktab.\n");
    fprintf(stderr, "After running FASTK to count k-mers on all of your species' transcriptomes,\n");
    fprintf(stderr, "   run this prgram to read FASTK's k-mer table files and output lists of\n");
    fprintf(stderr, "   k-mers unique to each species. These can then be used to demultiplex\n");
    fprintf(stderr, "   reads from multiple species without mapping to a reference genome.\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "  -k Two or more FASTK tables from reference genomes (specify -k\n");
    fprintf(stderr, "     multiple times). These are the files ending with .ktab.\n");
    fprintf(stderr, "  -n The names of the species for which you have provided k-mer tables,\n");
    fprintf(stderr, "     in the same order (specify -n multiple names). Names should not\n");
    fprintf(stderr, "     contain spaces or special characters.\n");
    fprintf(stderr, "  -N Sample this many unique k-mers per species (in decreasing order\n");
    fprintf(stderr, "     of frequency). -1 = include all unique k-mers (default).\n");
    fprintf(stderr, "  -o Where to write species-specific kmers. Files will be\n");
    fprintf(stderr, "     in the format {outprefix}.{index}.kmers (gzip compressed).\n"); 
    fprintf(stderr, "  -d Maximum allowable DUST score (measure of sequence complexity;\n");
    fprintf(stderr, "     higher means more repetitive. Default = 2\n");
    exit(code);
}

/**
 * Sort function for FASTK Kmers - needed to find kmers in one
 * but not another set
 */
int kcomp(const Kmer_Stream* k1, const Kmer_Stream* k2){
    if (k1->cpre < k2->cpre){
        return -1;
    }
    else if (k1->cpre > k2->cpre){
        return 1;
    }
    else{
        for (int idx = 0; idx < k1->kbyte - k1->ibyte; ++idx){
            if (k1->csuf[idx] < k2->csuf[idx]){
                return -1;
            }
            else if (k1->csuf[idx] > k2->csuf[idx]){
                return 1;
            }
        }
    }
    return 0;
}

/**
 * Suffix tree node to look up 3-mers for DUST
 */
struct suffixnode{
    struct suffixnode* A;
    struct suffixnode* C;
    struct suffixnode* G;
    struct suffixnode* T;
    int idx;
};

/**
 * Retrieves a 3-mer from the suffix tree and returns its
 * index
 */
int lookup(struct suffixnode* tree, const char* str){
    struct suffixnode* node;
    switch(str[0]){
        case 'A':
        case 'a':
            node = &tree[0];
            break;
        case 'C':
        case 'c':
            node = &tree[1];
            break;
        case 'G':
        case 'g':
            node = &tree[2];
            break;
        case 'T':
        case 't':
            node = &tree[3];
            break;
    }
    for (int i = 1; i < 3; ++i){
        switch(str[i]){
            case 'A':
            case 'a':
                node = node->A;
                break;
            case 'C':
            case 'c':
                node = node->C;
                break;
            case 'G':
            case 'g':
                node = node->G;
                break;
            case 'T':
            case 't':
                node = node->T;
                break;
        }
    }
    return node->idx;
}

/**
 * Destroy suffix tree (helper function)
 */
void free_suftree_aux(struct suffixnode* tree){
    if (tree->idx == -1){
        free_suftree_aux(tree->A);
        free_suftree_aux(tree->C);
        free_suftree_aux(tree->G);
        free_suftree_aux(tree->T);
        tree->A = NULL;
        tree->C = NULL;
        tree->G = NULL;
        tree->T = NULL;
    }
    free(tree);
}

/**
 * Destroy suffix tree
 */
void free_suftree(struct suffixnode* tree){
    for (int i = 0; i < 4; ++i){
        free_suftree_aux(tree[i].A);
        free_suftree_aux(tree[i].C);
        free_suftree_aux(tree[i].G);
        free_suftree_aux(tree[i].T);
    }
}

/**
 * Compute a measure of sequence complexity that we can use to filter
 * k-mers
 */
double dust(char* kmer, int k, struct suffixnode* tree){
    
    // Use DUST algorithm to compute trimer-based measure of 
    // sequence complexity
    
    // References:
    // https://www.researchgate.net/publication/6987052_A_Fast_and_Symmetric_DUST_Implementation_to_Mask_Low-Complexity_DNA_Sequences
    // https://academic.oup.com/bioinformatics/article/27/6/863/236283

    int counts[64];
    for (int i = 0; i < 64; ++i){
        counts[i] = 0;
    }

    for (int start_idx = 0; start_idx < k-3; ++start_idx){
        int idx = lookup(tree, &kmer[start_idx]);
        counts[idx]++;        
    }
    
    double num = 0;
    for (int i = 0; i < 64; ++i){
        num += (double)(counts[i]*(counts[i]-1))/2.0;
    }    
    int l = k-2;
    double S = num/(double)(l-1);
    return S;
}

int print_dat_aux(char* kmer, int k, struct suffixnode* tree, double maxdust){
    double d = dust(kmer, k, tree);

    if (d > maxdust){
        return 0;
    }    

    for (int i = 0; i < k; ++i){
        if (kmer[i] == 'a'){
            kmer[i] = 'A';
        }
        else if (kmer[i] == 'c'){
            kmer[i] = 'C';
        }
        else if (kmer[i] == 'g'){
            kmer[i] = 'G';
        }
        else if (kmer[i] == 't'){
            kmer[i] = 'T';
        }
        else if (kmer[i] != 'A' && kmer[i] != 'C' && 
            kmer[i] != 'G' && kmer[i] != 'T'){
            kmer[i] = 'N';
        }
    }
    return 1;
}

/**
 * Given a k-mer that is unique to one set, checks to see whether it
 * is sufficiently complex, and if it is, writes it to the appropriate
 * output file.
 *
 * Returns 1 if it passes the DUST filter, 0 otherwise
 */
int print_dat(FILE* outf, int k, char* kmer, int count, struct suffixnode* tree, double maxdust){
    if (print_dat_aux(kmer, k, tree, maxdust) == 0){
        return 0;
    }    
    fprintf(outf, "%d\t%s\n", count, &kmer[0]);
    return 1;
}

/**
 * Same as above, but writes to gzFile instead of FILE* 
 */
int print_dat_gz(gzFile* outf, int k, char* kmer, int count, struct suffixnode* tree, double maxdust){
    if (print_dat_aux(kmer, k, tree, maxdust) == 0){
        return 0;
    }
    gzwrite(*outf, &kmer[0], k);
    gzwrite(*outf, "\n", 1);
    return 1;
}

int intcomp(const void* elem1, const void* elem2){
    int x = *((int*)elem1);
    int y = *((int*)elem2);
    if (x > y){
        return 1;
    }
    if (x < y){
        return -1;
    }
    return 0;
}

struct counts{
    int* vals;
    int nvals;
    int size;
};

void counts_init(struct counts* c, int s){
    c->size = s;
    c->nvals = 0;
    c->vals = (int*)malloc(c->size*sizeof(int));
}

void counts_destroy(struct counts* c){
    free(c->vals);
    c->vals = NULL;
    c->size = 0;
    c->nvals = 0;
}

void counts_add(struct counts* c, int item){
    // Make negative so that sorting will be largest-smallest
    c->vals[c->nvals] = -item;
    c->nvals++;
    if (c->nvals == c->size){
        c->size *= 2;
        c->vals = realloc((void*)c->vals, c->size*sizeof(int));
    }
}

int counts_get(struct counts* c, int idx){
    if (idx >= c->nvals){
        fprintf(stderr, "ERROR: idx %d out of range (%d)\n", idx, c->nvals);
        return -1;
    }
    return c->vals[idx];
}

void counts_sort(struct counts* c){
    // Helpful resource: https://stackoverflow.com/a/1788048
    qsort(c->vals, c->nvals, sizeof(int), intcomp);
}

int main(int argc, char* argv[]){
    srand(time(NULL));

    struct suffixnode* tree = (struct suffixnode*)malloc(4*sizeof(struct suffixnode));
    for (int i = 0; i < 4; ++i){
        tree[i].idx = -1;
        tree[i].A = NULL;
        tree[i].C = NULL;
        tree[i].G = NULL;
        tree[i].T = NULL;
    }
    
    int word_idx = 0;
    for (int i = 0; i < 4; ++i){        
        for (int j = 0; j < 4; ++j){
            struct suffixnode* node = (struct suffixnode*)malloc(sizeof(struct suffixnode));
            node->idx = -1;
            switch(j){
                case 0:
                    tree[i].A = node;
                    break;
                case 1:
                    tree[i].C = node;
                    break;
                case 2:
                    tree[i].G = node;
                    break;
                case 3:
                    tree[i].T = node;
                    break;
            }
            for (int k = 0; k < 4; ++k){
                struct suffixnode* node2 = (struct suffixnode*)malloc(sizeof(struct suffixnode));
                node2->idx = word_idx;
                switch(k){
                    case 0:
                        node->A = node2;
                        break;
                    case 1:
                        node->C = node2;
                        break;
                    case 2:
                        node->G = node2;
                        break;
                    case 3:
                        node->T = node2;
                        break;
                }
                ++word_idx;
            }
        }
    }

    // Set defaults
    // Assume all strings will be under 500 char
    char kmers[50][500];
    int kmers_idx = 0;
    char names[50][500];
    int names_idx = 0;
    char readsfile[500];
    char proffile[500];
    char output_prefix[500];   
    readsfile[0] = '\0';
    output_prefix[0] = '\0';
    double maxdust = 2;
    //int num_samp = 10000000;
    // Default to retrieving all usable k-mers
    int num_samp = -1;

    int option_index = 0;
    int ch;
    if (argc == 1){
        help(0);
    }
    
    // Parse arguments
    while ((ch = getopt(argc, argv, "k:n:o:d:N:h")) != -1){
        switch(ch){
            case 'k':
                strcpy(&kmers[kmers_idx][0], optarg);
                kmers[kmers_idx][strlen(optarg)] = '\0';
                kmers_idx++;
                break;
            case 'n':
                strcpy(&names[names_idx][0], optarg);
                names[names_idx][strlen(optarg)] = '\0';
                names_idx++;
                break;
            case 'd':
                maxdust = atof(optarg);
                break;
            case 'N':
                num_samp = atoi(optarg);
                break;
            case 'h':
                help(0);
                break;
            case 0:
                break;
            case 'o':
                strcpy(&output_prefix[0], optarg);
                break;
            default:
                help(0);
                break;
        }
    }
    // Validate args
    if (strlen(output_prefix) == 0){
        fprintf(stderr, "ERROR: --output_prefix / -o required.\n");
        exit(1);
    }
    if (kmers_idx == 0){
        fprintf(stderr, "ERROR: --kmers / -k is required\n");
        exit(1);
    }
    if (names_idx == 0){
        fprintf(stderr, "ERROR: --names / -n is required\n");
        exit(1);
    }
    if (kmers_idx != names_idx){
        fprintf(stderr, "ERROR: unequal number of k-mers (-k) and names (-n) provided\n");
        exit(1);
    }
    if (kmers_idx == 1){
        fprintf(stderr, "ERROR: cannot demultiplex with only one species\n");
        exit(1);
    }
    for (int i = 0; i < names_idx; ++i){
        for (int j = 0; j < strlen(names[i]); ++j){
            if (names[i][j] == ' '){
                fprintf(stderr, "ERROR: spaces detected in given species name(s)\n");
                exit(1);
            }
        }
    }
    
    // Load k-mer tables
    int num_tables = kmers_idx;
    Kmer_Stream* tables[num_tables];
    // Ensure all k-mers are same
    int kmer = -1;
    for (int i = 0; i < num_tables; ++i){
        fprintf(stderr, "loading %s\n", kmers[i]);
        tables[i] = Open_Kmer_Stream(kmers[i]);
        if (tables[i] == NULL){
            fprintf(stderr, "ERROR opening %s\n", kmers[i]);
            exit(1);
        }
        else{
            if (kmer != -1 && tables[i]->kmer != kmer){
                fprintf(stderr, "ERROR: conflicting k-mer lengths: %d %d\n", kmer, tables[i]->kmer);
                exit(1);
            }
            kmer = tables[i]->kmer;
            fprintf(stderr, "Loaded %d-mer table %s\n", tables[i]->kmer, kmers[i]);
            
        }
    }
    
    // Store all counts of k-mers, which can later be sorted in descending order
    struct counts counts_tables[num_tables];
    for (int i = 0; i < num_tables; ++i){
        counts_init(&counts_tables[i], 1048576);
    }
    
    // Write out species names
    char namebuf[1024];
    sprintf(&namebuf[0], "%s.names", output_prefix);
    FILE* f = fopen(&namebuf[0], "w");
    for (int i = 0; i < num_tables; ++i){
        fprintf(f, "%s\n", names[i]);
    }
    fclose(f);
    
    // prepare output files
    gzFile outs[num_tables];
    FILE* outs_tmp[num_tables];
    for (int i = 0; i < num_tables; ++i){
        if (num_samp > 0){
            sprintf(&namebuf[0], "%s.%d.tmp", output_prefix, i);
            outs_tmp[i] = fopen(namebuf, "w");
            if (!outs_tmp[i]){
                fprintf(stderr, "ERROR opening %s for writing.\n", &namebuf[0]);
                exit(1);
            }
        }
        sprintf(&namebuf[0], "%s.%d.kmers", output_prefix, i);
        outs[i] = gzopen(namebuf, "w");
        if (!outs[i]){
            fprintf(stderr, "ERROR opening %s for writing.\n", &namebuf[0]);
            exit(1);
        }
    }

    for (int i = 0; i < num_tables; ++i){
        First_Kmer_Entry(tables[i]);
    }

    char kmer_text[kmer+1]; 
    int ntie = 0;
    int ties[num_tables];
    int finished = 0;
    while (!finished){
        // Each iteration: check for unique k-mers
        // Then increment only iterators with lowest sort-order k-mers
        // A k-mer is unique if it's lower sort order than the other streams
        // and does not match others
        int nvalid = 0;
        for (int i = 0; i < num_tables; ++i){
            if (tables[i]->csuf != NULL){
                nvalid++;
            }
        }
        if (nvalid < 2){
            // Not enough k-mers to compare.
            finished = 1;
            break;
        } 
        else{
            ntie = 0;
            int min_idx = 0;
            while (tables[min_idx]->csuf == NULL){
                ++min_idx;
            }
            // Compare all k-mers
            for (int i = min_idx + 1; i < num_tables; ++i){
                if (tables[i]->csuf != NULL){ 
                    int comp = kcomp(tables[i], tables[min_idx]);
                    if (comp < 0){
                        min_idx = i;
                        ntie = 0;
                    }
                    else if (comp == 0){
                        // tie
                        ties[ntie] = i;
                        ntie++;
                    }
                }
            }
            
            // If no tie, print the lowest.
            if (ntie == 0){
                char* b = Current_Kmer(tables[min_idx], &kmer_text[0]);
                int c = Current_Count(tables[min_idx]);
                if (num_samp <= 0){
                    if (print_dat_gz(&outs[min_idx], kmer, b, c, tree, maxdust) == 1){
                    
                    }
                }
                else{
                    if (print_dat(outs_tmp[min_idx], kmer, b, c, tree, maxdust) == 1){\
                        counts_add(&counts_tables[min_idx], c);
                    }
                }
            }
            // Only increment lowest-value iterators.
            int incremented = 0;
            if (tables[min_idx]->csuf != NULL){
                Next_Kmer_Entry(tables[min_idx]);
                incremented = 1;
            }
            for (int i = 0; i < ntie; ++i){
                if (tables[ties[i]]->csuf != NULL){
                    Next_Kmer_Entry(tables[ties[i]]);
                    incremented = 1;
                }   
            }
            if (!incremented){
                // Break out of the loop
                finished = 1;
                break;
            }
        }
    }
    
    // Any remaining entries in the last iterator are unique.
    for (int i = 0; i < num_tables; ++i){
        while (tables[i]->csuf != NULL){
            char* b = Current_Kmer(tables[i], &kmer_text[0]);
            int c = Current_Count(tables[i]);
            if (num_samp <= 0){
                if (print_dat_gz(&outs[i], kmer, b, c, tree, maxdust) == 1){

                }
            }
            else{
                if (print_dat(outs_tmp[i], kmer, b, c, tree, maxdust) == 1){
                    counts_add(&counts_tables[i], c);
                }
            }
            Next_Kmer_Entry(tables[i]);
        }
    }
    
    if (num_samp <= 0){
        for (int i = 0; i < num_tables; ++i){
            gzclose(outs[i]);
        }
    }
    else{
        char linebuf[1024];
        char intbuf[50];
        int intbuf_size = 50;

        // Now sort counts in decreasing order and print the top number specified to each final output file.
        for (int i = 0; i < num_tables; ++i){
            // Stop writing tmp file
            fclose(outs_tmp[i]);

            // Compute a cutoff, where all counts above the value will be chosen,
            // and a sample fraction for k-mers with the count at the cutoff.
            int cutoff = -1;
            double sampfrac = 0.0;
            if (num_samp != -1 && counts_tables[i].nvals > num_samp){
                fprintf(stderr, "Sampling %.2fM k-mers for species %d...\n", (double)num_samp/(double)1e6, i);
                counts_sort(&counts_tables[i]);
                cutoff = -counts_get(&counts_tables[i], num_samp);
                sampfrac = 1.0;
                // Find out how many other k-mers have this same count
                int n_at_val = 1;
                int n_include = 1;
                for (int z = num_samp-1; z >= 0; z--){
                    if (-counts_get(&counts_tables[i], z) == cutoff){
                        n_at_val++;
                        n_include++;
                    }
                    else{
                        break;
                    }
                }
                for (int z = num_samp + 1; z < counts_tables[i].nvals; ++z){
                    if (-counts_get(&counts_tables[i], z) == cutoff){
                        n_at_val++;
                    }
                    else{
                        break;
                    }
                }
                sampfrac = (double)n_include/(double)n_at_val;
                fprintf(stderr, "Choosing k-mers with count above %d plus %.2f%% of k-mers with count = %d\n", cutoff,
                    sampfrac*100.0, cutoff);
            }
            counts_destroy(&counts_tables[i]);

            fprintf(stderr, "Filtering tmp file %d...\n", i);
            // Take all counts passing filter 
            sprintf(&namebuf[0], "%s.%d.tmp", output_prefix, i);
            FILE* inf = fopen(namebuf, "r");
            ssize_t read;
            size_t len = 0;
            while(fgets(&linebuf[0], 1024*sizeof(char), inf)){
                int sep = -1;
                for (int x = 0; x < intbuf_size; ++x){
                    if (linebuf[x] == '\n' || linebuf[x] == '\0'){
                        break;
                    }
                    else if (linebuf[x] == '\t'){
                        intbuf[x] = '\0';
                        sep = x;
                        break;
                    }
                    else{
                        intbuf[x] = linebuf[x];
                    }
                }
                
                if (sep != -1){
                    int count = atoi(intbuf);
                    if (count > cutoff || 
                       (count == cutoff && (double)rand() / (double)RAND_MAX < sampfrac)){
                        // Accept.
                        strncpy(&kmer_text[0], &linebuf[sep+1], kmer);
                        kmer_text[kmer] = '\0';
                        gzwrite(outs[i], kmer_text, kmer);
                        gzwrite(outs[i], "\n", 1);
                    }
                }
            }
            gzclose(outs[i]);
            sprintf(&namebuf[0], "%s.%d.tmp", output_prefix, i);
            remove(namebuf);
            //remove(outs_tmp[i]);
        }
    }
    //free(entries);
    //free(bases);
    
    free_suftree(tree);
    free(tree);
    return 0;

}
