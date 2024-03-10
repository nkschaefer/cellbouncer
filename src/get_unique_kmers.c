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
    fprintf(stderr, "After running FASTK to count k-mers on different species' transcriptomes,\n");
    fprintf(stderr, "    run this prgram to read FASTK's k-mer table files and output lists of\n");
    fprintf(stderr, "    k-mers unique to each species. These can then be used to demultiplex\n");
    fprintf(stderr, "    reads from multiple species without mapping to a reference genome.\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "--kmers -k Two or more FASTK tables from reference genomes (specify -k\n");
    fprintf(stderr, "    multiple times). These are the files ending with .ktab\n");
    fprintf(stderr, "--names -n The names of the species for which you have provided k-mer tables,\n");
    fprintf(stderr, "   in the same order (specify -n multiple names). Names should not contain spaces\n");
    fprintf(stderr, "   or special characters.\n");
    fprintf(stderr, "--output_prefix -o Where to write species-specific kmers. Files will be\n");
    fprintf(stderr, "    in the format {outprefix}.{index}.kmers (gzip compressed).\n"); 
    exit(code);
}

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

void print_dat(gzFile* outf, char* kmer){
    
    for (int i = 0; i < strlen(kmer); ++i){
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
    gzwrite(*outf, &kmer[0], strlen(kmer));
    gzwrite(*outf, "\n", 1);
}

int main(int argc, char* argv[]){
    
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
    
    int option_index = 0;
    int ch;
    if (argc == 1){
        help(0);
    }
    
    // Parse arguments
    while ((ch = getopt(argc, argv, "k:n:o:h")) != -1){
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
        for (int j = 0; j < strlen(names[names_idx]); ++j){
            if (names[names_idx][i] == ' '){
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

    // Write out species names
    char namebuf[500];
    sprintf(&namebuf[0], "%s.names", output_prefix);
    FILE* f = fopen(&namebuf[0], "w");
    for (int i = 0; i < num_tables; ++i){
        fprintf(f, "%s\n", names[i]);
    }
    fclose(f);

    // prepare output files
    gzFile outs[num_tables];
    for (int i = 0; i < num_tables; ++i){
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
    char kmer_text[150]; 
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
                if (c == 1){
                    print_dat(&outs[min_idx], b);
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
            print_dat(&outs[i], b);
            Next_Kmer_Entry(tables[i]);
        }
    }
    for (int i = 0; i < num_tables; ++i){
        gzclose(outs[i]);
    }
    //free(entries);
    //free(bases);
    return 0;

}
