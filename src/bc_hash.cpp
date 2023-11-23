#include <cstdlib>
#include <string.h>
#include <utility>
#include <bitset>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <sys/time.h>
#include "bc_hash.h"

using namespace std;


bool str2bc(const char* str, bc& this_bc, int len){
    this_bc.reset();
    int bcbit = 0;
    for (int i = 0; i < len; ++i){
        switch(str[i]){
            case 'A':
            case 'a':
                bcbit += 2;
                break;
            case 'C':
            case 'c':
                this_bc.set(bcbit);
                bcbit += 2;
                break;
            case 'G':
            case 'g':
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            case 'T':
            case 't':
                this_bc.set(bcbit);
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            default:
                return false;
                break;
        }
    }
    return true;
}


bool str2kmer(const char* str, kmer& this_kmer, int k, int& start, bool& reset){
    if (start > strlen(str) - k){
        // impossible
        start = -1;
        return false;
    }
    if (start == 0 || reset){
        // Fill in the whole thing
        int kmerbit = 0;
        for (int i = start; i < start + k; ++i){
            switch(str[i]){
                case 'a':
                case 'A':
                    // Do nothing
                    kmerbit += 2;
                    break;
                case 'c':
                case 'C':
                    this_kmer.set(kmerbit);
                    kmerbit += 2;
                    break;
                case 'g':
                case 'G':
                    this_kmer.set(kmerbit+1);
                    kmerbit += 2;
                case 't':
                case 'T':
                    this_kmer.set(kmerbit);
                    this_kmer.set(kmerbit+1);
                    kmerbit += 2;
                    break;
                default:
                    // Can't build this k-mer. Instruct method to start at a position where we
                    // will miss this base next time.
                    reset = true;
                    start = i + 1;
                    return false;
                    break;
            }
        }
        start++;
        return true;
    }
    else{
        // Bit-shift and fill in the last one.
        this_kmer >>= 2;
        switch(str[start + k - 1]){
            case 'A':
            case 'a':
                // Set no bits
                break;
            case 'C':
            case 'c':
                this_kmer.set((k-1)*2);
                break;
            case 'G':
            case 'g':
                this_kmer.set((k-1)*2+1);
                break;
            case 't':
            case 'T':
                this_kmer.set((k-1)*2);
                this_kmer.set((k-1)*2 + 1);
                break;
            default:
                start = start + k;
                reset = true;
                return false;
                break;
        }
        start++;
        return true;
    }
}

bool str2bc_rc(const char* str, bc& this_bc, int len){
    this_bc.reset();
    int bcbit = 0;
    for (int i = len-1; i >= 0; --i){
        switch(str[i]){
            case 'T':
            case 't':
                bcbit += 2;
                break;
            case 'G':
            case 'g':
                this_bc.set(bcbit);
                bcbit += 2;
                break;
            case 'C':
            case 'c':
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            case 'A':
            case 'a':
                this_bc.set(bcbit);
                this_bc.set(bcbit+1);
                bcbit += 2;
                break;
            default:
                return false;
                break;
        }
    }
    return true;
}

bool str2kmer_rc(const char* str, kmer& this_kmer, int k, int& start, bool& reset){
    if (start > strlen(str) - k){
        // impossible
        start = -1;
        return false;
    }
    if (start == 0 || reset){
        // Fill in the whole thing
        int kmerbit = 0;
        for (int i = start; i < start + k; ++i){
            switch(str[strlen(str)-1-i]){
                case 't':
                case 'T':
                    // Do nothing
                    kmerbit += 2;
                    break;
                case 'G':
                case 'g':
                    this_kmer.set(kmerbit);
                    kmerbit += 2;
                    break;
                case 'c':
                case 'C':
                    this_kmer.set(kmerbit+1);
                    kmerbit += 2;
                case 'a':
                case 'A':
                    this_kmer.set(kmerbit);
                    this_kmer.set(kmerbit+1);
                    kmerbit += 2;
                    break;
                default:
                    // Can't build this k-mer. Instruct method to start at a position where we
                    // will miss this base next time.
                    reset = true;
                    start = i + 1;
                    return false;
                    break;
            }
        }
        start++;
        return true;
    }
    else{
        // Bit-shift and fill in the last one.
        this_kmer >>= 2;
        switch(str[strlen(str) - 1 - (start + k - 1)]){
            case 'T':
            case 't':
                // Set no bits
                break;
            case 'G':
            case 'g':
                this_kmer.set((k-1)*2);
                break;
            case 'C':
            case 'c':
                this_kmer.set((k-1)*2+1);
                break;
            case 'a':
            case 'A':
                this_kmer.set((k-1)*2);
                this_kmer.set((k-1)*2 + 1);
                break;
            default:
                start = start + k;
                reset = true;
                return false;
                break;
        }
        start++;
        return true;
    }
}
string bc2str(bc& this_bc, int len){
    char strbuf[len+1];
    strbuf[len] = '\0';
    for (int i = 0; i < len; ++i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[i] = 'T';
            }
            else{
                strbuf[i] = 'C';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[i] = 'G';    
            }
            else{
                strbuf[i] = 'A';
            }
        }
    }
    return string(strbuf);
}

string bc2str_rc(bc& this_bc, int len){
    char strbuf[len+1];
    strbuf[len] = '\0';
    int buf_idx = 0;
    for (int i = len-1; i >= 0; --i){
        if (this_bc.test(i*2)){
            if (this_bc.test(i*2+1)){
                strbuf[buf_idx] = 'A';
            }
            else{
                strbuf[buf_idx] = 'G';
            }
        }
        else{
            if (this_bc.test(i*2 + 1)){
                strbuf[buf_idx] = 'C';    
            }
            else{
                strbuf[buf_idx] = 'T';
            }
        }
        buf_idx += 1;
    }
    return string(strbuf);
}

void read_whitelist(string& filename, 
    bcset& bc2idx, 
    vector<string>& bcs,
    multimap<unsigned long, unsigned long>& kmer2bc){
    
    int k = 7;

    ifstream wlfile(filename);
    string line;
    int bclen = 0;
    bc cur_bc;
    unsigned int line_idx = 0;
    bool first = true;
    kmer cur_kmer;
    cur_kmer.reset();

    while(wlfile >> line){
        if (str2bc(line.c_str(), cur_bc, 16)){
            unsigned long hashed = cur_bc.to_ulong();
            //bc2idx.insert(make_pair(cur_bc.to_ulong(), line_idx));
            bc2idx.emplace(cur_bc.to_ulong(), line_idx);
            // Get all k-mers (F orientation only).
            int start = 0;
            bool reset = false;
            kmer cur_kmer;
            bool valid;
            while (start != -1){
                valid = str2kmer(line.c_str(), cur_kmer, k, start, reset);
                if (valid){
                    kmer2bc.insert(make_pair(cur_kmer.to_ulong(), cur_bc.to_ulong())); 
                }
            }        
        }
        bcs.push_back(line);
        line_idx++;
    }
}
void write_fastq(const char* id, 
    int idlen, 
    const char* seq, 
    int seqlen, 
    const char* qual, 
    gzFile& out, 
    const char* bc){

    int buflen = idlen + 3;
    if (bc != NULL){
        buflen += 6 + strlen(bc); 
    }
    // Allow for addition of -1 for 10X-style barcodes
    if (bc != NULL){
        buflen += 2;
    }
    char buf[buflen];
    if (bc != NULL){
        sprintf(buf, "@%s CB:Z:%s-1\n", id, bc);
    }
    else{
        sprintf(buf, "@%s\n", id);
    }
    gzwrite(out, buf, buflen-1);
    char buf2[seqlen + 1];

    sprintf(buf2, "%s\n", seq);
    gzwrite(out, buf2, seqlen + 1);
    
    char buf3[3];
    sprintf(buf3, "+\n");
    gzwrite(out, buf3, 2);
    
    sprintf(buf2, "%s\n", qual);
    gzwrite(out, buf2, seqlen + 1);
}

void revcomp(char* input, char* rc){
    for (int i = 0; i < strlen(input); ++i){
        switch(input[strlen(input)-1-i]){
            case 'A':
            case 'a':
                rc[i] = 'T';
                break;
            case 'C':
            case 'c':
                rc[i] = 'G';
                break;
            case 'G':
            case 'g':
                rc[i] = 'C';
                break;
            case 'T':
            case 't':
                rc[i] = 'A';
                break;
            default:
                rc[i] = 'N';
                break;
        }       
    }
    rc[strlen(input)] = '\0';
}

/**
 *  Given a candidate barcode, checks against the set of all legitimate
 *  barcodes. If found, returns true and sets result_bc_str to a string
 *  representation of the correct barcode.
 */
int match_bc(const char* cur_bc,   
    bool rc, 
    bcset& bc2idx,
    multimap<unsigned long, unsigned long>& kmer2bc){
    static bc bc_binary;
    static char bc_str[17]; // For replacing N characters

    bool success = false;
    if (rc){
        //struct timeval t1;
        //gettimeofday(&t1, NULL);

        if (str2bc_rc(cur_bc, bc_binary, 16)){ 
            /*
            struct timeval t2;
            gettimeofday(&t2, NULL);
            */
            unsigned long ul = bc_binary.to_ulong();
            /*
            struct timeval t3;
            gettimeofday(&t3, NULL);
            */
            if (bc2idx.count(ul) > 0){
                //struct timeval t4;
                //gettimeofday(&t4, NULL);
                //fprintf(stderr, "%ld %ld %ld\n", t2.tv_usec-t1.tv_usec, t3.tv_usec-t2.tv_usec, t4.tv_usec-t3.tv_usec);
                return bc2idx[ul];
            }
            else{
                //struct timeval t4;
                //gettimeofday(&t4, NULL);
                //fprintf(stderr, "%ld %ld %ld\n", t2.tv_usec-t1.tv_usec, t3.tv_usec-t2.tv_usec, t4.tv_usec-t3.tv_usec);
                return -1;
            }
        }
        else{
            // Contains N.
            // Only allow up to one.
            int npos = -1;
            for (int i = 0; i < 16; ++i){
                if (cur_bc[i] == 'N'){
                    if (npos != -1){
                        // Multiple Ns; bail.
                        return -1;
                    }
                    else{
                        npos = i;
                    }
                }
            }
            
            // Try to mutate to every possible letter.
            bc bc_A;
            bc bc_C;
            bc bc_G;
            bc bc_T;

            bool pass_A = false;
            bool pass_C = false;
            bool pass_G = false;
            bool pass_T = false;
            
            strncpy(&bc_str[0], cur_bc, 16);
            bc_str[16] = '\0';

            bc_str[npos] = 'A';
            if (str2bc_rc(bc_str, bc_A, 16) && bc2idx.count(bc_A.to_ulong()) > 0){
                pass_A = true;
            } 
            bc_str[npos] = 'C';
            if (str2bc_rc(bc_str, bc_C, 16) && bc2idx.count(bc_C.to_ulong()) > 0){
                pass_C = true;
            }
            bc_str[npos] = 'G';
            if (str2bc_rc(bc_str, bc_G, 16) && bc2idx.count(bc_G.to_ulong()) > 0){
                pass_G = true;
            }
            bc_str[npos] = 'T';
            if (str2bc_rc(bc_str, bc_T, 16) && bc2idx.count(bc_T.to_ulong()) > 0){
                pass_T = true;
            }
            if (pass_A && !pass_C && !pass_G && !pass_T){
                return bc2idx[bc_A.to_ulong()];
            }
            else if (pass_C && !pass_A && !pass_G && !pass_T){
                return bc2idx[bc_C.to_ulong()];
            }
            else if (pass_G && !pass_A && !pass_C && !pass_T){
                return bc2idx[bc_G.to_ulong()];
            }
            else if (pass_T && !pass_A && !pass_C && !pass_G){
                return bc2idx[bc_T.to_ulong()];
            }
            else{
                return -1;
            }
        }
    }
    else{
        if (str2bc(cur_bc, bc_binary, 16)){
            if (bc2idx.count(bc_binary.to_ulong()) > 0){
                return bc2idx[bc_binary.to_ulong()];
            }
            else{
                return -1;
            }
        }
        else{
            // Contains N.
            // Only allow up to one.
            int npos = -1;
            for (int i = 0; i < 16; ++i){
                if (cur_bc[i] == 'N'){
                    if (npos != -1){
                        // Multiple Ns; bail.
                        return -1;
                    }
                    else{
                        npos = i;
                    }
                }
            }
            
            strncpy(&bc_str[0], cur_bc, 16);
            bc_str[16] = '\0';

            // Try to mutate to every possible letter.
            bc bc_A;
            bc bc_C;
            bc bc_G;
            bc bc_T;

            bool pass_A = false;
            bool pass_C = false;
            bool pass_G = false;
            bool pass_T = false;

            bc_str[npos] = 'A';
            if (str2bc(bc_str, bc_A, 16) && bc2idx.count(bc_A.to_ulong()) > 0){
                pass_A = true;
            } 
            bc_str[npos] = 'C';
            if (str2bc(bc_str, bc_C, 16) && bc2idx.count(bc_C.to_ulong()) > 0){
                pass_C = true;
            }
            bc_str[npos] = 'G';
            if (str2bc(bc_str, bc_G, 16) && bc2idx.count(bc_G.to_ulong()) > 0){
                pass_G = true;
            }
            bc_str[npos] = 'T';
            if (str2bc(bc_str, bc_T, 16) && bc2idx.count(bc_T.to_ulong()) > 0){
                pass_T = true;
            }
            if (pass_A && !pass_C && !pass_G && !pass_T){
                return bc2idx[bc_A.to_ulong()];
            }
            else if (pass_C && !pass_A && !pass_G && !pass_T){
                return bc2idx[bc_C.to_ulong()];
            }
            else if (pass_G && !pass_A && !pass_C && !pass_T){
                return bc2idx[bc_G.to_ulong()];
            }
            else if (pass_T && !pass_A && !pass_C && !pass_G){
                return bc2idx[bc_T.to_ulong()];
            }
            else{
                return -1;
            }
        }
    }
    // If we've made it here, there are no Ns and the barcode does not fully match
    // an existing barcode.
    // Check for k-mer matches.

    return -1;
}

/**
 *  RNA-seq barcodes are at the beginning of R1, in forward
 *  orientation.
 */
int get_bc_rna(const char* seq, 
    int seqlen, 
    bcset& bc2idx){

    // Placeholder
    static multimap<unsigned long, unsigned long> kmer2bc;

    return match_bc(seq, false, bc2idx, kmer2bc); 
}

/**
 *  Given an index read for ATAC-seq data, check the beginning &
 *  end of the read in forward & reverse comp orientation to
 *  find a valid barcode. Returns success/failure
 */
int get_bc_atac(const char* seq, 
    int seqlen, 
    bcset& bc2idx){
 
    // Placeholder
    static multimap<unsigned long, unsigned long> kmer2bc;

    // seq_i should be R2 for ATAC data (index read)
    
    // Count each time we see a barcode in each orientation. Whichever orientation
    // is the most common will be the determinant of which is checked first.
    static int b_f = 0; // beginning of read; forward orientation
    static int b_r = 0; // beginning of read; reverse orientation
    static int e_f = 0; // end of read; forward orientation
    static int e_r = 0; // end of read; reverse orientation

    static int bclen = 16;

    bool success;
    int bc_idx = -1;

    // Try beginning & end of index read; forward & reverse comp orientation
    
    unsigned int corr_bc_idx;

    bool hasN = false;
    
    // After 5,000 reads, stop checking for alternate orientations
    int test_cutoff = 5000;
    bool reached_cutoff = b_f + b_r + e_f + e_r >= test_cutoff;

    if (b_f > (b_r + e_f + e_r) || b_r > (b_f + e_f + e_r)){
        if (b_r > (b_f + e_f + e_r)){
            // Beginning, reverse
            bc_idx = match_bc(seq, true, bc2idx, kmer2bc);
            if (bc_idx != -1){
                b_r++;
            }   
            else if (!reached_cutoff){
                bc_idx = match_bc(seq, false, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    b_f++;   
                }
            }
        }
        else{
            // Beginning, forward
            bc_idx = match_bc(seq, false, bc2idx, kmer2bc);
            if (bc_idx != -1){
                b_f++;
            }
            else if (!reached_cutoff){
                bc_idx = match_bc(seq, true, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    b_r++;
                }
            }
        }
        if (!reached_cutoff && bc_idx == -1){
            // End, forward
            bc_idx = match_bc(seq + (seqlen - bclen), false, bc2idx, kmer2bc);
            if (bc_idx != -1){
                e_f++;
            }
            else{
                // End, reverse
                bc_idx = match_bc(seq + (seqlen - bclen), true, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    e_r++;
                }
            }
        }
    }
    else{
        if (e_f > (b_f + b_r + e_r)){
            // End, forward
            bc_idx = match_bc(seq + (seqlen - bclen), false, bc2idx, kmer2bc);
            if (bc_idx != -1){
                e_f++;
            }
            else if (!reached_cutoff){
                bc_idx = match_bc(seq + (seqlen - bclen), true, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    e_r++;
                }
            }
        }
        else{
            // End, reverse
            bc_idx = match_bc(seq + (seqlen - bclen), true, bc2idx, kmer2bc);
            if (bc_idx != -1){
                e_r++;
            }
            else if (!reached_cutoff){
                // End, forward
                bc_idx = match_bc(seq + (seqlen - bclen), false, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    e_f++;
                }
            }
        }
        if (!reached_cutoff && bc_idx != -1){
            // Beginning, forward
            bc_idx = match_bc(seq, false, bc2idx, kmer2bc);
            if (bc_idx != -1){
                b_f++;
            }
            else{
                // Beginning, reverse
                bc_idx = match_bc(seq, true, bc2idx, kmer2bc);
                if (bc_idx != -1){
                    b_r++;
                }
            }
        }
    }
    return bc_idx;
}
