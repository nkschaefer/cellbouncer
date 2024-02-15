#include <zlib.h>
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
#include <unordered_map>
#include <set>
#include <deque>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <limits.h>
#include <thread>
#include <sys/stat.h>
#include <condition_variable>
#include <mutex>
#include <htswrapper/bc.h>
#include <htswrapper/serialize.h>
#include "common.h"
#include "reads_io.h"
#include "kmsuftree.h"

KSEQ_INIT(gzFile, gzread);

using std::cout;
using std::endl;
using namespace std;

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
 * Only using canonical k-mers = first in sort order relative to reverse complement.
 * Check if a given k-mer or its reverse complement comes first in alphabetical sort order.
 */
bool is_rc_first(char* kmer, int k){
    for (int i = 0; i < k; ++i){
        char char_f = kmer[k];
        switch(kmer[k]){
            case 'a':
            case 'A':
                char_f = 'A';
                break;
            case 'c':
            case 'C':
                char_f = 'C';
                break;
            case 'g':
            case 'G':
                char_f = 'G';
                break;
            case 't':
            case 'T':
                char_f = 'T';
                break;
            default:
                char_f = 'N';
                break;
        }
        char char_rc = 'N';
        switch(kmer[k-1-i]){
            case 'A':
            case 'a':
                char_rc = 'T';
                break;
            case 'c':
            case 'C':
                char_rc = 'G';
                break;
            case 'g':
            case 'G':
                char_rc = 'C';
                break;
            case 't':
            case 'T':
                char_rc = 'A';
                break;
            default:
                char_rc = 'N';
                break;
        }
        if (char_f == 'N' || char_rc == 'N'){
            // No matter; neither will be in the lookup table.
            return false;
        }
        if (char_rc < char_f){
            return true;
        }
        else if (char_f < char_rc){
            return false;
        }   
    }
    return false;
}

/**
 * Write a FASTQ record to disk.
 */
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
        if (str2bc_rc(cur_bc, bc_binary, 16)){ 
            unsigned long ul = bc_binary.to_ulong();
            if (bc2idx.count(ul) > 0){
                return bc2idx[ul];
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

void demux_atac_reads(string& r1filename, 
    string& r2filename, 
    string& r3filename, 
    bcset& bc2species,
    map<short, string>& idx2species,
    string& outdir,
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna,
    long int& atac_reads_tot,
    map<short, long int>& atac_reads_pass){
    
    // Trim directory paths off file names for output file names
    string r1filetrim = filename_nopath(r1filename);
    string r2filetrim = filename_nopath(r2filename);
    string r3filetrim = filename_nopath(r3filename);
    
    // If no "ATAC" identifier present in filenames, add it
    if (r1filetrim.length() < 5 || r1filetrim.substr(0, 5) != "ATAC_"){
        r1filetrim = "ATAC_" + r1filetrim;
    }
    if (r2filetrim.length() < 5 || r2filetrim.substr(0, 5) != "ATAC_"){
        r2filetrim = "ATAC_" + r2filetrim;
    }
    if (r3filetrim.length() < 5 || r3filetrim.substr(0, 5) != "ATAC_"){
        r3filetrim = "ATAC_" + r3filetrim;
    }

    // To keep the 10X pipeline happy, we'll create output files with the same name as
    // input files, separated in directories according to species of origin.
    if (outdir[outdir.length()-1] != '/'){
        outdir += "/";
    }
    gzFile outfps[idx2species.size() * 3];
    for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
        string dirn = outdir + spec->second;
        if (!mkdir(dirn.c_str(), 0775)){
            // Assume directory already exists
        }
        string fn = dirn + "/" + r1filetrim;
        outfps[spec->first * 3] = gzopen(fn.c_str(), "w");
        if (!outfps[spec->first * 3]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r2filetrim;
        outfps[spec->first * 3 + 1] = gzopen(fn.c_str(), "w");
        if (!outfps[spec->first * 3 + 1]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r3filetrim;
        outfps[spec->first * 3 + 2] = gzopen(fn.c_str(), "w");
        if (!outfps[spec->first * 3 + 2]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
    }
    
    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    int i_progress;
    gzFile f_fp;
    gzFile r_fp;
    gzFile i_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    kseq_t* seq_i;
    f_fp = gzopen(r1filename.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1filename.c_str());
        exit(1);
    }    
    r_fp = gzopen(r3filename.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r3filename.c_str());
        exit(1);
    }
    i_fp = gzopen(r2filename.c_str(), "r");
    if (!i_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2filename.c_str());
        exit(1);
    }
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    seq_i = kseq_init(i_fp);

    int bclen = 16;
    char bc_str[bclen+1];
    
    bc cur_bc;
    bc rna_bc;

    vector<string> bc_dummy;

    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R3 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        i_progress = kseq_read(seq_i);
        if (i_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R2 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        
        atac_reads_tot++;

        bool success = false;
        short species = -1;
        if (whitelist_rna.size() > 0 && atac_bc2idx.size() > 0){
            int bc_idx = get_bc_atac(seq_i->seq.s, seq_i->seq.l, atac_bc2idx);
            if (bc_idx != -1){
                // This is an index into the atac whitelist
                if (bc2species.count(whitelist_rna[bc_idx]) > 0){
                    species = bc2species[whitelist_rna[bc_idx]];
                    success = true;
                }       
            } 
        }
        else{
            int bc_idx = get_bc_atac(seq_i->seq.s, seq_i->seq.l, bc2species);
            if (bc_idx != -1){
                species = bc_idx;       
            }
        }
        if (success && species != -1){    
            atac_reads_pass[species]++;
            write_fastq(seq_f->name.s, seq_f->name.l, seq_f->seq.s, seq_f->seq.l, seq_f->qual.s, 
                outfps[species*3], NULL);
            write_fastq(seq_i->name.s, seq_i->name.l, seq_i->seq.s, seq_i->seq.l, seq_i->qual.s,
                outfps[species*3 + 1], NULL);
            write_fastq(seq_r->name.s, seq_r->name.l, seq_r->seq.s, seq_r->seq.l, seq_r->qual.s,
                outfps[species*3 + 2], NULL);                
        }
    }

    // Close output files.
    for (int i = 0; i < idx2species.size(); ++i){
        gzclose(outfps[i*3]);
        gzclose(outfps[i*3 + 1]);
        gzclose(outfps[i*3 + 2]);
    }
}

void demux_rna_reads(string& r1filename, 
    string& r2filename,
    string& file_prefix,
    bcset& bc2species,
    map<short, string>& idx2species,
    string& outdir,
    long int& rna_reads_tot,
    map<short, long int>& rna_reads_pass){
    
    // Trim paths from filenames for output filenames
    string r1filetrim = filename_nopath(r1filename);
    string r2filetrim = filename_nopath(r2filename);

    // Append prefix (i.e. "GEX") identifier to output files, if necessary
    string prefix = file_prefix + "_";
    if (r1filetrim.length() < prefix.length() || r1filetrim.substr(0, prefix.length()) != prefix){
        r1filetrim = prefix + r1filetrim;
    }
    if (r2filetrim.length() < prefix.length() || r2filetrim.substr(0, prefix.length()) != prefix){
        r2filetrim = prefix + r2filetrim;
    }

    // To keep the 10X pipeline happy, we'll create output files with the same name as
    // input files, separated in directories according to species of origin.
    if (outdir[outdir.length()-1] != '/'){
        outdir += "/";
    }
    
    gzFile outfps[2 * idx2species.size()];

    for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
        string dirn = outdir + spec->second;
        if (!mkdir(dirn.c_str(), 0775)){
            // Assume directory already exists
        }
        string fn = dirn + "/" + r1filetrim;
        outfps[spec->first * 2] = gzopen(fn.c_str(), "w");
        if (!outfps[spec->first * 2]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r2filetrim;
        outfps[spec->first * 2 + 1] = gzopen(fn.c_str(), "w");
        if (!outfps[spec->first * 2 + 1]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
    }
    
    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    f_fp = gzopen(r1filename.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1filename.c_str());
        exit(1);
    }    
    r_fp = gzopen(r2filename.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2filename.c_str());
        exit(1);
    }
   
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    
    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R3 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        rna_reads_tot++;

        int bc_idx = get_bc_rna(seq_f->seq.s, seq_f->seq.l, bc2species);
        if (bc_idx != -1){
            int species = bc_idx;
            rna_reads_pass[species]++;
            write_fastq(seq_f->name.s, seq_f->name.l, seq_f->seq.s, seq_f->seq.l, seq_f->qual.s, 
                outfps[species*2], NULL);
            write_fastq(seq_r->name.s, seq_r->name.l, seq_r->seq.s, seq_r->seq.l, seq_r->qual.s,
                outfps[species*2 + 1], NULL);
        }
    }

    // Close output files.
    for (int i = 0; i < idx2species.size(); ++i){
        gzclose(outfps[i*2]);
        gzclose(outfps[i*2 + 1]);
    }
}

/**
 * Read data about whitelisted 10X ATAC / RNA barcodes. Required to map ATAC to RNA barcodes.
 */
void parse_whitelists(string& whitelist_atac_filename,
    string& whitelist_rna_filename, 
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna,
    bcset& rna_bc2idx){

    if (whitelist_atac_filename != ""){
        // Create scope to free temporary data structures
        ifstream wlfile_atac(whitelist_atac_filename);
        int bc_idx = 0;
        string line;
        bc atac_bc;

        while(wlfile_atac >> line){
            if (!str2bc(line.c_str(), atac_bc, 16)){
                fprintf(stderr, "ERROR: invalid barcode %s in ATAC whitelist\n", line.c_str());
                exit(1);
            }
            //atac_bc2idx.insert(make_pair(atac_bc.to_ulong(), bc_idx));
            atac_bc2idx.emplace(atac_bc.to_ulong(), bc_idx);
            bc_idx++;
        }
        // Now map these to RNA barcodes
        ifstream wlfile_rna(whitelist_rna_filename);
        bc_idx = 0;
        while(wlfile_rna >> line){
            if (!str2bc(line.c_str(), atac_bc, 16)){
                fprintf(stderr, "ERROR: invalid barcode %s in RNA-seq whitelist\n", line.c_str());
                exit(1);
            }
            whitelist_rna.push_back(atac_bc.to_ulong());
            rna_bc2idx.emplace(atac_bc.to_ulong(), bc_idx);
            bc_idx++;
        }
    }
    else if (whitelist_rna_filename != ""){
        ifstream wlfile_rna(whitelist_rna_filename);
        int bc_idx = 0;
        string line;
        bc rna_bc;
        while (wlfile_rna >> line){
            if (!str2bc(line.c_str(), rna_bc, 16)){
                fprintf(stderr, "ERROR: invalid barcode %s in RNA whitelist\n", line.c_str());
                exit(1);
            }
            rna_bc2idx.emplace(rna_bc.to_ulong(), bc_idx);
            whitelist_rna.push_back(rna_bc.to_ulong());
            bc_idx++;
        }
    }    
}

void process_gex_files(string& r1filename, 
    string& r2filename,
    readsThreadPool& pool){

    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    f_fp = gzopen(r1filename.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1filename.c_str());
        exit(1);
    }    
    r_fp = gzopen(r2filename.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2filename.c_str());
        exit(1);
    }
   
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R3 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        rp_info params;
        params.seq_f = seq_f->seq.s;
        params.seq_r = seq_r->seq.s;
        
        pool.add_rp_job(params);        
    
    }
}

readsThreadPool::readsThreadPool(int nt, int ns, string& outbase){
    this->terminate_threads = false;
    this->num_threads = nt;
    this->num_species = ns;
    this->initialized = false;
    this->on = false;
    this->out_base = outbase;
}

readsThreadPool::~readsThreadPool(){
    if (this->on){
        this->close_pool();
        this->on = false;
    }
    // Erase all temporary files from disk
    for (int i = 0; i < this->num_threads; ++i){
        string fn = this->out_filenames[i];
        remove(fn.c_str());
    } 
}

void readsThreadPool::init(){
    this->initialized = true;
    this->terminate_threads = false;

}

// Stop all running threads and then destroy them
void readsThreadPool::close_pool(){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->terminate_threads = true;   
    }
    this->has_jobs.notify_all();
    for (int i = 0; i < this->num_threads; ++i){
        this->threads[i].join();
        outstream_close(this->out_files[i]);
    }
    //this->out_files.clear();
    this->threads.clear();
    this->on = false;
}

// Load all data from all output files into the consensus data structure.
void readsThreadPool::load_counts_from_tmp(
    robin_hood::unordered_map<unsigned long, map<short, double> >& bc_species_counts,
    robin_hood::unordered_map<unsigned long, double>& bc_tots){
    
    for (int i = 0; i < this->num_threads; ++i){
        fprintf(stderr, "Loading counts from file %s\n", this->out_filenames[i].c_str());
        instream_info is;
        instream_init(is, this->out_filenames[i], 1024);
        while (!is.finished()){
            unsigned long bc;
            unserialize_ul(is, bc);
            if (bc_tots.count(bc) == 0){
                bc_tots.emplace(bc, 1);
                map<short, double> m;
                bc_species_counts.emplace(bc, m);
            }
            else{
                bc_tots[bc]++;
            }
            for (int j = 0; j < this->num_species; ++j){
                short species_id = (short)j;
                int species_count;
                unserialize_int(is, species_count);
                if (bc_species_counts[bc].count(species_id) == 0){
                    bc_species_counts[bc].insert(make_pair(species_id, species_count));   
                }   
                else{
                    bc_species_counts[bc][species_id] += species_count;
                }
            }   
        }
        fprintf(stderr, "done\n");
    } 
}

/**
 * Add to the queue of jobs for parsing k-mer files
 */
void readsThreadPool::add_kmer_parse_job(string& filename, short species){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->kmer_parse_jobs.push_back(make_pair(filename, species));
    }
    this->has_jobs.notify_one();
}

/**
 * Add to the queue of jobs to read paired reads
 */
void readsThreadPool::add_rp_job(rp_info& info){
    {
        // Add job to processing queue
        unique_lock<mutex> lock(this->queue_mutex);
        this->rp_jobs.push_back(info);
    }
    this->has_jobs.notify_one();
}

/**
 * Add to the queue of jobs to read read triplets
 */
void readsThreadPool::add_rt_job(rt_info& info){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->rt_jobs.push_back(info);
    }
    this->has_jobs.notify_one();
}

void scan_seq_kmers(string& seq,
    kmer_tree_p kt,
    int n_species,
    unsigned long bc_key,
    outstream_info& outfile){
    
    vector<int> species_counts(n_species);
    for (int i = 0; i < n_species; ++i){
        species_counts[i] = 0;
    } 

    for (int i = 0; i < seq.length() - kt->k; ++i){
        int rc = 0;
        if (is_rc_first(&seq[i], kt->k)){
            rc = 1;
        }
        void* dat = kmer_tree_lookup(&seq[i], kt, rc);
        if (dat != NULL){
            short species = *(short*)dat;
            {
                species_counts[species]++;
            }              
        }
    }
    serialize_ul(outfile, bc_key);
    for (int i = 0; i < n_species; ++i){
        serialize_int(outfile, species_counts[i]);
    }    

}

void scan_gex_data(rp_info& params, 
    kmer_tree_p kt, 
    outstream_info& out_file,
    int n_species,
    bcset& rna_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
    int bclen = 16;

    int bc_idx = get_bc_rna(params.seq_f.c_str(), bclen, rna_bc2idx);
    if (bc_idx != -1){
        unsigned long bc_key = whitelist_rna[bc_idx];
        
        // In 10x scRNA-seq, only the reverse read contains information
        // Forward read is barcode
        scan_seq_kmers(params.seq_r, kt, 
            n_species,
            bc_key,
            out_file);

    }
}

void scan_atac_data(rt_info& params, kmer_tree_p kt,
    outstream_info& out_file,
    int n_species, 
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
    int bclen = 16;
    // Get barcode
    int bc_idx = get_bc_atac(params.seq_2.c_str(), params.seq_2.length(), atac_bc2idx);
    if (bc_idx != -1){
        unsigned long bc_key = whitelist_rna[bc_idx];
        
        // R1 and R3 are the forward and reverse reads (R2 is the barcode read)
        scan_seq_kmers(params.seq_1, kt, n_species,
            bc_key,
            out_file);
        scan_seq_kmers(params.seq_3, kt, n_species,
            bc_key,
            out_file);
    }    
}

void parse_kmer_counts(string& countsfilename, 
    short species_idx, 
    kmer_tree_p kt,
    recursive_mutex& ktmutex){

    ifstream infile(countsfilename);
    string line;
    while (infile >> line){
        short* si = (short*) malloc(sizeof(short));
        *si = species_idx;
        {
            lock_guard<recursive_mutex> ktlock(ktmutex);
            kmer_tree_add(line.c_str(), kt, (void*)si, 0); 
        }
    }
}

/** 
 * Non-multithreadable version of above.
 */
void parse_kmer_counts_serial(string& countsfilename,
    short species_idx,
    kmer_tree_p kt){
    
    ifstream infile(countsfilename);
    string line;
    while (infile >> line){
        short* si = (short*) malloc(sizeof(short));
        *si = species_idx;
        kmer_tree_add(line.c_str(), kt, (void*)si, 0); 
    }

}

/** 
 * Worker function for a kmer file parsing job
 */
void readsThreadPool::kmer_parse_thread(kmer_tree_p kt,
    recursive_mutex& ktmutex){
    while (true){
        if (this->kmer_parse_jobs.size() == 0 && this->terminate_threads){
            return;
        }
        pair<string, short> params;
        {
            unique_lock<mutex> lock(this->queue_mutex);
            this->has_jobs.wait(lock, [this]{ return this->kmer_parse_jobs.size() > 0 ||
                this->terminate_threads;});
            if (this->kmer_parse_jobs.size() == 0 && this->terminate_threads){
                return;
            }
            params = this->kmer_parse_jobs[0];
            this->kmer_parse_jobs.pop_front();
        }
        parse_kmer_counts(params.first, params.second, kt, ktmutex);
    }
}

/**
 * Worker function for a GEX read scanning job
 */
 void readsThreadPool::gex_thread(kmer_tree_p kt, 
    outstream_info& out_file,
    bcset& rna_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
     while(true){
        if (this->rp_jobs.size() == 0 && this->terminate_threads){
            return;
        }
        rp_info params;
        {
            unique_lock<mutex> lock(this->queue_mutex);
            this->has_jobs.wait(lock, [this]{ return rp_jobs.size() > 0 ||
                terminate_threads;});
            if (this->rp_jobs.size() == 0 && this->terminate_threads){
                return;
            }
            params = this->rp_jobs[0];
            this->rp_jobs.pop_front();
        }
        scan_gex_data(params, kt, 
            out_file, this->num_species,
            //bc_species_counts, bc_tots, 
            //counts_mutex, tots_mutex, 
            rna_bc2idx, whitelist_rna);
     }
 }

/**
 * Worker function for an ATAC read scanning job
 */
void readsThreadPool::atac_thread(kmer_tree_p kt, 
    outstream_info& out_file,
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
    while(true){
        if (this->rp_jobs.size() == 0 && this->terminate_threads){
            return;
        }
        rt_info params;
        {
            unique_lock<mutex> lock(this->queue_mutex);
            this->has_jobs.wait(lock, [this]{ return rt_jobs.size() > 0 ||
                terminate_threads;});
            if (this->rt_jobs.size() == 0 && this->terminate_threads){
                return;
            }
            params = this->rt_jobs[0];
            this->rt_jobs.pop_front();
        }
        scan_atac_data(params, kt, 
            out_file, this->num_species,
            //bc_species_counts, bc_tots,
            //counts_mutex, tots_mutex, 
            atac_bc2idx, whitelist_rna);
    }
}

void readsThreadPool::launch_kmer_parse_threads(kmer_tree_p kt,
    recursive_mutex& ktmutex){
    if (!this->initialized){
        this->init();
    }
    this->terminate_threads = false;
    for (int i = 0; i < this->num_threads; ++i){
        this->threads.push_back(thread(&readsThreadPool::kmer_parse_thread,
            this,
            kt,
            ref(ktmutex)));
    }
    this->on = true;
}

void readsThreadPool::launch_gex_threads(kmer_tree_p kt, 
    bcset& rna_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
    if (!this->initialized){
        this->init();
    }
    this->terminate_threads = false;

    for (int i = 0; i < this->num_threads; ++i){
        
        // Create output file
        char out_fn_buf[100];
        sprintf(out_fn_buf, "%s.%d.kc", this->out_base.c_str(), i);
        string out_fn = out_fn_buf;
        outstream_init(this->out_files[i], out_fn); 
        this->out_filenames.push_back(out_fn);

        this->threads.push_back(thread(&readsThreadPool::gex_thread,
            this,
            kt,
            ref(this->out_files[i]),
            ref(rna_bc2idx),
            ref(whitelist_rna)));
    }
    this->on = true;
}

void readsThreadPool::launch_atac_threads(kmer_tree_p kt,
    bcset& atac_bc2idx,
    vector<unsigned long>& whitelist_rna){
    
    if (!this->initialized){
        this->init();
    }
    this->terminate_threads = false;
    for (int i = 0; i < this->num_threads; ++i){
        this->threads.push_back(thread(&readsThreadPool::atac_thread,
            this,
            kt, 
            ref(this->out_files[i]),
            ref(atac_bc2idx),
            ref(whitelist_rna)));
    }
    this->on = true;
}


