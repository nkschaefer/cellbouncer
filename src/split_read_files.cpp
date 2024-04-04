#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "split_read_file [OPTIONS]\n");
    fprintf(stderr, "Given a single fastq file or a pair (F and R), and a number of chunks to create, splits the file into that many chunks, so \
processing can be run more efficiently on a cluster.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "    --r1 -1 The input file for forward, paired reads\n");
    fprintf(stderr, "    --r2 -2 The input file for reverse, paired reads\n");
    fprintf(stderr, "    --single -s The input file for unpaired reads\n");
    fprintf(stderr, "    --output_directory -o The output directory for split files\n");
    fprintf(stderr, "    --num_chunks -n The number of smaller read files to create\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

KSEQ_INIT(gzFile, gzread);

void write_fastq(kseq_t* seq, gzFile& out){
    int buflen = seq->name.l + 3;
    if (seq->comment.l > 0){
        buflen += seq->comment.l + 1;
    }
    char buf[buflen];
    if (seq->comment.l > 0){
        sprintf(buf, "@%s %s\n", seq->name.s, seq->comment.s);
    }
    else{
        sprintf(buf, "@%s\n", seq->name.s);
    }
    gzwrite(out, buf, buflen-1);
    char buf2[seq->seq.l + 1];
    sprintf(buf2, "%s\n", seq->seq.s);
    gzwrite(out, buf2, seq->seq.l + 1);
    char buf3[3];
    sprintf(buf3, "+\n");
    gzwrite(out, buf3, 2);
    sprintf(buf2, "%s\n", seq->qual.s);
    gzwrite(out, buf2, seq->qual.l + 1);
}

string filename_noext(const string& filename){
    size_t pos = filename.rfind(".fastq");
    if (pos == string::npos){
        pos = filename.rfind(".fq");
    }
    if (pos == string::npos){
        pos = filename.rfind(".gz");
    }
    if (pos != string::npos){
        return filename.substr(0, pos);
    }
    return filename;
}

int main(int argc, char *argv[]) {    
    
    /** Define arguments 
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html#Getopt-Long-Options
     * http://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     * Fields for each argument: name, has_arg (values: no_argument, required_argument,
     *     optional_argument)
     * flag = int value to store flag for the option, or NULL if option is string
     * val = short name for string option, or NULL
     */
     
    static struct option long_options[] = {
       {"r1", required_argument, 0, '1'},
       {"r2", required_argument, 0, '2'},
       {"single", required_argument, 0, 's'},
       {"output_directory", required_argument, 0, 'o'},
       {"num_chunks", required_argument, 0, 'n'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string r1file;
    string r2file;
    bool has_r1 = false;
    bool has_r2 = false;
    bool has_paired = false;
    string sfile;
    bool has_single = false;
    string output_directory;
    bool has_output_directory = false;
    int num_chunks = -1;

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "1:2:s:o:n:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case '1':
                r1file = optarg;
                has_r1 = true;
                break;
            case '2':
                r2file = optarg;
                has_r2 = true;
                break;
            case 's':
                sfile = optarg;
                has_single = true;
                break;
            case 'o':
                output_directory = optarg;
                has_output_directory = true;
                break;
            case 'n':
                num_chunks = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    has_paired = has_r1 && has_r2;
    if (!has_paired && !has_single){
        fprintf(stderr, "ERROR: at least one of either -1 and -2 or -s must be provided.\n");
        exit(1);
    }
    
    if (num_chunks <= 0){
        fprintf(stderr, "ERROR: num chunks must be a positive integer\n");
        exit(1);
    }    
    
    if (!has_output_directory){
        fprintf(stderr, "ERROR: output_directory / -o required\n");
        exit(1);
    }
    
    if (output_directory[output_directory.size()-1] == '/'){
        output_directory = output_directory.substr(0, output_directory.length()-1);
    }

    if (!mkdir(output_directory.c_str(), 0775)){
        // Assume directory already exists
    }
    
    // Define output files
    gzFile outfiles[num_chunks*2 + 1];
    
    string base1;
    string base2;
    string base_single;
    
    if (has_paired){
        base1 = filename_nopath(r1file);
        base1 = filename_noext(base1);
        base2 = filename_nopath(r2file);
        base2 = filename_noext(base2);
    }
    else{
        base_single = filename_nopath(sfile);
        base_single = filename_noext(base_single);
    }

    for (int i = 0; i < num_chunks; ++i){
        if (has_paired){
            char fn1[150];
            char fn2[150];
            sprintf(&fn1[0], "%s/%s.%d.fastq.gz", output_directory.c_str(), base1.c_str(), i+1);
            sprintf(&fn2[0], "%s/%s.%d.fastq.gz", output_directory.c_str(), base2.c_str(), i+1);
            outfiles[i*2] = gzopen(fn1, "w");
            outfiles[i*2+1] = gzopen(fn2, "w");
            if (!outfiles[i*2]){
                fprintf(stderr, "ERROR opening %s for writing.\n", fn1);
                exit(1);
            }
            if (!outfiles[i*2+1]){
                fprintf(stderr, "ERROR opening %s for writing.\n", fn2);
                exit(1);
            }
        }
        else{
            char fn[150];
            sprintf(&fn[0], "%s/%s.%d.fastq.gz", output_directory.c_str(), base_single.c_str(),
                i+1);
            outfiles[i*2] = gzopen(fn, "w");
            if (!outfiles[i*2]){
                fprintf(stderr, "ERROR opening %s for writing.\n", fn);
                exit(1);
            }
        }
    }

    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    gzFile s_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    if (has_paired){
        f_fp = gzopen(r1file.c_str(), "r");
        if (!f_fp){
            fprintf(stderr, "ERROR opening %s for reading\n", r1file.c_str());
            exit(1);
        }    
        r_fp = gzopen(r2file.c_str(), "r");
        if (!r_fp){
            fprintf(stderr, "ERROR opening %s for reading\n", r2file.c_str());
            exit(1);
        }
        seq_f = kseq_init(f_fp);
        seq_r = kseq_init(r_fp);
    }

    else{
        s_fp = gzopen(sfile.c_str(), "r");
        if (!s_fp){
            fprintf(stderr, "ERROR opening %s for reading\n", sfile.c_str());
            exit(1);
        }
        seq_f = kseq_init(s_fp);
    }
    
    int cur_chunk_idx = 0;
    while ((f_progress = kseq_read(seq_f)) >= 0){
        
        write_fastq(seq_f, outfiles[cur_chunk_idx*2]);

        if (has_paired){
            r_progress = kseq_read(seq_r);
            if (r_progress < 0){
                fprintf(stderr, "ERROR: read order no longer matching at seq %s in R1 file\n", seq_f->name.s);
                exit(1);
            }    
            write_fastq(seq_r, outfiles[cur_chunk_idx*2+1]);
        }
        
        cur_chunk_idx++;

        if (cur_chunk_idx >= num_chunks){
            cur_chunk_idx = 0;
        }
    }

    kseq_destroy(seq_f);

    // Close output files.
    if (has_paired){
        kseq_destroy(seq_r);
        for (int i = 0; i < num_chunks; ++i){
            gzclose(outfiles[i*2]);
            gzclose(outfiles[i*2 + 1]);
        }
    }
    else{
        for (int i = 0; i < num_chunks; ++i){
            gzclose(outfiles[i*2]);
        }
    }    
    return 0;
}
