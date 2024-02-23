#ifndef _CELLID_READS_DEMUX_H
#define _CELLID_READS_DEMUX_H
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
#include <cstdlib>
#include <htslib/kseq.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"

class reads_demuxer{
    private:
        // For barcode matching
        bc_whitelist whitelist;
        
        std::string file_prefix;

        // Filenames of reads
        std::string r1;
        std::string r2;
        std::string r3;
        
        // Is this ATAC-seq?
        bool is_atac;

        // Map barcode to species
        robin_hood::unordered_map<unsigned long, short> bc2species;
        
        // Map species to name
        std::map<short, std::string> idx2species;
        
        // Map species to output file
        std::vector<gzFile> outfiles;

        // Output directory
        std::string outdir;
        
        bool initialized;
        
        void close();
        void init_rna_or_custom(std::string prefix, std::string& r1, std::string& r2);

    public:
        
        reads_demuxer(bc_whitelist& whitelist,
            robin_hood::unordered_map<unsigned long, short>& bc2species,
            std::map<short, std::string>& idx2species,
            std::string& outdir);

        ~reads_demuxer();

        void init_atac(std::string& r1, std::string& r2, std::string& r3);
        void init_rna(std::string& r1, std::string& r2);
        void init_custom(std::string prefix, std::string& r1, std::string& r2);

        bool scan_atac();
        bool scan_rna();
        bool scan_custom();

        void write_fastq(const char* id, int idlen, const char* seq, int seqlen,
            const char* qual, int out_idx);
};
#endif

