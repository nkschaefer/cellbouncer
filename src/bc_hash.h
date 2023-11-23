#ifndef BCHASH_H
#define BCHASH_H
#include <utility>
#include <cstdlib>
#include <bitset>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <zlib.h>
#include "robin_hood.h"

struct ulong_hash_func {
    // https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
    unsigned long operator()(const unsigned long& x) const{
        unsigned long y = ((x >> 16) ^ x) * 0x45d9f3b;
        y = ((y >> 16) ^ y) * 0x45d9f3b;
        y  = (y >> 16) ^ y;
        return y;
    }
};

struct ulong_equal_func {
    bool operator()( const unsigned long& u1, const unsigned long& u2) const{
        return u1 == u2;
    }
};

//typedef std::unordered_map<unsigned long, unsigned int, ulong_hash_func, ulong_equal_func> bcset;
//typedef std::unordered_map<unsigned long, int> bcset;
typedef robin_hood::unordered_map<unsigned long, int> bcset;

typedef std::bitset<32> bc;
typedef std::bitset<14> kmer;

bool str2bc(const char* string, bc& this_bc, int len);
bool str2bc_rc(const char* string, bc& this_bc, int len);

bool str2kmer(const char* string, kmer& this_kmer, int k, int& start, bool& reset);
bool str2kmer_rc(const char* string, kmer& this_kmer, int k, int& start, bool& reset);

std::string bc2str(bc& this_bc, int len);
std::string bc2str_rc(bc& this_bc, int len);

void read_whitelist(std::string& filename,
    bcset& bc2idx, 
    std::vector<std::string>& bcs,
    std::multimap<unsigned long, unsigned long>& kmer2bc);

void revcomp(char* str, char* rc);

void write_fastq(const char* id, int idlen, const char* seq, int seqlen, const char* qual, gzFile& out, const char* bc);

int get_bc_rna(const char* seq, int seqlen, bcset& bc2idx);
int get_bc_atac(const char* seq, int seqlen, bcset& bc2idx);

#endif
