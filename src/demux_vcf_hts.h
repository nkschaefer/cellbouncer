#ifndef _DEMUX_VCF_HTS_H
#define _DEMUX_VCF_HTS_H
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
#include <sys/stat.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htswrapper/bam.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * Structure to represent population-level genotype information at a 
 * single SNP.
 */
struct var{
    char ref;
    char alt;
    
    // This limits us to 500 individuals in the input VCF
    bitset<500> haps1;
    bitset<500> haps2;
    bitset<500> haps_covered;
    vector<float> gqs;
    float vq;
    var(){
        this->ref = 'N';
        this->alt = 'N';
        this->haps1.reset();
        this->haps2.reset();
        this->haps_covered.reset();
        this->vq = 0.0;
    }
    var(const var& v){
        this->ref = v.ref;
        this->alt = v.alt;
        this->haps1 = v.haps1;
        this->haps2 = v.haps2;
        this->haps_covered = v.haps_covered;
        this->gqs = v.gqs;
        this->vq = v.vq;
    }
};

int read_vcf(std::string& filename, 
    bam_reader& reader,
    std::vector<std::string>& samples,
    std::map<int, std::map<int, var> >& snps,
    int min_vq,
    bool hdr_only,
    bool skip_seq2tid,
    bool allow_missing=true);
 
void get_conditional_match_fracs(std::map<int, std::map<int, var> >& snpdat,
    std::map<std::pair<int, int>, std::map<int, float> >& conditional_match_fracs, 
    int n_samples);

void process_bam_record(bam_reader& reader,
    int snppos,
    var& vardat,
    std::map<int, robin_hood::unordered_map<unsigned long, 
        std::pair<float, float> > >& var_counts,
    bool has_bc_list,
    std::set<unsigned long>& bcs_valid);

void process_bam_record_bulk(bam_reader& reader,
    int snppos,
    var& vardat,
    std::map<int, std::map<int, std::pair<float, float> > >& snp_ref_alt,
    std::map<int, std::map<int, float> >& snp_err);

void dump_vcs_counts(robin_hood::unordered_map<unsigned long, 
        std::pair<float, float> >& varcounts_site,
    robin_hood::unordered_map<unsigned long, std::map<std::pair<int, int>, 
        std::map<std::pair<int, int>, std::pair<float, float> > > >& indv_allelecounts,
    var& snpdat,
    int n_samples);


 


#endif 
