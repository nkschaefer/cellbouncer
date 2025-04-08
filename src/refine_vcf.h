#ifndef _CELLBOUNCER_REFINE_VCF_H
#define _CELLBOUNCER_REFINE_VCF_H
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
#include <utility>
#include <math.h>
#include <float.h>
#include <zlib.h>
#include <mutex>
#include <htswrapper/robin_hood/robin_hood.h>

struct job_dat{
    std::map<int, std::pair<float, float> >* data;
    int tid;
    int pos;
    var* v;

    job_dat(std::map<int, std::pair<float, float> >* d, int c, int p, var* vv){
        this->data = d;
        this->tid = c;
        this->pos = p;
        this->v = vv;
    }
};

struct job_status{
    bool finished;
    bool rm;
    bool changed;
};

struct vcf_line{
    int tid;
    int pos;
    char ref;
    char alt;
    double qual;
    std::vector<int> gts;
   
    vcf_line(){
        this->tid = -1;
        this->pos = -1;
        this->qual = 0;
        this->ref = 'N';
        this->alt = 'N';
    } 
    vcf_line(int t, int p, int q, char r, char a){
        this->tid = t;
        this->pos = p;
        this->qual = q;
        this->ref = r;
        this->alt = a;
    }
    void add_gt(int g){
        this->gts.push_back(g);
    }
    void write_record(htsFile* outf, bcf_hdr_t* header, bcf1_t* record){
        if (tid < 0 || pos < 0){
            return;
        }
        // Update chrom/pos
        record->rid = this->tid;
        record->pos = this->pos;
        record->qual = this->qual;

        char alleles[4];
        alleles[0] = this->ref;
        alleles[1] = ',';
        alleles[2] = this->alt;
        alleles[3] = '\0';
        
        // set alleles
        bcf_update_alleles_str(header, record, &alleles[0]);

        // update GTs
        int32_t* gts2 = (int*)malloc(this->gts.size()*2*sizeof(int));
        for (int i = 0; i < this->gts.size(); ++i){
            if (gts[i] == 0){
                gts2[i*2] = bcf_gt_unphased(0);
                gts2[i*2+1] = bcf_gt_unphased(0);
            }
            else if (gts[i] == 1){
                gts2[i*2] = bcf_gt_unphased(0);
                gts2[i*2+1] = bcf_gt_unphased(1);
            }
            else if (gts[i] == 2){
                gts2[i*2] = bcf_gt_unphased(1);
                gts2[i*2+1] = bcf_gt_unphased(1);
            }
            else{
                gts2[i*2] = bcf_gt_missing;
                gts2[i*2+1] = bcf_gt_missing;
            }
        }
        bcf_update_genotypes(header, record, gts2, gts.size()*2);
        free(gts2);
        int ret = bcf_write1(outf, header, record);
    }
};


// This class exists solely to make re-genotyping possible in
// multiple threads (it handles thread management).

class regenotyper{
    private:
        
        int n_samples;
        std::map<int, std::string> tid2chrom;
        std::map<int, double> weights;

        bool threads_init;
        bool terminate_threads;
        int nthread;
        std::mutex queue_mutex;
        std::deque<job_dat> jobs;
        
        double p_thresh;

        std::condition_variable has_jobs;
        std::vector<std::thread> threads;

        void create_threads();
        void worker();


    public:
        regenotyper(int ns, std::map<int, std::string>& t2c, std::map<int, double>& w, int nt, double p);
        
        void add_job(int tid, int pos, var* v, map<int, pair<float, float> >* m);
        std::mutex output_mutex;
        std::map<std::pair<int, int>, vcf_line> output_lines;
        std::map<std::pair<int, int>, job_status> output_success;
        
        void launch_threads();
        
        void close_pool();
};

#endif

