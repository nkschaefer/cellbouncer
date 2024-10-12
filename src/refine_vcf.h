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
        std::map<std::pair<int, int>, std::string> output_lines;
        std::map<std::pair<int, int>, job_status> output_success;
        
        void launch_threads();
        
        void close_pool();
};

#endif

