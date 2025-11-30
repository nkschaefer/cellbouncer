#ifndef _CELLBOUNCER_DOWNSAMPLE_VCF_H
#define _CELLBOUNCER_DOWNSAMPLE_VCF_H
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
#include <deque>
#include <list>

// Make bitsets hashable

namespace std {
    template <>
    struct hash<std::bitset<NBITS>> {
        size_t operator()(const std::bitset<NBITS>& bs) const {
            const size_t numBlocks = NBITS / 64 + 1;
            size_t h = 0;

            for (size_t i = 0; i < numBlocks; i++) {

                uint64_t block = 0;
                size_t base = i * 64;

                for (size_t bit = 0; bit < 64; ++bit) {
                    size_t idx = base + bit;
                    if (idx >= NBITS) break;
                    if (bs.test(idx))
                        block |= (uint64_t(1) << bit);
                }

                // Same mixing step you originally had:
                h ^= std::hash<uint64_t>{}(block)
                     + 0x9e3779b97f4a7c15ULL
                     + (h << 6)
                     + (h >> 2);
            }
            return h;
        }
    };
}

namespace std {
    template <>
    struct hash<std::pair<std::bitset<NBITS>, std::bitset<NBITS>>> {
        size_t operator()(const std::pair<std::bitset<NBITS>, std::bitset<NBITS>>& p) const {
            size_t h1 = std::hash<std::bitset<NBITS>>{}(p.first);
            size_t h2 = std::hash<std::bitset<NBITS>>{}(p.second);
            // Combine hashes (Boost's method)
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };
}

void pr(const std::bitset<NBITS>& cl){
    for (int i = 0; i < NBITS; ++i){
        if (cl.test(i)){
            fprintf(stderr, " %d", i);
        }
    }
    fprintf(stderr, "\n");
}

// Class to store nodes in the hierarchy of clades

// This is too memory-inefficient for use currently, as it creates every possible path and
// thus every possible bitset leading from one to another. 

// This was abandoned in favor of limiting the search space based on sizes of bitsets.

class cladeNode{
    private:
        std::map<int, cladeNode*> nodes;
        int N;
    public:
        bool visited;
        bool observed;
        std::bitset<NBITS> clade;

        cladeNode(int n){
            N = n;
            observed = false;
            visited = false;
        }
        ~cladeNode(){
        }

        void add(const std::bitset<NBITS>& cl, 
            std::bitset<NBITS>& cl_progress,
            std::bitset<NBITS>& cl_path,
            int depth,
            std::unordered_map<std::bitset<NBITS>, cladeNode* >& by_clade,
            std::deque<cladeNode>& nodes_all,
            bool verbose = false){
            
            if (visited){
                return;
            }
            visited = true;

            bool addedNode = false;
            for (int i = 0; i < N; ++i){
                if (cl_progress.test(i)){
                    bitset<NBITS> cl2 = cl_progress;
                    cl2.reset(i);
                    bitset<NBITS> cl3 = cl_path;
                    cl3.set(i);
                    if (nodes.count(i) > 0){
                        if (verbose){
                            for (int x = 0; x < depth+1; ++x){
                                fprintf(stderr, "  ");
                            }
                            fprintf(stderr, "EXISTS1 %p DEPTH %d CL ", nodes[i], depth+1);
                            pr(nodes[i]->clade);
                        }
                        nodes[i]->add(cl, cl2, cl3, depth+1, by_clade, nodes_all, verbose);
                    }
                    else if (by_clade.count(cl3) > 0){
                        cladeNode* n = by_clade[cl3];
                        nodes.emplace(i, n);
                        if (verbose){
                            for (int x = 0; x < depth+1; ++x){
                                fprintf(stderr, "  ");
                            }
                            fprintf(stderr, "EXISTS2 %p DEPTH %d CL ", n, depth+1);
                            pr(n->clade);
                        }
                        n->add(cl, cl2, cl3, depth+1, by_clade, nodes_all, verbose);
                    }
                    else{
                        nodes_all.emplace_back(N);
                        cladeNode* newnode = &(*nodes_all.rbegin());
                        newnode->observed = false;
                        newnode->visited = false;
                        nodes.emplace(i, newnode);
                        by_clade.emplace(cl3, newnode);
                        newnode->clade = cl3;
                        if (verbose){
                            for (int x = 0; x < depth+1; ++x){
                                fprintf(stderr, "  ");
                            }
                            fprintf(stderr, "INS %p DEPTH %d CL ", newnode, depth+1);
                            pr(cl3);
                        }
                        newnode->add(cl, cl2, cl3, depth+1, by_clade, nodes_all, verbose);
                    }
                    addedNode = true;
                }
            }
            if (!addedNode){
                // This is the end of the path.
                this->observed = true;
                this->clade = cl;
                if (verbose){
                    for (int x = 0; x < depth; ++x){
                        fprintf(stderr, "  ");
                    }
                    fprintf(stderr, "STORE %p CL ", this);
                    pr(this->clade);
                }
            }
            visited = false;
        }
        
        void traverse_parents(std::vector<std::bitset<NBITS>* >& results, 
            std::deque<cladeNode* >& path, bool verbose=false){
            if (this->visited){
                return;
            }
            else{
                if (verbose){
                    fprintf(stderr, "TRAVERSE ");
                    pr(this->clade);
                }
                this->visited = true;
                path.push_back(this);
                if (this->observed){
                    results.push_back(&this->clade);
                }
                for (typename std::map<int, cladeNode* >::iterator n = nodes.begin(); n != nodes.end(); ++n){
                    n->second->traverse_parents(results, path, verbose);
                }
            }
        }
};

// Class to store a hierarchy of bitset "clades"
class clademap{
    private:
        std::unordered_map<std::bitset<NBITS>, cladeNode* > by_clade;
        std::deque<cladeNode> nodes_all;
        int N;
        bool verbose;
    public:
        
        static bool isparent(const std::bitset<NBITS>& cl1, const std::bitset<NBITS>& cl2){
            if (cl2.count() > cl1.count()){
                if ((cl2 & cl1).count() >= cl1.count()){
                    return true;
                }
            }
            return false;
        }
        
        clademap(int n){
            N = n;
            // Create root
            nodes_all.emplace_back(n);
            verbose = false;
        }

        clademap(){
            N = NBITS;
            nodes_all.emplace_back(NBITS);
            verbose = false;
        }

        ~clademap(){
            by_clade.clear();
            nodes_all.clear();
        }

        void set_verbose(){
            verbose = true;
        }
        
        void insert(const std::bitset<NBITS>& cl){
            if (by_clade.count(cl) > 0){
                by_clade[cl]->observed = true;
            }
            else{
                bitset<NBITS> cpy = cl;
                bitset<NBITS> clpath;
                if (verbose){
                    fprintf(stderr, "(INSERT) ");
                    pr(cl);
                }
                nodes_all[0].add(cl, cpy, clpath, 0, by_clade, nodes_all, verbose);
            }
        }
        
        void lookup(const std::bitset<NBITS>& cl, std::vector<std::bitset<NBITS>* >& results){
            // Must start from exact clade, then look for parents
            if (by_clade.count(cl) > 0){
                if (verbose){
                    fprintf(stderr, "looked up ");
                    pr(cl);
                }
                std::deque<cladeNode* > path;
                by_clade[cl]->traverse_parents(results, path, verbose);
                for (typename std::deque<cladeNode* >::iterator p = path.begin(); p != path.end(); ++p){
                    (*p)->visited = false;
                }
                path.clear();
            }
            else{
                return;
            }
        }
};

#endif
