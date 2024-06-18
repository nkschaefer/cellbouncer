#! /usr/bin/env python3
import sys
import os
import gzip
import argparse
from collections import Counter
import re
"""
Given a set of 3 MEX-format single cell matrix files, finds all unique
library names appended to cell barcodes and splits the MEX-format input 
into three new output files per library.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mtx", "-m", help="input MTX file", required=True)
    parser.add_argument("--features", "-f", help="input features/genes file", required=True)
    parser.add_argument("--barcodes", "-b", help="input barcodes file", required=True)
    parser.add_argument("--output_prefix", "-o", help="Base name for output files", required=True)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    libs = set([])
    libcounts = Counter()
    lib2idx = {}
    idx2lib = {}
    idx2idx = {}

    f = None
    bc_gz = False
    if options.barcodes[-3:] == '.gz':
        bc_gz = True
        f = gzip.open(options.barcodes, 'r')
    else:
        f = open(options.barcodes, 'r')
    
    bc_outs = []

    bc_idx = 1
    for line in f:
        if bc_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        
        lib = None
        bc = None

        match = re.match(r'^([A-Za-z0-9_\-\.]+)(_|-)([ACGT]+)$', line)
        if match is not None:
            lib = match.group(1)
            bc = match.group(3)
        else:
            match = re.match(r'^([ACGT]+)(_|-)([A-Za-z0-9_\-\.]+)$', line)
            if match is not None:
                lib = match.group(3)
                bc = match.group(1)
        
        if lib is not None:
            
            libidx = None
            
            if lib in libs:
                libidx = lib2idx[lib]
            else:
                libidx = len(libs)
                lib2idx[lib] = libidx
                bcf = gzip.open('{}_{}.barcodes.tsv.gz'.format(options.output_prefix, \
                    lib), 'wt')
                bc_outs.append(bcf)
                libs.add(lib)
            print(bc, file=bc_outs[libidx])
            idx2idx[bc_idx] = libcounts[lib] + 1
            libcounts[lib] += 1
            idx2lib[bc_idx] = lib
        
        bc_idx += 1
    
    f.close()
    
    for bcf in bc_outs:
        bcf.close()
    
    f = None
    features_gz = False
    if options.features[-3:] == '.gz':
        features_gz = True
        f = gzip.open(options.features, 'r')
    else:
        f = open(options.features, 'r')
    
    features_outs = []
    for lib in libs:
        ff = gzip.open('{}_{}.features.tsv.gz'.format(options.output_prefix, lib), 'wt')
        features_outs.append(ff)

    nfeatures = 0
    for line in f:
        if features_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        
        for outf in features_outs:
            print(line, file=outf)
        nfeatures += 1

    f.close()

    for ff in features_outs:
        ff.close()
    
    mtx_outs = {}
    mtx_out_lines = {}
    for lib in libs:
        mtx_outs[lib] = \
            gzip.open('{}_{}.matrix.mtx.gz'.format(options.output_prefix, lib), 'wt')
        mtx_out_lines[lib] = []

    # Print MTX header to each file
    for lib in mtx_outs:
        print("%%MatrixMarket matrix coordinate integer general", file=mtx_outs[lib])
        print("%", file=mtx_outs[lib])

    f = None
    mtx_gz = False
    if options.mtx[-3:] == '.gz':
        f = gzip.open(options.mtx, 'r')
        mtx_gz = True
    else:
        f = open(options.mtx, 'r')
    
    
    first = True
    bc_col1 = True
    for line in f:
        if mtx_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] == '%':
            continue
        else:
            if first:
                nums = line.split(' ')
                if int(nums[0]) == bc_idx - 1:
                    bc_col1 = True
                elif int(nums[1]) == bc_idx - 1:
                    bc_col1 = False
                else:
                    print("ERROR: MTX header does not match the number of barcodes.", file=sys.stderr)
                    exit(1)
                first = False
            else:
                i, j, count = line.split()
                i = int(i)
                j = int(j)
                count = float(count)
                if bc_col1:
                    mtx_out_lines[idx2lib[i]].append("{} {} {:.2f}".format(idx2idx[i], j, count))    
                else:
                    mtx_out_lines[idx2lib[j]].append("{} {} {:.2f}".format(i, idx2idx[j], count))
    
    f.close()

    for lib in mtx_outs:
        if bc_col1:
            print("{} {} {}".format(libcounts[lib], nfeatures, len(mtx_out_lines[lib])), \
                file=mtx_outs[lib])
        else:
            print("{} {} {}".format(nfeatures, libcounts[lib], len(mtx_out_lines[lib])), \
                file=mtx_outs[lib])

        for line in mtx_out_lines[lib]:
            print(line, file=mtx_outs[lib])

        mtx_outs[lib].close()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
