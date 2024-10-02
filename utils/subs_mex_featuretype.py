#! /usr/bin/env python3
import sys
import os
import gzip
import argparse
from collections import Counter
import re
"""
Given a set of 3 MEX-format single cell matrix files containing multiple
data types, subsets to a single chosen data/feature type.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mtx", "-m", help="input MTX file", required=True)
    parser.add_argument("--features", "-f", help="input features/genes file", required=True)
    parser.add_argument("--barcodes", "-b", help="input barcodes file", required=True)
    parser.add_argument("--type", "-t", help="Feature type to retain", required=False, default="Gene Expression")
    parser.add_argument("--output_dir", "-o", help="Directory for output files", required=True)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    if options.output_dir[-1] == '/':
        options.output_dir = options.output_dir[0:-1]

    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)

    f = None
    bc_gz = False
    if options.barcodes[-3:] == '.gz':
        bc_gz = True
        f = gzip.open(options.barcodes, 'r')
    else:
        f = open(options.barcodes, 'r')
    
    bc_idx = 1
    bc_out = gzip.open('{}/barcodes.tsv.gz'.format(options.output_dir), 'wt')

    for line in f:
        if bc_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        print(line, file=bc_out)
        bc_idx += 1
    
    bc_out.close()
    f.close()
    
    f = None
    features_gz = False
    if options.features[-3:] == '.gz':
        features_gz = True
        f = gzip.open(options.features, 'r')
    else:
        f = open(options.features, 'r')
    
    features_out = gzip.open('{}/features.tsv.gz'.format(options.output_dir), 'wt')
    
    feature_idx_new = 1
    feature_idx = 1
    feature_idx2idx = {}

    for line in f:
        if features_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        
        dat = line.split('\t')
        if len(dat) < 3:
            print("ERROR: no feature type; only {} fields".format(len(dat)), file=sys.stderr)
            exit(1)
        feature = dat[2]
        if feature == options.type:
            print(line, file=features_out)
            feature_idx2idx[feature_idx] = feature_idx_new
            feature_idx_new += 1
        
        feature_idx += 1

    f.close()

    mtx_out = gzip.open('{}/matrix.mtx.gz'.format(options.output_dir), 'wt')

    # Print MTX header to each file
    print("%%MatrixMarket matrix coordinate integer general", file=mtx_out)
    print("%", file=mtx_out)

    mtx_out_lines = []

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
                    if j in feature_idx2idx:
                        mtx_out_lines.append("{} {} {}".format(i, feature_idx2idx[j], count))    
                elif i in feature_idx2idx:
                    mtx_out_lines.append("{} {} {}".format(feature_idx2idx[i], j, count))
    
    f.close()

    if bc_col1:
        print("{} {} {}".format(bc_idx-1, feature_idx_new-1, len(mtx_out_lines)), \
            file=mtx_out)
    else:
        print("{} {} {}".format(feature_idx_new-1, bc_idx-1, len(mtx_out_lines)), \
            file=mtx_out)

    for line in mtx_out_lines:
        print(line, file=mtx_out)

    mtx_out.close()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
