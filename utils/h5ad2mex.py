#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
import anndata as ad
import scipy.io
import pandas as pd
"""
Given an h5ad format output file containing counts of tags or
sgRNAas per cell, converts it into MEX format, which can be
read by demux_tags.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", "-H", help="The h5ad file to load", \
        required=True)
    parser.add_argument("--out", "-o", help="The directory for output files", \
        required=True)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    # Ensure output directory exists
    if not os.path.isdir(options.out):
        os.mkdir(options.out)
    
    print("Loading {}...".format(options.h5), file=sys.stderr)
    adata = ad.read_h5ad(options.h5)
    
    print("Writing barcodes...", file=sys.stderr)
    f = gzip.open('{}/barcodes.tsv.gz'.format(options.out), 'wt')
    for bc in adata.obs.index:
        print(bc, file=f)
    f.close()

    print("Writing features...", file=sys.stderr)
    f = gzip.open('{}/features.tsv.gz'.format(options.out), 'wt')
    for gene in adata.var.index:
        print(gene, file=f)
    f.close()

    print("Writing mtx...", file=sys.stderr)
    f = gzip.open('{}/matrix.mtx.gz'.format(options.out), 'wb')
    scipy.io.mmwrite(f, adata.X, field='integer')
    f.close()



if __name__ == '__main__':
    sys.exit(main(sys.argv))
