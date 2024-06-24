#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
import anndata as ad
import scipy.io
import pandas as pd
"""
Given an h5 format (Cellranger or scanpy) output file containing
counts of tags or sgRNAas per cell, converts it into MEX format, 
which can be read by demux_tags.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", "-H", help="The h5ad file to load", \
        required=True)
    parser.add_argument("--out", "-o", help="The directory for output files", \
        required=True)
    parser.add_argument("--scanpy", "-s", help="Specify input is .h5ad format from \
AnnData/scanpy. Default is to assume 10X Genomics format h5", \
        action="store_true", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    # Ensure output directory exists
    if not os.path.isdir(options.out):
        os.mkdir(options.out)
    
    print("Loading {}...".format(options.h5), file=sys.stderr)
    adata = None
    try:
        if not options.scanpy:
            import scanpy as sc
            adata = sc.read_10x_h5(options.h5, gex_only=False)
        else:
            adata = ad.read_h5ad(options.h5)
    except:
        print("ERROR loading input h5 file. Did you specify (or omit) --scanpy correctly?", \
            file=sys.stderr)
        exit(1)

    print("Writing barcodes...", file=sys.stderr)
    f = gzip.open('{}/barcodes.tsv.gz'.format(options.out), 'wt')
    for bc in adata.obs.index:
        print(bc, file=f)
    f.close()

    print("Writing features...", file=sys.stderr)
    f = gzip.open('{}/features.tsv.gz'.format(options.out), 'wt')
    if options.tenX:
        adata.var.to_csv(f, header=False, sep='\t')
    else:
        for gene in adata.var.index:
            print(gene, file=f)
    f.close()

    print("Writing mtx...", file=sys.stderr)
    f = gzip.open('{}/matrix.mtx.gz'.format(options.out), 'wb')
    scipy.io.mmwrite(f, adata.X, field='integer')
    f.close()



if __name__ == '__main__':
    sys.exit(main(sys.argv))
