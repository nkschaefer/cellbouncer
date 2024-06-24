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
        required=False)
    parser.add_argument("--feature_type", "-t", help="If not --scanpy/-s, and the \
file contains more than one feature type, subset to a specific feature type here. \
You should wrap it in double quotes (\"\"). This is useful, for example, if you \
processed gene expression data together with sgRNA capture data using CellRanger. \
In this case, you can extract the sgRNA data using -t \"CRISPR Guide Capture\".",
        required=False, default=None)
    parser.add_argument("--list_features", "-l", help="For help choosing a feature \
for the --feature_type/-t argument, choose this argument to print a list of all \
feature types in the input data file and exit", action="store_true", required=False)
    parser.add_argument("--scanpy", "-s", help="Specify input is .h5ad format from \
AnnData/scanpy. Default is to assume 10X Genomics format h5", \
        action="store_true", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    if options.feature_type is not None:
        if options.list_features:
            print("ERROR: cannot specify both --list_features/-l and --feature_type/-t", \
                file=sys.stderr)
            exit(1)
        if options.scanpy:
            print("ERROR: cannot specify a --feature_type/-t with the --scanpy/-s option.", \
                file=sys.stderr)
            exit(1)
    if options.list_features and options.scanpy:
        print("ERROR: cannot specify --list_features/-l with option --scanpy/-s.", file=sys.stderr)
        exit(1)
    if options.out is None and not options.list_features:
        print("ERROR: you must provide an --out/-o directory unless listing features (-l).", \
            file=sys.stderr)
        exit(1)

    # Ensure output directory exists
    if options.out is not None and not os.path.isdir(options.out):
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
    adata.var_names_make_unique()
    
    if options.list_features or options.feature_type is not None:
        if 'feature_types' in adata.var.columns:
            if options.list_features:
                for elt in adata.var['feature_types'].unique():
                    print(elt)
                # Done here
                exit(0)
            else:
                # Subset anndata instead of listing features.
                adata = adata[:,list(adata.var.loc[adata.var['feature_types'] == options.feature_type,:].index)] 
        else:
            print("ERROR: no feature_types column in adata.var. Is this a --scanpy .h5ad file?", \
                file=sys.stderr)
            exit(1)


    print("Writing barcodes...", file=sys.stderr)
    f = gzip.open('{}/barcodes.tsv.gz'.format(options.out), 'wt')
    for bc in adata.obs.index:
        print(bc, file=f)
    f.close()

    print("Writing features...", file=sys.stderr)
    f = gzip.open('{}/features.tsv.gz'.format(options.out), 'wt')
    if not options.scanpy:
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
