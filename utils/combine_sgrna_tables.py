#! /usr/bin/env python3
import sys
import os
"""
If you have run demux_tags on multiple sgRNA libraries and want to load all 
.table files together (i.e. the libraries are part of one big high-MOI CRISPR
screen), run this program to combine them.
"""
def main(args):
    if len(args) < 3:
        print("USAGE: combine_sgrna_tables.py [output_prefix1].table [output_prefix2].table...", \
            file=sys.stderr)
        print("This program assumes you have already kept barcode names unique by specifying", file=sys.stderr)
        print("  a unique --libname\\-n for each run of demux_tags.", file=sys.stderr)
        print("The output file will be printed to stdout.", file=sys.stderr)
        exit(1)
    

    # Load headers.
    names_all = set([])
    
    idx2name = []

    for fn in args[1:]:
        f = open(fn, 'r')
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            idx2name_row = []
            for idx, elt in enumerate(dat[1:]):
                names_all.add(elt)
                idx2name_row.append(elt)
            idx2name.append(idx2name_row)
            break
        f.close()

    name2idx = {}
    idx = 0
    
    print_hdr = ['barcode']
    for name in sorted(list(names_all)):
        name2idx[name] = idx
        idx += 1
        print_hdr.append(name)
    print("\t".join(print_hdr))

    for file_idx, fn in enumerate(args[1:]):
        f = open(fn, 'r')
        first = True
        for line in f:
            if first:
                first = False
            else:
                line = line.rstrip()
                dat = line.split('\t')
                println = [dat[0]] + ['0'] * len(names_all)
                for idx in range(0, len(dat)-1):
                    if int(dat[idx + 1]) != 0:
                        name = idx2name[file_idx][idx]
                        println[name2idx[name] + 1] = dat[idx+1]
                print("\t".join(println))


    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
