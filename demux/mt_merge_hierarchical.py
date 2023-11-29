#! /usr/bin/env python3
import sys
import os
import subprocess
import argparse
"""
After running demux_mt hierarchically, this program can combine
information from the separate runs to produce one set of haplotypes
and assignments.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", "-r", help="Output file prefix for first/top-level \
run", required=True)
    parser.add_argument("--leaves", "-l", help="Output file prefixes for each \
second/leaf-level run. Must be provided in same order as assigned cluster number \
(i.e. out.0 should come before out.1). If for any top-level cluster, you do not \
wish to sub-cluster, provide \"NA\" in place of the output prefix for that cluster. \
So if you had 3 root-level clusters, and you sub-clustered the first and last \
but not the second, you could pass --leaves out.0 NA out.1.", nargs="+", required=True)
    parser.add_argument("--bam", "-b", help="The BAM file containing alignment data \
for all clusters", required=True)
    parser.add_argument("--out", "-o", help="The output file prefix for the \
combined/merged data.", required=True)
    parser.add_argument("--doublet_rate", "-D", type=float, required=False, \
        default=0.5, help="Doublet rate prior for combined demux_mt run")
    parser.add_argument("--barcodes", "-B", required=False, default=None, \
        help="Limit final output file to this list of valid barcodes")
    return parser.parse_args()

def parse_vars(filename, vars_combined):
    f = open('{}.vars'.format(filename), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        # Let earlier vars definitions supersede later, so root ref/alt bases take precedence
        key = (dat[0], int(dat[1]))
        if key not in vars_combined:
            vars_combined[key] = dat[2:]
    f.close()

def main(args):
    options = parse_args()
    script_dir = '/'.join(os.path.abspath(__file__).split('/')[0:-1])

    print("Combine sites and map barcodes to identities...", file=sys.stderr)

    # First, combine all vars into one file.
    vars_combined = {}
    parse_vars(options.root, vars_combined)
    for fn in options.leaves:
        if fn != 'NA' and fn != 'na':
            parse_vars(fn, vars_combined)
    
    # Now write out vars to output file
    vars_out = open('{}.vars'.format(options.out), 'w')
    
    # Detect what the mitochondrion is called in the BAM file
    chrM = None
    for chrom, pos in sorted(vars_combined.keys(), key=lambda x: x[1]):
        if chrM is None:
            chrM = chrom
        elif chrom != chrM:
            print("ERROR: multiple names of mitochondrial seqs detected: {} vs {}".\
                format(chrom, chrM), file=sys.stderr)
            exit(1)
        print("{}\t{}\t{}".format(chrom, pos, "\t".join(vars_combined[(chrom, pos)])), \
            file=vars_out)
    vars_out.close()

    # Next, create a "bcmap" file mapping cells to haplotypes.
    bcmap_out = open('{}.bcmap'.format(options.out), 'w')
    
    idx2name = {}
    if os.path.isfile('{}.ids'.format(options.root)):
        f = open('{}.ids'.format(options.root), 'r')
        for idx, line in enumerate(f):
            line = line.rstrip()
            idx2name[idx] = line
        f.close()
    else:
        f = open('{}.haps'.format(options.root), 'r')
        for idx, line in enumerate(f):
            idx2name[idx] = str(idx)
        f.close()

    # Determine which clusters to keep from the root clustering
    keep_root = []
    for leafIdx, fn in enumerate(options.leaves):
        if fn == 'NA' or fn == 'na':
            keep_root.append(idx2name[leafIdx])
    if len(keep_root) > 0:
        keep_root = set(keep_root)
        f = open('{}.assignments'.format(options.root), 'r')
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            if dat[2] == 'S' and dat[1] in keep_root:
                print("{}\t{}".format(dat[0], dat[1]), file=bcmap_out)
        f.close()
    # Add in cluster assignments from all sub-clustered clusters
    for leafIdx, fn in enumerate(options.leaves):
        if fn != 'NA' and fn != 'na':
            f = open('{}.assignments'.format(fn), 'r')
            for line in f:
                line = line.rstrip()
                dat = line.split('\t')
                if dat[2] == 'S':
                    clust_id = '{}_{}'.format(idx2name[leafIdx], dat[1])
                    print("{}\t{}".format(dat[0], clust_id), file=bcmap_out)
            f.close()
    bcmap_out.close()
    
    print("Obtain allelic state for all cells at all sites...", \
        file=sys.stderr)

    # Now, dump haplotypes using the chosen sites by re-running demux_mt with 
    # the dump option.
    subprocess.call(['{}/demux_mt'.format(script_dir), '-b', options.bam, '-m', chrM, \
        '-v', '{}.vars'.format(options.out), '-d', '-o', options.out])
    
    print("Create consensus haplotypes...", file=sys.stderr)

    # Next, create consensus haplotypes from what we just obtained
    subprocess.call(['Rscript', '{}/../utils/mt_consensus_haps.R'.format(script_dir), \
        options.out, '{}.bcmap'.format(options.out), options.out])
    
    print("Create final barcode -> identity assignments...", file=sys.stderr)

    # Finally, re-run demux_mt to get assignments
    cmd = ['{}/demux_mt'.format(script_dir), \
        '-b', options.bam, \
        '-m', chrM, \
        '-v', '{}.vars'.format(options.out), \
        '-H', '{}.haps'.format(options.out), \
        '-i', '{}.ids'.format(options.out), \
        '-o', options.out, \
        '-D', '{}'.format(options.doublet_rate)]
    if options.barcodes is not None:
        cmd.append("-B")
        cmd.append(options.barcodes)

    subprocess.call(cmd)
    
    # Clean up
    os.unlink('{}.bcmap'.format(options.out))

    # Make plots
    subprocess.call(['Rscript', '{}/../plot/demux_mt_clust.R'.format(script_dir), \
        '{}'.format(options.out)])
    subprocess.call(['Rscript', '{}/../plot/demux_mt_unclust.R'.format(script_dir), \
        '{}'.format(options.out)])

if __name__ == '__main__':
    sys.exit(main(sys.argv))
