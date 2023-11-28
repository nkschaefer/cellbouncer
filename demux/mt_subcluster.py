#! /usr/bin/env python3
import sys
import os
import subprocess
import argparse
from collections import defaultdict
"""
If a data set contains very divergent groups of individuals (i.e. multiple species),
running demux_mt once may only capture the deepest divergences among haplotypes. 

If this is the case, you can run this program on the output of a demux_mt run
and tell it which clusters to subcluster. It will run demux_mt on each desired
subcluster, and then merge the output hierarchically into one output data set.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", "-b", required=True, 
        help="The BAM file containing alignment data for all clusters")
    parser.add_argument("--inbase", "-i", required=True,
        help="The demux_mt output file prefix for the run to be sub-clustered")
    parser.add_argument("--out", "-o", 
        help="The output file prefix for the final (subclustered, merged) data set.", required=True)
    parser.add_argument("--clusts", "-c", nargs="+", required=True, \
        help="The cluster indices (or names) to subcluster")
    parser.add_argument("--doublet_rate", "-D", type=float, required=False, default=0.5, \
        help="The doublet rate to use in the final, merged re-clustering.")
    return parser.parse_args()

def main(args):
    options = parse_args()
    script_dir = '/'.join(os.path.abspath(__file__).split('/')[0:-1])
    
    # Get chromosome name
    chrM = None
    f = open('{}.vars'.format(options.inbase), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        chrM = dat[0]
        break

    # Note: here, we allow doublets and ignore doublet calls, to make 
    # sure we get the cleanest possible haplotypes.
    # The user will need to decide whether to allow doublet calls in 
    # the final run, using the merge_hierarchical.py
    # program.

    print("Split barcodes by cluster...", file=sys.stderr)
    clust2bc = defaultdict(list)
    f = open('{}.assignments'.format(options.inbase), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if dat[2] == 'S':
            clust2bc[dat[1]].append(dat[0])
    f.close()
    
    for clust in options.clusts:
        if clust not in clust2bc:
            print("ERROR: cluster {} not found in assignments {}.assignments".format(\
                clust, options.out), file=sys.stderr)
            exit(1)
    
    clustset = set(options.clusts)
    merge_str = []

    # Write out barcode files
    for clust in clust2bc:
        if clust in clustset:
            f = open('{}.{}.bcs'.format(options.out, clust), 'w')
            for bc in clust2bc[clust]:
                print(bc, file=f)
            f.close()
            # Replace root cluster with subclusters for this cluster
            merge_str.append('{}.{}'.format(options.out, clust))
        else:
            # Keep root when merging hierarchically
            merge_str.append('NA')

    fnroot_all = []
    for clust in sorted(list(clustset)):
        fnroot_all.append('{}.{}'.format(options.out, clust))
        print("Subcluster root cluster {}...".format(clust), file=sys.stderr)
        cmd = ['{}/demux_mt'.format(script_dir), '-b', options.bam, '-m', \
            chrM, '-f', '{}.{}.bcs'.format(options.out, clust), \
            '-o', '{}.{}'.format(options.out, clust), '-D', '0.0']
        subprocess.call(cmd)
    
    # Merge hierarchically
    cmd = ['{}/mt_merge_hierarchical.py'.format(script_dir), '-b', options.bam, \
        '-r', options.inbase, '-D', '{}'.format(options.doublet_rate), \
        '-o', options.out]
    cmd.append('-l')
    cmd += merge_str
    subprocess.call(cmd)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
