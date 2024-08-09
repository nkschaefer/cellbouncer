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
    parser.add_argument("--keep_intermediate", "-k", action="store_true", \
        help="Do not remove temporary files from subcluster runs", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    script_dir = '/'.join(os.path.abspath(__file__).split('/')[0:-1])
    
    if not os.path.isfile(options.bam):
        print("ERROR: BAM file {} does not exist".format(options.bam), file=sys.stderr)
        exit(1)

    # Get chromosome name
    chrM = None
    if not os.path.isfile('{}.vars'.format(options.inbase)):
        print("ERROR: .vars file for {} does not exist".format(options.inbase), file=sys.stderr)
        exit(1)

    f = open('{}.vars'.format(options.inbase), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        chrM = dat[0]
        break
    
    # Translate haplotype names to numeric IDs, so we have everything
    # in the order it appears in the file.
    clustname2idx = {}
    clustname_sorted = []
    if os.path.isfile('{}.ids'.format(options.inbase)):
        f = open('{}.ids'.format(options.inbase))
        for idx, line in enumerate(f):
            line = line.rstrip()
            if line != "":
                clustname2idx[line] = idx
                clustname_sorted.append(line)
        f.close()
    else:
        if not os.path.isfile('{}.haps'.format(options.inbase)):
            print("ERROR: haps file does not exist for {}".format(options.inbase), file=sys.stderr)
            exit(1)
        f = open('{}.haps'.format(options.inbase))
        for idx, line in enumerate(f):
            line = line.rstrip()
            if line != "":
                clustname2idx['{}'.format(idx)] = idx
                clustname_sorted.append('{}'.format(idx))
        f.close()
    
    for clust in options.clusts:
        if clust not in clustname2idx:
            print("ERROR: cluster {} not found in cluster IDs.".format(clust), file=sys.stderr)
            exit(1)

    # Note: here, we allow doublets and ignore doublet calls, to make 
    # sure we get the cleanest possible haplotypes.
    # The user will need to decide whether to allow doublet calls in 
    # the final run, using the merge_hierarchical.py
    # program.

    print("Split barcodes by cluster...", file=sys.stderr)
    clust2bc = defaultdict(list)
    if not os.path.isfile('{}.assignments'.format(options.inbase)):
        print("ERROR: assignments file does not exist for {}".format(options.inbase), file=sys.stderr)
        exit(1)
    f = open('{}.assignments'.format(options.inbase), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if dat[2] == 'S':
            clust2bc[dat[1]].append(dat[0])
    f.close()
    
    for clust in options.clusts:
        if clust not in clust2bc:
            print("ERROR: cluster {} not found in assignments file {}.assignments".format(\
                clust, options.inbase), file=sys.stderr)
            exit(1)
    
    clustset = set(options.clusts)
    merge_str = []

    # Write out barcode files
    for clust in clustname_sorted:
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
    
    outs_cleanup = []
    for ci, clust in enumerate(clustname_sorted):
        if clust in clustset:
            print("Subcluster root cluster {}...".format(clust), file=sys.stderr)
            cmd = ['{}/demux_mt'.format(script_dir), '-b', options.bam, '-m', \
                chrM, '-f', '{}.{}.bcs'.format(options.out, clust), \
                '-o', '{}.{}'.format(options.out, clust), '-D', '0.0']
            subprocess.call(cmd)
            outs_cleanup.append('{}.{}'.format(options.out, clust))
            if not os.path.isfile('{}.{}.haps'.format(options.out, clust)):
                print("  cluster {} had no subclusters".format(clust), file=sys.stderr)
                # Clustering failed.
                # Swap out for root-level assignments.
                merge_str[ci] = 'NA'

    # Merge hierarchically
    cmd = ['{}/utils/mt_merge_hierarchical.py'.format(script_dir), '-b', options.bam, \
        '-r', options.inbase, '-D', '{}'.format(options.doublet_rate), \
        '-o', options.out]
    cmd.append('-l')
    cmd += merge_str
    subprocess.call(cmd)
    
    # Clean up
    if not options.keep_intermediate:
        for outbase in outs_cleanup:
            if os.path.isfile('{}.assignments'.format(outbase)):
                os.unlink('{}.assignments'.format(outbase))
            if os.path.isfile('{}.bcs'.format(outbase)):
                os.unlink('{}.bcs'.format(outbase))
            if os.path.isfile('{}.cellhaps'.format(outbase)):
                os.unlink('{}.cellhaps'.format(outbase))
            if os.path.isfile('{}.haps'.format(outbase)):
                os.unlink('{}.haps'.format(outbase))
            if os.path.isfile('{}.summary'.format(outbase)):
                os.unlink('{}.summary'.format(outbase))
            if os.path.isfile('{}.vars'.format(outbase)):
                os.unlink('{}.vars'.format(outbase))

if __name__ == '__main__':
    sys.exit(main(sys.argv))
