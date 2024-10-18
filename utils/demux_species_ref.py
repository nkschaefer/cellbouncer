#! /usr/bin/env python3
import sys
import os
import glob
import subprocess
import argparse
import shutil
"""
Builds reference data for demux_species, given transcriptomes from multiple
species. Requires FASTK to be available in $PATH.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", "-k", help="Length of k-mers", type=int,
        default=32)
    parser.add_argument("--num", "-N", help="Number of k-mers to sample. To \
include all, set to -1.", type=int, default=10000000)
    parser.add_argument("--out", "-o", help="Base name for output files. This \
will be what you provide as the -k argument to demux_species.", required=True)
    parser.add_argument("--names", "-n", help="Names of species, space separated. \
--names, --fasta, (and optionally --gtf) must be the same number and in the same \
order.", required=True, nargs="+")
    parser.add_argument("--fasta", "-f", help="Two or more FASTA files for \
species to demultiplex. If you also provide matching --gtf arguments, it will \
attempt to extract transcriptome sequences from these FASTA files. If you do not, \
(e.g. if you already have FASTA files of transcriptomes), it will index these \
FASTA files directly.", required=True, nargs="+")
    parser.add_argument("--gtf", "-g", help="If you want to extract transcripts \
from a reference genome, provide all species' genomes via the --fasta argument \
and provide a GTF for each, in the same order, via this argument.", 
        required=False, nargs="+")
    
    parser.add_argument("--FastK", "-F", help="If you do not have FastK available \
in your $PATH, you can specify the path to it here", required=False)
    parser.add_argument("--gffread", "-G", help="If you do not have gffread available \
in your $PATH, you can specify the path to it here", required=False)

    options = parser.parse_args()
    if len(options.names) < 2:
        print("ERROR: at least two species must be provided.", file=sys.stderr)
        exit(1)
    if len(options.names) != len(options.fasta):
        print("ERROR: you must provide one name and one FASTA per species.", 
            file=sys.stderr)
        exit(1)
    if options.gtf is not None and \
        len(options.gtf) > 0 and len(options.gtf) != len(options.fasta):
        print("ERROR: if providing GTF files, you must provide one for each \
FASTA.", file=sys.stderr)
        exit(1)
    
    if options.FastK is not None and not os.path.isfile(options.FastK):
        print("ERROR: {} does not exist".format(options.FastK), file=sys.stderr)
        exit(1)
    if options.gffread is not None and not os.path.isfile(options.gffread):
        print("ERROR: {} does not exist".format(options.gffread), file=sys.stderr)
        exit(1)

    return options

def is_gz_file(filepath):
    with open(filepath, 'rb') as test:
        return test.read(2) == b'\x1f\x8b'

def extract_tx(fasta, gtf, out, gffread=None):
    """
    Use gffread to extract transcript FASTAs from a reference genome (FASTA) 
    and an annotation (GTF).
    """
    gffread_bin = 'gffread'
    if gffread is not None:
        gffread_bin = gffread
    
    if gffread is None and shutil.which("gffread") is None:
        print("ERROR: gffread is not installed/available in $PATH", file=sys.stderr)
        exit(1)
    
    to_zip = []
    txfiles = []
    for idx in range(0, len(fasta)):
        this_fasta = fasta[idx]
        this_gtf = gzip[idx]
        
        # Note: gffread can't handle gzipped FASTA or GTF.
        # Check for this
        if is_gz_file(fasta[idx]):
            print("FASTA is gzipped. Unzipping...", file=sys.stderr)
            if shutil.which("gunzip") is None:
                print("ERROR: gzip not installed/available in $PATH", file=sys.stderr)
                exit(1)
            else:
                subprocess.call(['gunzip', this_fasta])
                fn_new = this_fasta.split('.gz')[0]
                if os.path.isfile(fn_new):
                    to_zip.append(fn_new)
                    this_fasta = fn_new
                elif os.path.isfile(this_fasta):
                    to_zip.append(this_fasta)
                else:
                    print("ERROR: can't find unzipped version of {}".format(this_fasta), 
                        file=sys.stderr)
                    exit(1)

        if is_gz_file(gtf[idx]):
            print("Annotation is gzipped. Unzipping...", file=sys.stderr)
            if shutil.which("gunzip") is None:
                print("ERROR: gzip is not installed/available in $PATH", file=sys.stderr)
                exit(1)
            else:
                subprocess.call(['gunzip', this_gtf])
                fn_new = this_gtf.split('.gz')[0]
                if os.path.isfile(fn_new):
                    to_zip.append(fn_new)
                    this_gtf = fn_new
                elif os.path.isfile(this_gtf):
                    to_zip.append(this_gtf)
                else:
                    print("ERROR: can't find unzipped version of {}".format(this_gtf),
                        file=sys.stderr)
                    exit(1)
        
        # Extract transcripts
        print("Extracting transcripts from {} and {}...".format(this_fasta, this_gtf), 
            file=sys.stderr)
        
        subprocess.call([gffread_bin, '-F', '-w', '{}.tx.{}.fasta'.format(out, idx), \
            '-g', this_fasta, this_gtf])
        txfiles.append('{}.tx.{}.fasta'.format(out, idx))

    return (to_zip, txfiles)

def count_kmers(fasta, out, k, fastk=None):
    fastk_bin = 'FastK'
    if fastk is not None:
        fastk_bin = fastk

    if fastk is None and shutil.which('FastK') is None:
        print("ERROR: FastK is not installed or unavailable in $PATH", file=sys.stderr)
        exit(1)

    print("Counting k-mers in {}...".format(fasta), file=sys.stderr)
    subprocess.call([fastk_bin, '-N{}'.format(out), '-k{}'.format(k), '-t1', fasta])

def re_gzip(files):
    if len(files) > 0:
        print("Re-zipping decompressed files...", file=sys.stderr)
    for fn in files:
        if shutil.which('bgzip') is not None:
            subprocess.call(['bgzip', fn])
        elif shutil.which('gzip') is not None:
            subprocess.call(['gzip', fn])
        else:
            print("ERROR: gzip/bgzip not installed or unavailable in $PATH", 
                file=sys.stderr)
            return

def get_unique_kmers(names, ktabs, out, num):
    
    # Get directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    if script_dir[-1] == '/':
        script_dir = script_dir[0:-1]

    cmd = ['{}/get_unique_kmers'.format(script_dir), '-N', str(num), '-o', out]
    for idx in range(0, len(names)):
        cmd.append('-n')
        cmd.append(names[idx])
        cmd.append('-k')
        cmd.append(ktabs[idx])

    print("Getting unique kmers...", file=sys.stderr)
    subprocess.call(cmd)


def main(args):
    options = parse_args()
    
    files_cleanup = []
    to_zip = []
    tx_fa = []

    if options.gtf is not None and len(options.gtf) > 0:
        # Extract transcripts from FASTA
        to_zip, tx_fa = extract_tx(options.fasta, options.gtf, options.gffread)    
        files_cleanup = tx_fa
    else:
        tx_fa = options.fasta
    
    # Count k-mers
    ktabs = []
    for idx, tx in enumerate(tx_fa):
        count_kmers(tx, '{}.fastk.{}'.format(options.out, idx), options.k, options.FastK)
        ktabs.append('{}.fastk.{}.ktab'.format(options.out, idx))
    
    # Re-gzip any unzipped files
    re_gzip(to_zip)

    # Create unique k-mer lists
    get_unique_kmers(options.names, ktabs, options.out, options.num)

    # Clean up
    for fn in files_cleanup:
        os.unlink(fn)
    
    # Remove FastK temporary files / data tables
    for idx in range(0, len(options.names)):
        for fn in glob.glob('{}.fastk.{}.*'.format(options.out, idx)):
            os.unlink(fn)
        for fn in glob.glob('.{}.fastk.{}.*'.format(options.out, idx)):
            os.unlink(fn)



if __name__ == '__main__':
    sys.exit(main(sys.argv))
