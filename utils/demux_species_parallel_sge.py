#! /usr/bin/env python3
import sys
import os
import glob
import subprocess
import re
import argparse
"""
Given raw RNA-seq (required), ATAC-seq (optional), and/or feature barcode (optional) data,
runs demux_species in parallel on an SGE cluster.
"""

def parse_libs(filename):
    """
    Parse a file containing a list of (RNA-seq)
    library names, one per line.
    """
    libs = []
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            libs.append(line)
    f.close()
    return libs

def parse_map(filename):
    """
    Parse a file mapping ATAC lib name -> GEX lib name
    If not provided, assume names are identical
    """
    m = {}
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            custom, gex = line.split('\t')
            m[gex] = custom
    f.close()
    return m

def parse_map_custom(filename):
    """
    Parse a lib name map file for custom-type seq data.
    Lines should be lib name -> data type -> (OPTIONALLY) GEX lib name
    tab separated
    """
    m = {}
    names = {}
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            dat = line.split('\t')
            if len(dat) == 2:
                names[dat[0]] = dat[1]
            elif len(dat) > 2:
                custom, datatype, gex = dat
                names[custom] = datatype
                m[gex] = custom
            else:
                print("ERROR parsing custom lib name map file {}".format(filename), file=sys.stderr)
                print("Make sure lines have 2 or 3 fields, tab-separated.", file=sys.stderr)
                exit(1)

    f.close()
    return (m, names)

def parse_args():
    """
    Handle command-line arguments.
    """
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    reads = parser.add_argument_group()
    reads.add_argument("--rna_dir", "-rd", help="Directory containing RNA-seq reads", 
        required=True)
    reads.add_argument("--libs", "-l", help="Text file listing names of libraries", 
        type=parse_libs, required=True)
    reads.add_argument("--atac_dir", "-ad", help="Directory containing ATAC-seq reads (optional)", 
        required=False, default=None)
    reads.add_argument("--atac_map", "-am", help="If ATAC & GEX libraries have different names, \
pass a file mapping ATAC to GEX names, tab separated, one per line.",
        required=False, default=None, type=parse_map)
    reads.add_argument("--custom_dir", "-cd", help="Directory containing custom-type reads (optional)", 
        required=False, default=None)
    reads.add_argument("--custom_map", "-cm", help="Custom libraries must be given with a file \
mapping library name -> data type. If library names do not match library names for GEX data, then \
a third column should be included, listing corresponding GEX library name. Data types should \
match existing types recognized by CellRanger, i.e. \"CRISPR Guide Capture\" or \"Antibody Capture.\" \
One mapping should be given per line, and fields should be tab-separated.",
    required=False, default=None, type=parse_map_custom)

    parser.add_argument("--kmers", "-k", help="Path to root name for unique k-mers", 
        required=True)
    parser.add_argument("--whitelist", "-w", help="Path to 10X RNA-seq whitelist", 
        required=True)
    parser.add_argument("--whitelist_atac", "-W", help="Path to 10X ATAC-seq whitelist (Required \
if using ATAC data)", required=False)
    parser.add_argument("--output_directory", "-o", help="Final output path (REQUIRED)", 
        required=True)
    parser.add_argument("--tmp_directory", "-t", help="Output path for split files (REQUIRED)", 
        required=True)
    
    parser.add_argument("--max_mem", "-m", type=int, help="Maximum memory (in GB) to request for k-mer counting \
jobs. Default = 60. This number should be higher if you use longer k-mers (i.e. 120-150GB for 45-mers).", \
        default=60, required=False)

    parser.add_argument("--num_chunks", "-n", help="How many pieces should reads be split into? \
Default = 100", type=int, required=False, default=100)

    parsed = parser.parse_args()
    
    if parsed.rna_dir == parsed.tmp_directory or parsed.atac_dir == parsed.tmp_directory or parsed.custom_dir == parsed.tmp_directory:
        print("ERROR: please make --tmp_directory different from read file directories", file=sys.stderr)
        exit(1)
    
    if parsed.output_directory == parsed.tmp_directory:
        print("ERROR: please make --tmp_directory different from --output_directory", file=sys.stderr)
        exit(1)

    if parsed.atac_dir is not None and parsed.whitelist_atac is None:
        print("ERROR: ATAC whitelist required if using ATAC data", file=sys.stderr)
        exit(1)
    
    if parsed.atac_map is not None and parsed.atac_dir is None:
        print("ERROR: ATAC mapping given without an ATAC dir", file=sys.stderr)
        exit(1)

    if parsed.custom_map is not None and parsed.custom_dir is None:
        print("ERROR: custom mapping given without a custom dir", file=sys.stderr)
        exit(1)
    
    if parsed.custom_dir is not None and parsed.custom_map is None:
        print("ERROR: custom read data must be given with a custom mapping file, to determine data types",
            file=sys.stderr)
        exit(1)
    return parsed

def get_script_base():
    """
    Return the beginning of a basic shell script, as a list of lines.
    """
    return ['#! /usr/bin/env bash', '#$ -V', '#$ -S /bin/bash', '#$ -cwd']

def get_split_script(tmpdir):
    lines = get_script_base()
    lines.append('#$ -N demux_species__split_read_file')
    lines.append('#$ -l h_rt=24:00:00')
    lines.append('#$ -o {}/{}/split.out'.format(os.getcwd(), tmpdir))
    lines.append('#$ -e {}/{}/split.err'.format(os.getcwd(), tmpdir))
    return '\n'.join(lines)

def get_count_script(libdir, jid, max_mem, num_chunks):
    lines = get_script_base()
    lines.append('#$ -N demux_species__count')
    # Hopefully this is liberal enough: shouldn't need more than 10 hours per chunk
    # after splitting into many pieces
    lines.append('#$ -l h_rt=10:00:00')
    # Ensure it gets a lot of memory
    lines.append('#$ -l mem_free={}G'.format(max_mem))
    lines.append('#$ -o {}/{}/count.out'.format(os.getcwd(), libdir))
    lines.append('#$ -e {}/{}/count.err'.format(os.getcwd(), libdir))
    lines.append('#$ -hold_jid {}'.format(jid))
    lines.append('#$ -t 1-{}'.format(num_chunks))
    return '\n'.join(lines)

def get_join_script(libdir, jids):
    lines = get_script_base()
    lines.append('#$ -N demux_species__combine_species_counts')
    lines.append('#$ -l h_rt=01:00:00')
    lines.append('#$ -o {}/{}/join.out'.format(os.getcwd(), libdir))
    lines.append('#$ -e {}/{}/join.err'.format(os.getcwd(), libdir))
    lines.append('#$ -hold_jid {}'.format(','.join(jids)))
    return '\n'.join(lines)

def get_assn_script(libdir, jid):
    lines = get_script_base()
    lines.append('#$ -N demux_species__assn')
    # This should hopefully not take nearly this long
    lines.append('#$ -l h_rt=02:00:00')
    lines.append('#$ -l mem_free=25G')
    lines.append('#$ -o {}/{}/assn.out'.format(os.getcwd(), libdir))
    lines.append('#$ -e {}/{}/assn.err'.format(os.getcwd(), libdir))
    lines.append('#$ -hold_jid {}'.format(jid))
    return '\n'.join(lines)

def get_demux_script(libdir, jid):
    lines = get_script_base()
    lines.append('#$ -N demux_species__demux')
    lines.append('#$ -l h_rt=48:00:00')
    lines.append('#$ -l mem_free=25G')
    lines.append('#$ -o {}/{}/demux.out'.format(os.getcwd(), libdir))
    lines.append('#$ -e {}/{}/demux.err'.format(os.getcwd(), libdir))
    lines.append('#$ -hold_jid {}'.format(jid))
    return '\n'.join(lines)

def get_plot_script(libdir, jid):
    lines = get_script_base()
    lines.append('#$ -N demux_species__plot')
    lines.append('#$ -l h_rt=00:10:00')
    lines.append('#$ -o {}/{}/plot.out'.format(os.getcwd(), libdir))
    lines.append('#$ -e {}/{}/plot.err'.format(os.getcwd(), libdir))
    lines.append('#$ -hold_jid {}'.format(jid))
    return '\n'.join(lines)

def launch_job(script):
    """
    Given a shell script in text form, submits it as a job using qsub and
    returns the resulting job ID
    """
    #print(script, file=sys.stderr)
    #print("", file=sys.stderr)
    
    p = subprocess.Popen(['qsub', '-terse'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate(input=script.encode())

    # Get job ID
    jid = out.decode().strip('\n').split('.')[0]
    print("Queued job {}".format(jid), file=sys.stderr)
    return jid

def main(args):
    """
    Main method. Launches a bunch of SGE jobs.
    """
    options = parse_args()
    
    if options.rna_dir[-1] == '/':
        options.rna_dir = options.rna_dir[0:-1]
    if options.atac_dir is not None and options.atac_dir[-1] == '/':
        options.atac_dir = options.atac_dir[0:-1]
    if options.custom_dir is not None and options.custom_dir[-1] == '/':
        options.custom_dir = options.custom_dir[0:-1]
    
    # Make sure stuff exists
    # Directories of reads
    if not os.path.isdir(options.rna_dir):
        print("ERROR: directory {} not found".format(options.rna_dir), file=sys.stderr)
        exit(1)
    if options.atac_dir is not None and not os.path.isdir(options.atac_dir):
        print("ERROR: directory {} not found".format(options.atac_dir), file=sys.stderr)
        exit(1)
    if options.custom_dir is not None and not os.path.isdir(options.custom_dir):
        print("ERROR: directory {} not found".format(options.custom_dir), file=sys.stderr)
        exit(1)
    
    # Unique k-mer files
    if not os.path.isfile("{}.0.kmers".format(options.kmers)) or \
        not os.path.isfile("{}.names".format(options.kmers)):
        print("ERROR: kmers not found at {}".format(options.kmers), file=sys.stderr)
        exit(1)
    
    # CellBouncer programs
    cbpath = '/'.join(__file__.split('/')[:-2])
    if not os.path.isfile('{}/demux_species'.format(cbpath)) or \
        not os.path.isfile('{}/utils/split_read_files'.format(cbpath)):
        print("ERROR: cellbouncer not found at {}".format(cbpath), file=sys.stderr)
        exit(1)
    
    # RNA-seq barcode whitelist
    if not os.path.isfile(options.whitelist):
        print("ERROR: whitelist not found at {}".format(options.whitelist), file=sys.stderr)
        exit(1)
    
    # ATAC-seq barcode whitelist
    if options.whitelist_atac is not None and not os.path.isfile(options.whitelist_atac):
        print("ERROR: ATAC whitelist not found at {}".format(options.whitelist_atac), file=sys.stderr)
        exit(1)

    # Create output directories if needed
    if not os.path.isdir(options.output_directory):
        os.mkdir(options.output_directory)
    if not os.path.isdir(options.tmp_directory):
        os.mkdir(options.tmp_directory)
    
    for lib in options.libs:
        
        libdir = '{}/{}'.format(options.tmp_directory, lib)
        if not os.path.isdir(libdir):
            os.mkdir(libdir)
        libdir_out = '{}/{}'.format(options.output_directory, lib)
        if not os.path.isdir(libdir_out):
            os.mkdir(libdir_out)

        lib_idx = 0

        lib_jids = []

        # Get all reads & split
        gex_r1 = []
        gex_r2 = []
        
        atac_r1 = []
        atac_r2 = []
        atac_r3 = []

        custom_r1 = []
        custom_r2 = []
        custom_names = []

        for fn in glob.glob('{}/{}*_S*_R1_001.fastq.gz'.format(options.rna_dir, lib)):
            fn1 = fn
            fn2 = fn1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
            
            gex_r1.append('-r {}'.format(fn1))
            gex_r2.append('-R {}'.format(fn2))

            script = get_split_script(options.tmp_directory) + '\n' + \
                "{}/utils/split_read_files -1 {} -2 {} -n {} -o {}".format(cbpath,
                fn1, fn2, options.num_chunks, libdir)
            
            jid = launch_job(script)

            # Create a job to run demux_species on this chunk
            script2 = get_count_script(libdir, jid, options.max_mem, options.num_chunks) + '\n' + \
                '{}/demux_species -r {}/{}.${{SGE_TASK_ID}}.fastq.gz -R {}/{}.${{SGE_TASK_ID}}.fastq.gz -w {} -k {} -b $(( $SGE_TASK_ID + {} )) -o {}'.format(\
                cbpath, libdir, fn1.split('/')[-1].split('.fastq.gz')[0], \
                libdir, fn2.split('/')[-1].split('.fastq.gz')[0], options.whitelist, options.kmers, \
                lib_idx * options.num_chunks, libdir_out)
            
            jid = launch_job(script2)
            lib_jids.append(jid)
            
            lib_idx += 1
        
        # Search for ATAC data, if it exists
        if options.atac_dir is not None:
            lib_atac = lib
            if options.atac_map is not None and lib in options.atac_map:
                lib_atac = options.atac_map[lib]
            for fn in glob.glob('{}/{}*_S*_R1_001.fastq.gz'.format(options.atac_dir, lib_atac)):
                fn1 = fn
                fn2 = fn1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                fn3 = fn1.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")
                atac_r1.append('-1 {}'.format(fn1))
                atac_r2.append('-2 {}'.format(fn2))
                atac_r3.append('-3 {}'.format(fn3))
        
        # Search for custom read type data, if it exists
        if options.custom_dir is not None:
            lib_custom = lib
            if options.custom_map[0] is not None and lib in options.custom_map[0]:
                lib_custom = options.custom_map[0][lib]
            lib_datatype = options.custom_map[1][lib_custom]
            for fn in glob.glob('{}/{}*_S*_R1_001.fastq.gz'.format(options.custom_dir, lib_custom)):
                fn1 = fn
                fn2 = fn1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                custom_r1.append('-x {}'.format(fn1))
                custom_r2.append('-X {}'.format(fn2))
                custom_names.append('-N {}'.format(lib_datatype))

        # Queue up script to combine counts
        script = get_join_script(libdir, lib_jids) + '\n' + \
            '{}/utils/combine_species_counts -o {}/{}'.format(cbpath, os.getcwd(), libdir_out)
        
        jid = launch_job(script)

        # Queue up job to assign species & create library files
        
        script_extra = []
        if len(gex_r1) > 0:
            script_extra.append(' '.join(gex_r1))
            script_extra.append(' '.join(gex_r2))
        if len(atac_r1) > 0:
            script_extra.append(' '.join(atac_r1))
            script_extra.append(' '.join(atac_r2))
            script_extra.append(' '.join(atac_r3))
        if len(custom_r1) > 0:
            script_extra.append(' '.join(custom_r1))
            script_extra.append(' '.join(custom_r2))
            script_extra.append(' '.join(custom_names))

        script = get_assn_script(libdir, jid) + '\n' + \
            '{}/demux_species -o {} -d -w {} '.format(cbpath, libdir_out, options.whitelist) + \
            ' ' + ' '.join(script_extra)
        
        jid = launch_job(script)

        # Create plots
        script = get_plot_script(libdir, jid) + '\n' + \
            '{}/plot/species.R {}/{}'.format(cbpath, os.getcwd(), libdir_out)
        
        jid2 = launch_job(script)

        # Queue up final jobs to finish demultiplexing
        for idx in range(0, len(gex_r1)):
            fn1 = gex_r1[idx]
            fn2 = gex_r2[idx]
            script = get_demux_script(libdir, jid) + '\n' + \
                '{}/demux_species -w {} -o {}'.format(cbpath, options.whitelist, libdir_out) + \
                ' ' + fn1 + ' ' + fn2
            jid2 = launch_job(script)

        for idx in range(0, len(atac_r1)):
            fn1 = atac_r1[idx]
            fn2 = atac_r2[idx]
            fn3 = atac_r3[idx]
            script = get_demux_script(libdir, jid) + '\n' + \
                '{}/demux_species -w {} -W {} -o {}'.format(cbpath, options.whitelist, \
                options.whitelist_atac, libdir_out) + \
                ' ' + fn1 + ' ' + fn2 + ' ' + fn3
            jid2 = launch_job(script)

        for idx in range(0, len(custom_r1)):
            fn1 = custom_r1[idx]
            fn2 = custom_r2[idx]
            name = custom_names[idx]
            script = get_demux_script(libdir, jid) + '\n' + \
                '{}/demux_species -w {} -o {}'.format(cbpath, options.whitelist, libdir_out) + \
                ' ' + fn1 + ' ' + fn2 + ' ' + name
            jid2 = launch_job(script)
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
