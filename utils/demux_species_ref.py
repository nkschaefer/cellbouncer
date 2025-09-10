#! /usr/bin/env python3
import sys
import os
import glob
import subprocess
import argparse
import shutil
from collections import deque, OrderedDict, defaultdict, Counter
import time
"""
Builds reference data for demux_species, given transcriptomes from multiple
species. Requires FASTK to be available in $PATH.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k", "-k", help="Length of k-mers", type=int,
        default=32)
    parser.add_argument("--num", "-N", help="Number of k-mers to sample. To \
include all, set to -1 (which is default). Sampling fewer (i.e. 10 or 20 million) can save \
memory and time at the cost of missing some cells.", type=int, default=-1)
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
    parser.add_argument("--no_orthology", "-O", help="Default behavior: filter \
transcriptomes for sets of homologous transcripts. This is done by aligning each \
pair of transcriptomes using minimap2, clustering whole transcriptomes using UPGMA, \
and then iteratively finding homologous transcript sets at each node in the tree, \
progressing to the root. Transcripts not belonging to a homologous set are discarded \
before counting k-mers. Setting this option will disable this procedure.", action="store_true")

    parser.add_argument("--mm2", "-M", help="If filtering for orthologous transcript \
and you do not have minimap2 in your $PATH, you can specify the path \
to it here", required=False)
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
    
    txfiles = []

    for idx in range(0, len(fasta)):
        this_fasta = fasta[idx]
        this_gtf = gtf[idx]
        
        to_rm = []

        # Note: gffread can't handle gzipped FASTA or GTF.
        # Check for this
        if is_gz_file(fasta[idx]):
            print("FASTA is gzipped. Unzipping...", file=sys.stderr)
            if shutil.which("gunzip") is None:
                print("ERROR: gzip not installed/available in $PATH", file=sys.stderr)
                exit(1)
            else:
                fa_base = this_fasta.split('/')[-1]
                fa_base = fa_base.split('.gz')[0]
                fa_tmp = '{}.{}'.format(out, fa_base)
                fa_f = open(fa_tmp, 'w')
                p = subprocess.Popen(['gunzip', '-c', this_fasta], stdout=fa_f)
                stdout, err = p.communicate()
                fa_f.close()
                to_rm.append(fa_tmp)
                this_fasta = fa_tmp

        if is_gz_file(gtf[idx]):
            print("Annotation is gzipped. Unzipping...", file=sys.stderr)
            if shutil.which("gunzip") is None:
                print("ERROR: gzip is not installed/available in $PATH", file=sys.stderr)
                exit(1)
            else:
                gtf_base = this_gtf.split('/')[-1]
                gtf_base = gtf_base.split('.gz')[0]
                gtf_tmp = '{}.{}'.format(out, gtf_base)
                gtf_f = open(gtf_tmp, 'w')
                p = subprocess.Popen(['gunzip', '-c', this_gtf], stdout=gtf_f)
                stdout, err = p.communicate()
                gtf_f.close()
                to_rm.append(gtf_tmp)
                this_gtf = gtf_tmp

        # Extract transcripts
        print("Extracting transcripts from {} and {}...".format(this_fasta, this_gtf), 
            file=sys.stderr)
        
        subprocess.call([gffread_bin, '-F', '-w', '{}.tx.{}.fasta'.format(out, idx), \
            '-g', this_fasta, this_gtf])
        txfiles.append('{}.tx.{}.fasta'.format(out, idx))
        
        for fn in to_rm:
            os.unlink(fn)

    return txfiles

def count_kmers(fasta, out, k, fastk=None):
    fastk_bin = 'FastK'
    if fastk is not None:
        fastk_bin = fastk

    if fastk is None and shutil.which('FastK') is None:
        print("ERROR: FastK is not installed or unavailable in $PATH", file=sys.stderr)
        exit(1)

    print("Counting k-mers in {}...".format(fasta), file=sys.stderr)
    subprocess.call([fastk_bin, '-N{}'.format(out), '-k{}'.format(k), '-t1', fasta])

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

def get_clique_aux(targets, i, graphs, clique):
    for source in targets[i]:
        targets2 = []
        for xi in range(0, len(graphs)):
            targets2.append(OrderedDict())
        for j in range(i+1, len(graphs)):
            for dest in targets[j]:
                if source in graphs[i] and (j, dest) in graphs[i][source] and\
                    dest in graphs[j] and (i, source) in graphs[j][dest]:
                    # Keep
                    targets2[j][dest] = graphs[i][source][(j, dest)]
        
        targets2_pass = True
        for j in range(i+1, len(targets2)):
            if len(targets2[j]) == 0:
                targets2_pass = False
                break

        if targets2_pass:
            clique.append(source)
            if i < len(graphs)-1:
                cl2 = clique[:]
                ret = get_clique_aux(targets2, i+1, graphs, cl2)
                if ret is not None:
                    return ret
            else:
                # Base case
                return clique
        else:
            return None

def get_clique(graphs, tx_chosen):
    # Algorithm: propose all nodes that have connections to all graphs
    # Require connections between all nodes
    
    clique = None

    for source in graphs[0]:
        if (0, source) not in tx_chosen[0]:
            targets = []
            for x in range(0, len(graphs)):
                targets.append(OrderedDict())
            for dest in graphs[0][source]:
                if dest not in tx_chosen[dest[0]] and dest[1] in graphs[dest[0]] and\
                    (0, source) in graphs[dest[0]][dest[1]]:
                    targets[dest[0]][dest[1]] = graphs[0][source][dest]

            targets_pass = True
            for i in range(1, len(targets)):
                if len(targets[i]) == 0:
                    targets_pass = False
                    break
            if targets_pass:
                # source is a plausible clique member. Check recursively.
                clique_test = [source]
                ret = get_clique_aux(targets, 1, graphs, clique_test)
                if ret is not None:
                    return ret
                    #exit(1)

            if clique is not None:
                break
        if clique is not None:
            break

    
    return clique

def elim_node_recursive(key, graphs, graphs_filt, edges_rev_all):
    idx, cl = key
    if cl in graphs[idx]:
        del graphs[idx][cl]
        if cl in graphs_filt[idx]:
            del graphs_filt[idx][cl]
    
    elim_extra = []

    if key in edges_rev_all:
        for k2 in edges_rev_all[key]:
            if k2[1] in graphs[k2[0]] and key in graphs[k2[0]][k2[1]]:
                del graphs[k2[0]][k2[1]][key]
                filt_del = False
                if k2[1] in graphs_filt[k2[0]] and key in graphs_filt[k2[0]][k2[1]]:
                    del graphs_filt[k2[0]][k2[1]][key]
                    filt_del = True
                if len(graphs[k2[0]][k2[1]]) == 0:
                    del graphs[k2[0]][k2[1]]
                    if filt_del and len(graphs_filt[k2[0]][k2[1]]) == 0:
                        del graphs_filt[k2[0]][k2[1]]
                else:
                    # If a node is unconnected to at least one species, it can be
                    # removed
                    species = set([])
                    for edge in graphs[k2[0]][k2[1]]:
                        species.add(edge[0])
                    if len(species) < len(graphs):
                        # Eliminate it
                        elim_extra.append(k2)

    for k in elim_extra:
        elim_node_recursive(k, graphs, graphs_filt, edges_rev_all)

def proc_graphs(graphs, out):

    clique_outf = open('{}.txgroups'.format(out), 'w')

    edges_rev_all = defaultdict(list)

    # First, eliminate edges that cannot be part of a clique.
    for i in range(0, len(graphs)):
        sdels = []
        for source in graphs[i]:
            dels = []
            species_targets = set([])
            for dest in graphs[i][source]:
                if dest[1] not in graphs[dest[0]] or (i, source) not in graphs[dest[0]][dest[1]]:
                    dels.append(dest)
                else:
                    edges_rev_all[dest].append((i, source))
            for d in dels:
                del graphs[i][source][d]
            if len(graphs[i][source]) == 0:
                sdels.append(source)
        for sd in sdels:
            del graphs[i][sd]

    sorts = []
    for i, g in enumerate(graphs):
        sorts.append(deque())
        sort_tmp = []
        for source in g:
            for dest in g[source]:
                sort_tmp.append((source, dest[0], dest[1], g[source][dest]))
        for tup in sorted(sort_tmp, key=lambda x: -x[3]):
            sorts[i].append(tup)

    # Algorithm:
    # Choose max weight edge from each graph to include
    # Look for fully connected cliques
    # Mark these tx as real and remove from graph
    # Use OrderedDict so earlier keys are higher-scoring
    tx_chosen = []
    graphs_filt = []
    for i in range(0, len(graphs)):
        tx_chosen.append(set([]))
        graphs_filt.append(OrderedDict())

    done_count = 0

    cliquecount = 0
    
    while done_count == 0:
        a = time.time()
        for i in range(0, len(graphs)):
            chosen = None
            while chosen is None and len(sorts[i]) > 0:
                chosen = sorts[i].popleft()
                source, desti, dest, score = chosen
                #if source in tx_chosen[i] or dest in tx_chosen[desti]:
                #    # Can't add
                #    chosen = None
                if False:
                    pass
                elif source not in graphs[i]:
                    chosen = None
                elif dest not in graphs[desti]:
                    chosen = None
                else:
                    # Add to filtered graph
                    if source not in graphs_filt[i]:
                        graphs_filt[i][source] = OrderedDict()
                    graphs_filt[i][source][(desti, dest)] = score

            if chosen is None:
                done_count += 1
        
        b = time.time()
        if done_count == 0:
            # Look for fully connected cliques
            x = time.time()
            clique = get_clique(graphs_filt, tx_chosen)
            y = time.time()
            
            #if cliquecount > 15000:
            if False:
                print("Choose next: {:.4f} seconds".format(b-a))
                print("Get clique: {:.4f} seconds".format(y-x))

            # Deal with choosing this clique
            if clique is not None:
                cliquecount += 1
                print("\t".join(clique), file=clique_outf)
                c = time.time()
                # Eliminate all edges involving a chosen node
                for idx, cl in enumerate(clique):
                    tx_chosen[idx].add(cl)
                    elim_node_recursive((idx, cl), graphs, graphs_filt, edges_rev_all)

                d = time.time()
                #if cliquecount > 15000:
                if False:
                    print("Clean up: {:.4f} seconds".format(d-c))

            if cliquecount % 100 == 0:
                print("Found {} homologous transcript sets".format(cliquecount), file=sys.stderr, end="\r")
    
    clique_outf.close()
    print("Found {} homologous transcript sets".format(cliquecount), file=sys.stderr)
    return tx_chosen

def tx_upgma(distmat, graphs, out, edges_rev):
    tx_chosen = []
    tx_groups = []
    tx_set = {}
    for i in range(0, len(graphs)):
        tx_chosen.append(set([]))
    
    n_orig = len(graphs)

    while len(distmat) > 0:
        tx_groups.clear()

        distmat2 = Counter()
        # Find max
        maxdist = 0
        maxkey = None
        for k in distmat:
            if distmat[k] > maxdist:
                maxdist = distmat[k]
                maxkey = k

        k = maxkey

        idx = len(graphs)
        for key in distmat:
            if not ((key[0] == k[0] or key[0] == k[1]) and (key[1] == k[0] or key[1] == k[1])):
                if (key[0] == k[0] or key[0] == k[1]):
                    distmat2[(key[1], idx)] += 0.5*distmat[key]
                elif (key[1] == k[0] or key[1] == k[1]):
                    distmat2[(key[0], idx)] += 0.5*distmat[key]
                else:
                    distmat2[key] = distmat[key]

        distmat.clear()
        distmat = distmat2

        graphs.append({})
        tx_chosen.append(set([]))

        txidx = 1
        # Join these two.
        txpairs = Counter()
        for tx1 in graphs[k[0]]:
            if k[1] in graphs[k[0]][tx1]:
                for tx2 in graphs[k[0]][tx1][k[1]]:
                    if tx2 in graphs[k[1]] and k[0] in graphs[k[1]][tx2] and \
                        tx1 in graphs[k[1]][tx2][k[0]]:
                        txpairs[(tx1, tx2)] += graphs[k[0]][tx1][k[1]][tx2]
                        txpairs[(tx1, tx2)] += graphs[k[1]][tx2][k[0]][tx1]
           
        for tx, score in sorted(txpairs.items(), key=lambda x: -x[1]):
            if tx[0] not in tx_chosen[k[0]] and tx[1] not in tx_chosen[k[1]]:
                tx_chosen[k[0]].add(tx[0])
                tx_chosen[k[1]].add(tx[1])
                
                # Choose this pair.
                txname = 'TX.{}.{}'.format(idx, txidx)
                txidx += 1
                tx_set[txname] = []
                if tx[0] in tx_set:
                    for elt in tx_set[tx[0]]:
                        tx_set[txname].append(elt)
                else:
                    tx_set[txname].append((k[0], tx[0]))
                
                if tx[1] in tx_set:
                    for elt in tx_set[tx[1]]:
                        tx_set[txname].append(elt)
                else:
                    tx_set[txname].append((k[1], tx[1]))
                
                tx_groups.append(tx_set[txname])

                graphs[idx][txname] = {}
                
                for dest in graphs[k[0]][tx[0]]:
                    if dest != k[1]:
                        for dest_tx in graphs[k[0]][tx[0]][dest]:
                            score1A = graphs[k[0]][tx[0]][dest][dest_tx]
                            if dest_tx in graphs[dest] and k[0] in graphs[dest][dest_tx] and \
                                tx[0] in graphs[dest][dest_tx][k[0]]:
                                score1B = graphs[dest][dest_tx][k[0]][tx[0]]
                                if dest in graphs[k[1]][tx[1]] and dest_tx in graphs[k[1]][tx[1]][dest]:
                                    score2A = graphs[k[1]][tx[1]][dest][dest_tx]
                                    if k[1] in graphs[dest][dest_tx] and tx[1] in graphs[dest][dest_tx][k[1]]:
                                        score2B = graphs[dest][dest_tx][k[1]][tx[1]]
                                        
                                        # Add edges.
                                        scoreA = (score1A + score2A)/2.0
                                        scoreB = (score1B + score2B)/2.0
                                        if dest not in graphs[idx][txname]:
                                            graphs[idx][txname][dest] = {}
                                        if idx not in graphs[dest][dest_tx]:
                                            graphs[dest][dest_tx][idx] = {}
                                        graphs[idx][txname][dest][dest_tx] = scoreA
                                        graphs[dest][dest_tx][idx][txname] = scoreB
                                        
                                        edges_rev[dest].append((idx, txname))
                                        edges_rev[idx].append((dest, dest_tx))

        # Delete old edges from graph.
        graphs[k[0]].clear()
        graphs[k[1]].clear()
        for er in set(edges_rev[k[0]]):
            if er[1] in graphs[er[0]] and k[0] in graphs[er[0]][er[1]]:
                del graphs[er[0]][er[1]][k[0]]
        del edges_rev[k[0]]
        for er in set(edges_rev[k[1]]):
            if er[1] in graphs[er[0]] and k[1] in graphs[er[0]][er[1]]:
                del graphs[er[0]][er[1]][k[1]]
    
    print("Found {} homologous gene sets".format(len(tx_groups)), file=sys.stderr)

    outfname = '{}.txgroups'.format(out)
    outf = open(outfname, 'w')
    for grp in tx_groups:
        grp2 = sorted(list(grp), key=lambda x: x[0])
        print(grp2[0][1], end="", file=outf)
        for i in range(1, len(grp2)):
            print("\t{}".format(grp2[i][1]), end="", file=outf)
        print("\n", end="", file=outf)
    
    outf.close()

    return tx_chosen

def filter_orthology(tx_fa, out, mm2_path=None):
    dbs = [] 
    tx_filt = []
    tx_filt_name = '{}.txgroups'.format(out)
    if os.path.isfile(tx_filt_name):
        for x in range(0, len(tx_fa)):
            tx_filt.append(set([]))

        print("Loading homologous transcript groups from {}...".format(tx_filt_name), file=sys.stderr)
        f = open(tx_filt_name, 'r')
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            for idx, elt in enumerate(dat):
                tx_filt[idx].add(elt)
        f.close()
    else:
        mm2 = 'minimap2'

        if mm2_path is not None:
            if mm2_path[-1] == '/':
                mm2_path = mm2_path[0:-1]
            if 'minimap2' in mm2_path:
                mm2 = mm2_path
            else:
                mm2 = '{}/minimap2'.format(mm2_path)
        
        mm2_args = ['--for-only', '--mask-level', '0.5', '-k15', '-w5', '-p0.5', '-N10', '-x', 'asm20']
        #mm2_args = ['--for-only', '--mask-level', '0.7', '-k12', '-w5', '-p0.5', '-N20', '-x', 'asm20']
        for tx in tx_fa:
            idxname = '{}.{}.idx'.format(out, tx.split('/')[-1])
            if not os.path.isfile(idxname):
                cmd = [mm2] + mm2_args + ['-d', idxname, tx]
                print("Minimap2 idx {}...".format(tx), file=sys.stderr)
                subprocess.call(cmd)
            dbs.append(idxname)
        
        graphs = []
        for i in range(0, len(tx_fa)):
            graphs.append({})

        # Do all pairwise comparisons
        print("Running alignments...", file=sys.stderr)

        # Store mean pairwise distances between transcriptomes
        distmat = Counter()
        distmat_count = Counter()
        
        edges_rev = defaultdict(list)

        for i in range(0, len(tx_fa)-1):
            dbi = dbs[i]
            txi = tx_fa[i]
            for j in range(i+1, len(tx_fa)):
                dbj = dbs[j]
                txj = tx_fa[j]
                cmd1 = [mm2] + mm2_args + [dbi, txj]
                cmd2 = [mm2] + mm2_args + [dbj, txi]
                
                outfn1_this = '{}.{}.{}.paf'.format(out, j, i)
                dbs.append(outfn1_this)
                if not os.path.isfile(outfn1_this):
                    out1_this = open(outfn1_this, 'w')
                    p1 = subprocess.Popen(cmd1, stdout=out1_this)
                    print("  {} vs {}".format(tx_fa[j].split('/')[-1], tx_fa[i].split('/')[-1]), file=sys.stderr)
                    stdout, stderr = p1.communicate()
                    out1_this.close()
                else:
                    print("Loading alignments from {}...".format(outfn1_this), file=sys.stderr)

                in1_this = open(outfn1_this, 'r')
                for line in in1_this:
                    line = line.rstrip()
                    dat = line.split('\t')
                    if len(dat) >= 12:
                        query = dat[0]
                        ref = dat[5]
                        nmatch = int(dat[9])
                        tot = int(dat[10])
                        lquery = int(dat[1])
                        lref = int(dat[6])
                        denom = lquery
                        if lref < lquery:
                            denom = lref
                        score = nmatch/denom
                        distmat[(i, j)] += score
                        distmat_count[(i, j)] += 1

                        if query not in graphs[j]:
                            graphs[j][query] = {}
                        if i not in graphs[j][query]:
                            graphs[j][query][i] = {}
                        graphs[j][query][i][ref] = score
                        edges_rev[i].append((j, query))
                        #key = (i, ref)
                        #if key not in graphs[j][query] or score > graphs[j][query][key]:
                        #    graphs[j][query][key] = score
                in1_this.close()
                outfn2_this = '{}.{}.{}.paf'.format(out, i, j)
                dbs.append(outfn2_this)
                if not os.path.isfile(outfn2_this):
                    out2_this = open(outfn2_this, 'w')
                    print("  {} vs {}".format(tx_fa[i].split('/')[-1], tx_fa[j].split('/')[-1]), file=sys.stderr)
                    p2 = subprocess.Popen(cmd2, stdout=out2_this)
                    stdout, stderr = p2.communicate()
                    out2_this.close()
                else:
                    print("Loading alignments from {}...".format(outfn2_this), file=sys.stderr)

                in2_this = open(outfn2_this, 'r')
                for line in in2_this:
                    line = line.rstrip()
                    dat = line.split('\t')
                    if len(dat) >= 12:
                        query = dat[0]
                        ref = dat[5]
                        nmatch = int(dat[9])
                        tot = int(dat[10])
                        lquery = int(dat[1])
                        lref = int(dat[6])
                        denom = lquery
                        if lref < lquery:
                            denom = lref
                        score = nmatch / denom
                        distmat[(i,j)] += score
                        distmat_count[(i,j)] += 1
                        
                        if query not in graphs[i]:
                            graphs[i][query] = {}
                        if j not in graphs[i][query]:
                            graphs[i][query][j] = {}
                        graphs[i][query][j][ref] = score
                        edges_rev[j].append((i, query))
                        #key = (j, ref)
                        #if key not in graphs[i][query] or score > graphs[i][query][key]:
                        #    graphs[i][query][key] = score
                in2_this.close()
        
        for key in distmat:
            distmat[key] /= distmat_count[key]
        tx_filt = tx_upgma(distmat, graphs, out, edges_rev)

        #tx_filt = proc_graphs(graphs, out)

    # Write out new FASTA files with only filtered transcripts.
    tx_fa_filt = []
    for i, fn in enumerate(tx_fa):
        fn_new = '{}.{}.filt.fa'.format(out, fn.split('/')[-1])
        tx_fa_filt.append(fn_new)
        fout = open(fn_new, 'w')
        fin = None
        gz = False
        if is_gz_file(fn):
            fin = gzip.open(fn, 'r')
            gz = True
        else:
            fin = open(fn, 'r')
        cur_id = None
        cur_print = False
        for line in fin:
            if gz:
                line = line.decode().rstrip()
            else:
                line = line.rstrip()
            if line[0] == '>':
                cur_id = line[1:].split()[0]
                if cur_id in tx_filt[i]:
                    cur_print = True
                else:
                    cur_print = False
            if cur_print:
                print(line, file=fout)
        fin.close()
        fout.close()
    
    return (tx_fa_filt, dbs)

def main(args):
    options = parse_args()
    
    files_cleanup = []
    tx_fa = []

    if options.gtf is not None and len(options.gtf) > 0:
        # Extract transcripts from FASTA
        tx_fa = extract_tx(options.fasta, options.gtf, options.out, options.gffread)    
        files_cleanup = tx_fa
    else:
        tx_fa = options.fasta
    
    tx_fa_filt = []
    if not options.no_orthology:
        tx_fa_filt, extrafiles = filter_orthology(tx_fa, options.out, options.mm2)
        files_cleanup += tx_fa_filt
        files_cleanup += extrafiles
    else:
        tx_fa_filt = tx_fa

    # Count k-mers
    ktabs = []
    for idx, tx in enumerate(tx_fa_filt):
        count_kmers(tx, '{}.fastk.{}'.format(options.out, idx), options.k, options.FastK)
        ktabs.append('{}.fastk.{}.ktab'.format(options.out, idx))
    
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
