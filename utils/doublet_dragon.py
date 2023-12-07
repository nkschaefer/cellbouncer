#! /usr/bin/env python3
import sys
import os
import math
import argparse
from collections import defaultdict, Counter
import numpy as np
import gzip
"""
Given one or more barcode assignment files from species or individual demultiplexing,
estimates a single overall doublet rate in the data set.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indiv_assignments", "-i", nargs="+", help="Barcode assignment files \
from multiple runs of demux_mt or demux_vcf on data from the same library. For example, if \
you pooled individuals from two species, you could provide demux_mt results from demultiplexing \
individuals from species 1 and species 2.", required=False)
    parser.add_argument("--species_assignments", "-s", help="Barcode assignment file from species \
demultiplexing. If you provide individual demultiplexing data from the same library, then only \
barcodes present in the individual files (which are assumed to have passed filtering) will be \
considered.", required=False, default=None)
    parser.add_argument("--name", "-n", help="The name of the library", required=False, \
        default=None)
    parsed = parser.parse_args()
    if (parsed.indiv_assignments is None or len(parsed.indiv_assignments) == 0) and \
        parsed.species_assignments is None:
        print("ERROR: you must provide at least one input file", file=sys.stderr)
        parser.print_help()
        exit(1) 
    return parser.parse_args()

"""
Create system of equations in matrix format, to be solved using
Newton-Raphson method
"""
def get_params(inds, sing_counts, doub_counts, tot_counts):
    # Set initial guesses
    # Last element is doublet rate
    params = [ 1/len(inds) ] * len(inds)
    params.append(0.25)

    # Store rhs of equations
    rhs = []

    for i in range(0, len(inds)):
        id1 = inds[i]
        frac = sing_counts[id1]/tot_counts
        rhs.append(frac)
    
    # Now add in doublets
    for i in range(0, len(inds)-1):
        id1 = inds[i]
        for j in range(i+1, len(inds)):
            id2 = inds[j]
            doubkey = id1 + "+" + id2
            if id2 < id1:
                doubkey = id2 + "+" + id1
            frac = doub_counts[doubkey] / tot_counts
            rhs.append(frac)
    
    return (params, rhs)

"""
Use Newton-Raphson method to solve system of equations.
"""
def solve(params, rhs, inds):
    delta_thresh = 0.0001
    delta = 999
    its = 0
    while delta > delta_thresh:
        its += 1
        # Compute the Jacobian
        J = []
        
        # What is the current estimate of doublet rate?
        doub_rate = params[-1]
        
        # Add "singlet" rows to Jacobian
        for i in range(0, len(inds)):
            J_row = [0.0] * (len(inds) + 1)

            # Current guess for self matching parameter
            self_rate = params[i]

            # Partial derivative with respect to self
            J_row[i] = 2.0 * doub_rate * self_rate - doub_rate + 1

            # Partial derivative with respect to doublet rate
            J_row[-1] = self_rate * self_rate - self_rate
            
            J.append(J_row)
        
        # Add doublet rows to Jacobian
        for i in range(0, len(inds)-1):
            i_rate = params[i]
            for j in range(i+1, len(inds)):
                j_rate = params[j]
                
                J_row = [0.0] * (len(inds) + 1)
                J_row[i] = doub_rate * j_rate
                J_row[j] = doub_rate * i_rate
                J_row[-1] = i_rate * j_rate
                
                J.append(J_row)
        
        # Add a row that ensures all terms other than doublet rate
        #  sum to 1
        J_row = []
        for idx in range(0, len(params)-1):
            # Partial derivative wrt each variable is 1, except doublet rate (0)
            J_row.append(1.0)
        J_row.append(0)
        J.append(J_row)

        J = np.array(J)

        # Solve negative F vector 
        neg_f = []
        k = 0
        for i in range(0, len(inds)):
            i_rate = params[i]
            f = doub_rate*i_rate*i_rate - doub_rate*i_rate + i_rate - rhs[i]
            neg_f.append(-f)
            k += 1
        for i in range(0, len(inds)-1):
            i_rate = params[i]
            for j in range(i+1, len(inds)):
                j_rate = params[j]
                f = doub_rate*i_rate*j_rate - rhs[k]
                neg_f.append(-f)
                k += 1
        
        # Handle last term that ensures proportions (other than doublet rate)
        # sum to 1
        ratesum = 0.0
        for idx in range(0, len(inds)):
            ratesum += params[idx]
        neg_f.append(-(ratesum - 1.0))
        k += 1

        neg_f = np.array(neg_f).reshape(-1,1)
        
        res = np.linalg.lstsq(J, neg_f, rcond=None)
        delta_vec = res[0].reshape(1,-1)[0]
        
        for idx, elt in enumerate(delta_vec):
            # This is the delta on the current coefficient
            
            if idx == len(res[0])-1:
                # Track the delta on the doublet rate, which is the main parameter of interest.
                delta = abs(elt / params[idx])
            
            params[idx] = params[idx] + delta_vec[idx]
     
"""
Read in a file of barcode -> individual (or species) assignments

bc_filter option tells whether to use valid_bc is a filter when reading 
    the file (True) or whether to populate valid_bc with bcs in the 
    file (False).
"""
def parse_assignments(filename, sing_counts, doub_counts, valid_bc, \
    bc_filter=False):
    llrsum = 0.0 
    tot_counts = 0
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        bc = dat[0].split('-')[0]
        line_pass = False
        if not bc_filter:
            valid_bc.add(bc)
            line_pass = True
        else:
            if bc in valid_bc:
                line_pass = True
        if line_pass:
            llr = float(dat[3])
            llrsum += llr
            if dat[2] == 'S':
                sing_counts[dat[1]] += 1
            elif dat[2] == 'D':
                doub_counts[dat[1]] += 1
            tot_counts += 1
    f.close()
    return (tot_counts, llrsum)

"""
Given Counters storing counts of each type of individual and doublet
combination, infer the overall doublet rate using the Newton-Raphson 
method
"""
def infer_doublet_rate(sing_counts, doub_counts, tot_counts, indv_p):
    if len(sing_counts) == 0 or len(doub_counts) == 0:
        # Can't compute
        return -1
    inds = sorted(list(sing_counts.keys()))
    
    # Each "singlet" count is actually a combination of singlets and 
    # within-individual or within-species doublets.

    # If the fraction of the data set that is "singlets" of ID x is f
    #   and the probability of sampling a cell from ID x is x, and the 
    #   doublet rate is d, then 
    # f = (1-d)*x + d*x*x = d*x^2 - d*x + x

    # If the fraction of the data set that is doublets of IDs x and y is
    #   g and the probability of sampling a cell from ID y is y, then
    # g = d*x*y
    
    # We will use these to set up a system of equations and solve for all
    #   x, y as well as d using the Newton-Raphson method.
    
    params, rhs = get_params(inds, sing_counts, doub_counts, tot_counts)        
    solve(params, rhs, inds)
    
    fsum = 0
    for idx in range(0, len(params)-1):
        fsum += params[idx]
        indv_p[inds[idx]] = params[idx]
    doub_rate = params[-1]
    
    return doub_rate

def main(args):
    options = parse_args()
    
    # If the user provided species and individual data, store all individual
    # barcodes to use as a filtered list for species assignments
    valid_bc = set([])
    
    llrsums = []
    llrtot = 0.0
    doub_rates = []
    doub_rates_keyed = {}

    indv_props = {}

    for assnfile in options.indiv_assignments:
        
        # Read file
        sing_counts = Counter()
        doub_counts = Counter()
        tot_counts, llrsum = parse_assignments(assnfile, sing_counts, doub_counts, valid_bc, \
            bc_filter=False)
        
        # Each data set can be weighted by the sum of its LLRs
        llrsums.append(llrsum)
        llrtot += llrsum

        # Estimate the doublet rate
        indv_p = {}
        doub_rate = infer_doublet_rate(sing_counts, doub_counts, tot_counts, indv_p)
        
        assnfile_clean = assnfile.split('/')[-1].split('.assignments')[0]
        indv_props[assnfile_clean] = indv_p
        doub_rates.append(doub_rate)
        doub_rates_keyed[assnfile_clean] = doub_rate

    if options.species_assignments is not None:
        
        sing_counts = Counter()
        doub_counts = Counter()
        tot_counts, llrsum = parse_assignments(options.species_assignments, sing_counts, doub_counts, \
            valid_bc, bc_filter=True)

        llrsums.append(llrsum)
        llrtot += llrsum

        indv_p = {}
        doub_rate = infer_doublet_rate(sing_counts, doub_counts, tot_counts, indv_p)
        
        assnfile_clean = assnfile.split('/')[-1].split('.assignments')[0]
        indv_props[assnfile_clean] = indv_p
        doub_rates.append(doub_rate)
        doub_rates_keyed[assnfile_clean] = doub_rate

    # Compute overall doublet rate as a weighted average
    doub_rate = 0.0
    for idx in range(0, len(doub_rates)):
        doub_rate += (llrsums[idx] / llrtot) * doub_rates[idx]

    name = 'all'
    if options.name is not None:
        name = options.name
    print("{}\tdoublet_rate\t{:.3f}".format(name, doub_rate))
    
    for fn in sorted(doub_rates_keyed.keys()):
        print("{}\tdoublet_rate\t{:.3f}".format(fn, doub_rates_keyed[fn]))
        
        # Sort other rates
        for tup in sorted(indv_props[fn].items(), key=lambda x: -x[1]):
            print("{}\t{}\t{:.3f}".format(fn, tup[0], tup[1]))

   
if __name__ == '__main__':
    sys.exit(main(sys.argv))

