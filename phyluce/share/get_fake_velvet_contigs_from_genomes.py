#!/usr/bin/env python
# encoding: utf-8

"""
get_fake_velvet_contigs_from_genomes.py

Created by Brant Faircloth on 11 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import glob
import argparse
import textwrap
import bx.seq.twobit
from seqcap.lib import lastz
from operator import itemgetter
from collections import defaultdict
from tools.sequence import transform
from uce_helpers import get_name
from uce_helpers import get_dupes

import pdb



def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('lastz', help='The lastz output')
    parser.add_argument('genome', help='The 2bit genome file')
    parser.add_argument('--fasta', help='The file in which to store the fasta output',
        type=argparse.FileType('w'))
    parser.add_argument('--bed', help = 'Return a BED file of the matching regions ± 500 bp flank',
        type=argparse.FileType('w'))
    parser.add_argument('--flank', help = 'The flank length to add', default = 500, type = int)
    parser.add_argument('--name-components', dest = 'components', help = 'number of parts in the name', default = 2, type = int)
    parser.add_argument('--splitchar', help = 'The name character on which to split', default = "_", type = str)
    parser.add_argument('--dupefile', help='The path to a lastz file of lastz-against-self results')
    parser.add_argument('--fish', help='If running against fish probes', action='store_true', default = False)
    parser.add_argument('--verbose', action='store_true', default = False)
    return parser.parse_args()

def run_checks(k, v, probes, verbose = True):
    try:
        # make sure we have as many matches as probes
        assert len(v) == probes[k]              # more matches than expected
        assert len(set([i[0] for i in v])) == 1 # matches to more than one chromosome
        assert len(set([i[1] for i in v])) == 1 # multiple match directions
        return True
    except AssertionError:
        if verbose:print "Probe: ", k
        result_string = ["{0}:{1}-{2}|{3}".format(i[0], i[2], i[3], i[1]) for i in v]
        if verbose:print "\t\tExpected hits: {0}\n\t\tObserved Hits: {1}\n\t\t\t{2}".format(probes[k], len(v), '\n\t\t\t'.join(result_string))
        return False

def get_matches(lastz_file, splitchar, components, fish = False):
    matches = defaultdict(list)
    probes = defaultdict(int)
    for lz in lastz.Reader(lastz_file, long_format = True):
        # skip silly hg19 mhc haplotypes
        if "hap" in lz.name1:
            print "Skipping: ", lz.name1
        else:
            if not fish:
                uce_name = get_name(lz.name2, "|", 1)
                probe_number = int(lz.name2.split(':')[-1])
            else:
                uce_name = get_name(lz.name2, "_", 1)
                # add 1 because fish probe indexing starts @ 0
                probe_number = int(lz.name2.split('|')[1].split('_')[1]) + 1
            #pdb.set_trace()
            if probe_number > probes[uce_name]:
                probes[uce_name] = probe_number
            matches[uce_name].append([get_name(lz.name1, splitchar = splitchar, items = components), lz.strand2, lz.zstart1, lz.end1])
    return matches, probes
    
def main():
    args = get_args()
    if args.dupefile:
        dupes = get_dupes(args.dupefile)
    else:
        dupes = None
    matches, probes = get_matches(args.lastz, args.splitchar, args.components, args.fish)
    #unique_matches = sum([1 for uce, map_pos in matches.iteritems() if len(map_pos) == probes[uce]])
    if args.fasta:
        tb = bx.seq.twobit.TwoBitFile(file(args.genome))
    count = 0
    for k,v in matches.iteritems():
        skip = False
        if len(v) > 1:
            if run_checks(k, v, probes):
                # sort by match position
                v_sort = sorted(v, key = itemgetter(2))
                start, end = v_sort[0][2], v_sort[-1][3]
                diff = end - start
                # ensure our range is less than N(probes) * probe_length - this
                # still gives us a little wiggle room because probes are ~ 2X tiled
                if diff > (probes[k] * 140):
                    skip = True
                    if args.verbose:
                        print "range longer than expected"
                else:
                    chromo = v[0][0]
                    strand = v[0][1]
            else:
                skip = True
        elif k in dupes:
            skip = True
            print "{0} is in dupefile".format(k)
        else:
            chromo, strand, start, end = v[0]
        if not skip and args.fasta:
            # slice out region + flank
            try:
                slc = tb[chromo][start - args.flank:end + args.flank]
            except:
                pdb.set_trace()
            # strip Ns from both ends
            slc = slc.strip('N')
            # reverse any strands where necessary
            if not strand == '+':
                slc = transform.DNA_reverse_complement(slc)
            if len(slc) != 0:
                args.fasta.write(">Node_{0}_length_{1}_cov_100\n{2}\n".format(count, len(slc), '\n'.join(textwrap.wrap(slc))))
        if not skip and args.bed:
            args.bed.write("{0} {1} {2} {3} 1000 {4}\n".format(chromo, start - args.flank, end + args.flank, k, strand))
        count += 1
        #pdb.set_trace()
        

if __name__ == '__main__':
    main()