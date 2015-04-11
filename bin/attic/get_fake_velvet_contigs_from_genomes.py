#!/usr/bin/env python
# encoding: utf-8

"""
get_fake_velvet_contigs_from_genomes.py

Created by Brant Faircloth on 11 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import re
import argparse
import textwrap
import bx.seq.twobit
from operator import itemgetter
from seqtools.sequence import transform
from phyluce.helpers import get_dupes, get_matches, run_checks

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('lastz',
            help='The lastz output')
    parser.add_argument('genome',
            help='The 2bit genome file')
    parser.add_argument('--fasta',
            help='The file in which to store the fasta output',
            type=argparse.FileType('w')
        )
    parser.add_argument('--bed',
            help='Return a BED file of the matching regions ± 500 bp flank',
            type=argparse.FileType('w')
        )
    parser.add_argument('--flank',
            help='The flank length to add',
            default=500,
            type=int
        )
    parser.add_argument('--name-components',
            dest='components',
            help='number of parts in the name',
            default=2,
            type=int
        )
    parser.add_argument('--splitchar',
            default="_",
            type=str,
            help='The name character on which to split'
        )
    parser.add_argument('--dupefile',
            help='The path to a lastz file of lastz-against-self results')
    parser.add_argument('--fish',
            default=False,
            help='If running against fish probes', action='store_true')
    parser.add_argument('--verbose',
            default=False,
            action='store_true')
    return parser.parse_args()


def quality_control_matches(matches, probes, dupes, k, v, verbose=False):
    """check to make sure we don't get any more matches than expected
    and that matches are reasonably close to each other on their respective
    chromos"""
    skip = False
    chromo, strand, start, end = None, None, None, None
    if len(v) > 1:
        if run_checks(k, v, probes):
            # sort by match position
            v_sort = sorted(v, key=itemgetter(2))
            start, end = v_sort[0][2], v_sort[-1][3]
            diff = end - start
            # ensure our range is less than N(probes) * probe_length - this
            # still gives us a little wiggle room because probes are ~ 2X tiled
            if diff > (probes[k] * 140):
                skip = True
                if verbose:
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
    return chromo, strand, start, end, skip


def snip_if_many_N_bases(regex, chromo, seq, uce):
    """Some genome builds contain long runs of Ns.  Since we're
    slicing reads from these genomes, sometimes these slices contains
    giant runs of Ns.  Remove these by finding the UCE and trimming out
    from the middle to retain the UCE while removing the Ns"""
    # find uce in seq
    uce_start = seq.find(uce)
    uce_end = uce_start + len(uce)
    # slice front
    seq_slice = seq[:uce_start]
    # reverse it - we want first occurence moving 5'
    # from uce start (so first going backwards)
    seq_slice = seq_slice[::-1]
    # search for Ns
    r = regex.search(seq_slice)
    if r:
        new_start = len(seq_slice) - r.start()
    else:
        new_start = 0
    # slice rear
    seq_slice = seq[uce_end:]
    r = regex.search(seq_slice)
    if r:
        new_end = uce_end + r.start()
    else:
        new_end = len(seq)
    seq = seq[new_start:new_end]
    print "{0} trimmed for > 20 N bases".format(chromo)
    return seq


def prep_and_write_fasta(tb, regex, fasta, chromo, strand, start, end, count, flank=500):
    """write out our sliced sequences to fasta"""
    # slice out region + flank
    chromo = chromo.lstrip('>')
    try:
        slc = tb[chromo][start - flank:end + flank]
    except:
        pdb.set_trace()
    # strip Ns from both ends
    slc = slc.strip('N')
    # deal with large N insertions
    result = regex.search(slc)
    if result:
        uce_slice = tb[chromo][start:end]
        slc = snip_if_many_N_bases(regex, chromo, slc, uce_slice)
    # reverse any strands where necessary
    if not strand == '+':
        slc = transform.DNA_reverse_complement(slc)
    if len(slc) != 0:
        fasta.write(">Node_{0}_length_{1}_cov_100\n{2}\n".format(count, len(slc), '\n'.join(textwrap.wrap(slc))))


def main():
    args = get_args()
    regex = re.compile("[N,n]{20,}")
    if args.dupefile:
        dupes = get_dupes(args.dupefile, longfile=False)
    else:
        dupes = None
    matches, probes = get_matches(args.lastz, args.splitchar, args.components, args.fish)
    #unique_matches = sum([1 for uce, map_pos in matches.iteritems() if len(map_pos) == probes[uce]])
    if args.fasta:
        tb = bx.seq.twobit.TwoBitFile(file(args.genome))
    count = 0
    for k, v in matches.iteritems():
        chromo, strand, start, end, skip = quality_control_matches(matches, probes, dupes, k, v, args.verbose)
        if not skip and args.fasta:
            prep_and_write_fasta(tb, regex, args.fasta, chromo, strand, start, end, count, args.flank)
        if not skip and args.bed:
            args.bed.write("{0} {1} {2} {3} 1000 {4}\n".format(chromo, start - args.flank, end + args.flank, k, strand))
        count += 1
        #pdb.set_trace()
    args.fasta.close()


if __name__ == '__main__':
    main()
