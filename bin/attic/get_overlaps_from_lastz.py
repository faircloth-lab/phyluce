#!/usr/bin/env python
# encoding: utf-8
"""
File: get_bed_from_lastz.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 June 2012 15:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Create a BED file from LASTZ alignments

"""

import os
import random
import argparse
import ConfigParser
from phyluce import lastz
from collections import defaultdict
from phyluce.helpers import FullPaths, is_file
from seqtools.sequence import fasta
from operator import itemgetter

from bx.intervals.intersection import Interval, IntervalTree

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Output a BED file from a LASTZ alignment""")
    parser.add_argument(
            "lastz",
            type=is_file,
            action=FullPaths,
            help="""The input lastz file."""
        )
    parser.add_argument(
            "fasta",
            type=is_file,
            action=FullPaths,
            help="""The input fasta file (used for lastz alignments)."""
        )
    parser.add_argument(
            "output",
            type=argparse.FileType('wb'),
            help="""The output config file"""
        )
    return parser.parse_args()


def add_new_locus(match, overlappers, chromo):
    intersecter = IntervalTree()
    locus = match.name2.split('|')[1]
    intersecter.insert(match.zstart1, match.end1, locus)
    overlappers[chromo][locus] = intersecter
    return overlappers


def get_longest_of_overlapping_loci(all_groups, seq_lengths):
    ''' when selecting only a single locus of the cluster, get the longest.  in case of a tie,
    make a random selection of the longest loci'''
    keep = []
    for group in all_groups:
        lengths = {locus:seq_lengths[locus] for locus in group}
        max_length = max(lengths.values())
        longest = {k:v for k,v in lengths.iteritems() if v == max_length}
        keep.append(random.choice(longest.keys()))
    return keep


def main():
    args = get_args()
    uce_loci = []
    # get lengths of loci
    seq_lengths = {}
    for seq in fasta.FastaReader(args.fasta):
        name = seq.identifier.split('|')[1]
        uce_loci.append(name)
        seq_lengths[name] = len(seq.sequence)
    overlappers = defaultdict(dict)
    names = defaultdict(list)
    coords = {}
    for match in lastz.Reader(args.lastz, long_format=True):
        locus = match.name2.split('|')[1]
        chromo = match.name1
        coords[locus] = (match.zstart1, match.end1)
        for pmatch, span in overlappers[chromo].iteritems():
            if locus == 'chr5_10696_s' and pmatch == 'chr13_710_s':
                pdb.set_trace()
            overlap = span.find(match.zstart1, match.end1)
            if overlap:
                overlappers[chromo][pmatch].insert(match.zstart1, match.end1, locus)
                names[pmatch].append(locus)
                break
        else:
            overlappers = add_new_locus(match, overlappers, chromo)
    overlapping_loci = []
    all_groups = []
    for k, v in names.iteritems():
        # group loci into overlapping clusters
        base = [k]
        base.extend(v)
        all_groups.append(base)
        # get list of "bad loci" so we can determine non-overlappers
        overlapping_loci.append(k)
        overlapping_loci.extend(v)
    pdb.set_trace()
    non_overlapping_loci = set(uce_loci).difference(set(overlapping_loci))
    # generate output in config-file format:
    config = ConfigParser.RawConfigParser()
    config.add_section('Non-overlapping loci')
    for locus in list(non_overlapping_loci):
        config.set('Non-overlapping loci', locus, seq_lengths[locus])
    longest_of_overlapping = get_longest_of_overlapping_loci(all_groups, seq_lengths)
    config.add_section('Longest loci of group')
    for locus in longest_of_overlapping:
        config.set('Longest loci of group', locus, seq_lengths[locus])
    config.add_section('Superlocus groups')
    for c, group in enumerate(all_groups):
        # order loci by start position
        starts = [(name,coords[name][0], coords[name][1]) for name in group]
        starts = sorted(starts, key=itemgetter(1))
        sorted_names = [n[0] for n in starts]
        print starts
        #pdb.set_trace()
        config.set('Superlocus groups', "Group{0}".format(c), ','.join(sorted_names))
    config.write(args.output)


if __name__ == '__main__':
    main()

    [('chr1_1302_s', 36426950), ('chr1_1303_s', 36426950), ('chr3_9755_s', 36426950)]

#chr9_13130_s
#chr17_2784_s