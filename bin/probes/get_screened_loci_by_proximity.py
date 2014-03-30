#!/usr/bin/env python
# encoding: utf-8
"""
File: get_screened_loci_by_proximity.py
Author: Brant Faircloth

Created by Brant Faircloth on 23 September 2013 13:09 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import argparse
import random
from Bio import SeqIO
from collections import defaultdict
from phyluce.helpers import is_file
from bx.intervals.cluster import ClusterTree

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Given a FASTA file of properly formatted sequences (and header), screen for distance""")
    parser.add_argument(
            "--input",
            required=True,
            type=is_file,
            help="""The FASTA file of input sequences"""
        )
    parser.add_argument(
            "--output",
            required=True,
            help="""The output FASTA file of filtered sequences"""
        )
    parser.add_argument(
            "--distance",
            type=int,
            default=10000,
            help="""The minimum distance on which to filter"""
        )
    parser.add_argument(
            "--length",
            type=int,
            default=None,
            help="""The minimum length of sequences to filter"""
        )
    return parser.parse_args()


def get_positions_from_coords(s):
    chromo, coords = s.split(':')
    start, end = [int(i) for i in coords.split('-')]
    return chromo, start, end


def main():
    args = get_args()
    # build a cluster tree for the loci, allowing the user to set the clustering distance
    # loci within args.distance
    trees = defaultdict(lambda:ClusterTree(args.distance, 2))
    # fill cluster tree with data.  name_int converts the int part of the locus name to an int
    with open(args.input, 'rU') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            name_int = int(record.id.split('_')[1])
            chromo, start, end = get_positions_from_coords(record.description.split('|')[1])
            trees[chromo].insert(start, end, name_int)
    bad = []
    # for each cluster, get the member names
    for chromo, cluster in trees.iteritems():
        clustered = cluster.getlines()
        if clustered != []:
            # if there are a cluster of loci, randomly choose one - ignore the rest
            choice = random.choice(clustered)
            clustered.remove(choice)
            bad.extend(clustered)
    bad = set(bad)
    with open(args.input, 'rU') as infile:
        with open(args.output, 'w') as outf:
            for record in SeqIO.parse(infile, 'fasta'):
                name_int = int(record.id.split('_')[1])
                chromo, start, end = get_positions_from_coords(record.description.split('|')[1])
                # Drop anything with an ambiguous base
                if 'N' in record.seq.upper():
                    pass
                # Filter by locus length, if so desired
                elif args.length is not None and not name_int in bad:
                    if end - start >= args.length:
                        outf.write(record.format('fasta'))
                # If we're not filtering by locus length, output everything that is
                # Not in a cluster
                elif args.length is None and not name_int in bad:
                    outf.write(record.format('fasta'))

if __name__ == '__main__':
    main()
