#!/usr/bin/env python
# encoding: utf-8
"""
File: explode_get_fastas_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 12 June 2013 23:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description: Given the fasta file produced by get_fastas_from_match_counts.py,
explode that concatenated fasta file into many separate files, 1 per locus.

"""

import os
import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Explode the fasta file produced by get_fastas_from_match_counts into single files"""
    )
    parser.add_argument(
        "--input",
        required=True,
        action=FullPaths,
        help="""The input fasta file to explode"""
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        action=FullPaths,
        help="""The output directory to create and in which to store the fastas""",
    )
    parser.add_argument(
        "--by-taxon",
        action="store_true",
        default=False,
        help="""Split file by taxon and not by locus""",
    )
    return parser.parse_args()


def main():
    args = get_args()
    seqdict = defaultdict(list)
    if os.path.isdir(args.output_dir):
        response = raw_input("{} exists.  Add results to directory [Y/n]? ".format(args.output_dir))
        if response == 'Y':
            pass
        else:
            sys.exit()
    else:
        os.makedirs(args.output_dir)
    print "Reading fasta..."
    if not args.by_taxon:
        with open(args.input, 'rU') as input:
            for seq in SeqIO.parse(input, 'fasta'):
                uce = seq.id.split('_')[0]
                seqdict[uce].append(seq)
    elif args.by_taxon:
        with open(args.input, 'rU') as input:
            for seq in SeqIO.parse(input, 'fasta'):
                taxon = '-'.join(seq.id.split('_')[1:])
                seqdict[taxon].append(seq)
    print "Writing fasta..."
    for uce, sequences in seqdict.iteritems():
        with open(os.path.join(args.output_dir, '{}.unaligned.fasta'.format(uce)), 'w') as outf:
            for sequence in sequences:
                outf.write(sequence.format('fasta'))


if __name__ == '__main__':
    main()
