#!/usr/bin/env python
# encoding: utf-8
"""
File: get_informative_sites.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 August 2012 21:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import glob
import sqlite3
import argparse
import multiprocessing
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from collections import Counter
from phyluce.helpers import is_dir, FullPaths, get_file_extensions


import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Count the number of informative sites in a given set of alignment"""
        )
    parser.add_argument(
            'input',
            type=is_dir,
            action=FullPaths,
            help="""The directory containing the alignment files"""
        )
    parser.add_argument(
            '--output',
            type=str,
            default='output',
            help="""The output filename (without extension - code will add .sqlite)"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of cores to use.""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def get_informative_sites(count):
    #pdb.set_trace()
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    sufficient_sites = len(count)
    if sufficient_sites >= 2:
        sufficient_sequences = sum([1 for i in count.values() if i >= 2])
        if sufficient_sequences >= 2:
            return True
    return False


def main():
    args = get_args()
    results = {}
    for f in get_files(args.input, args.input_format):
        aln = AlignIO.read(f, args.input_format)
        name = os.path.basename(f)
        informative_sites = []
        for idx in xrange(aln.get_alignment_length()):
            col = aln[:, idx]
            count = Counter(col)
            if get_informative_sites(count):
                informative_sites.append(1)
            else:
                informative_sites.append(0)
        #pdb.set_trace()
        results[name] = {
                'length': aln.get_alignment_length(),
                'informative': sum(informative_sites)
            }
    for locus, v in results.iteritems():
        print "{0}\t{1}\t{2}".format(locus, v['informative'], v['length'])

if __name__ == '__main__':
    main()
