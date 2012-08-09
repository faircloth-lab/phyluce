#!/usr/bin/env python
# encoding: utf-8
"""
File: get_all_locus_lengths.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 August 2012 19:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import shutil
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

import pdb

def get_args():
    parser = argparse.ArgumentParser(
            description="""Output the lengths of alignments"""
        )
    parser.add_argument(
            'input',
            type=is_dir,
            action=FullPaths,
            help='The directory containing the alignment files'
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument(
            "--containing-data-for",
            dest="containing",
            type=str,
            help="""Output alignments that contain data for a taxon"""
        )
    parser.add_argument(
            "--min-length",
            dest="min_length",
            type=int,
            default=0,
            help="""Filter out alignments longer than --min-length"""
        )
    parser.add_argument(
            "--output",
            type=is_dir,
            action=FullPaths,
            help="""Place alignments meeting criteria in an output folder"""
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def main():
    args = get_args()
    files = get_files(args.input, args.input_format)
    for f in files:
        try:
            aln = AlignIO.read(f, args.input_format)
            if args.containing:
                good = False
                for taxon in aln:
                    if taxon.id == args.containing:
                        if (set(taxon.seq) == set('-')) or (set(taxon.seq) == set('?')):
                            good = False
                        else:
                            good = True
            else:
                good = True
            if good and aln.get_alignment_length() >= args.min_length:
                print "{0}\t{1}".format(os.path.basename(f), aln.get_alignment_length())
            if args.output and good and aln.get_alignment_length() >= args.min_length:
                name = os.path.basename(f)
                shutil.copy(f, os.path.join(args.output, name))
        except ValueError, e:
            if e.message == 'No records found in handle':
                print 'No records found in {0}'.format(os.path.basename(f))
            else:
                raise ValueError('Something is wrong with alignment {0}'.format(os.path.basename(f)))


if __name__ == '__main__':
    main()

