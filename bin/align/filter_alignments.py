#!/usr/bin/env python
# encoding: utf-8
"""
File: filter_alignments.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 August 2012 19:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Given a folder of alignments, output the name and length.
Filter on presence of data for a taxon with --containing-data-for and/or
length with --min-length.  If --output, copy resulting files into a new
directory.

Usage: python filter_alignments phylip-with-gaps \
    --input-format phylip \
    --containing-data-for pol_sen \
    --output phylip-with-gaps-polypterus \
    --min-length 50
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
            nargs='+',
            type=str,
            help="""Output only alignments --containing-data-for a taxon"""
        )
    parser.add_argument(
            "--min-length",
            dest="min_length",
            type=int,
            default=0,
            help="""Filter out alignments shorter than --min-length"""
        )
    parser.add_argument(
            "--min-taxa",
            dest="min_taxa",
            type=int,
            default=0,
            help="""Filter out alignments with fewer than --min-taxa"""
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


def align_contains_taxa(args, aln):
    containing = False
    #pdb.set_trace()
    for taxon in aln:
        if taxon.id in args.containing:
            # if there's no seq data, then the focus taxon is not really present
            if (set(taxon.seq) == set('-')) or (set(taxon.seq) == set('?')):
                containing = False
            else:
                containing = True
            break
        else:
            pass
    return containing


def align_min_length(args, aln):
    if aln.get_alignment_length() >= args.min_length:
        length = True
    else:
        length = False
    return length


def align_min_taxa(args, aln):
    count = 0
    # remove taxa having only missing data designators
    for taxon in aln:
        if (set(taxon.seq) == set('-')) or (set(taxon.seq) == set('?')):
            pass
        else:
            count += 1
    if count >= args.min_taxa:
        taxa = True
    else:
        taxa = False
    return taxa


def main():
    args = get_args()
    files = get_files(args.input, args.input_format)
    print "Good Alignments\n"
    for f in files:
        try:
            aln = AlignIO.read(f, args.input_format)
            if args.containing:
                containing = align_contains_taxa(args, aln)
            else:
                containing = True
            if args.min_length:
                length = align_min_length(args, aln)
            else:
                length = True
            if args.min_taxa:
                taxa = align_min_taxa(args, aln)
            else:
                taxa = True
            if containing and taxa and length:
                print "{0}".format(os.path.basename(f))
            if containing and taxa and length and args.output:
                name = os.path.basename(f)
                shutil.copy(f, os.path.join(args.output, name))
        except ValueError, e:
            if e.message == 'No records found in handle':
                print 'No records found in {0}'.format(os.path.basename(f))
            else:
                raise ValueError('Something is wrong with alignment {0}'.format(os.path.basename(f)))

if __name__ == '__main__':
    main()
