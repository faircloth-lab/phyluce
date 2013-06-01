#!/usr/bin/env python
# encoding: utf-8

"""

remove_taxa_from_alignments.py

Created by Brant Faircloth on 21 September 2011 11:44 PDT (-0700).
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

PURPOSE:  Remove taxa passed on CLI from a folder of nexus alignment
files.

"""

import pdb
import os
import sys
import glob
import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment
from phyluce.helpers import is_dir

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument('input', help="The input directory containing nexus files", type=is_dir)
    parser.add_argument('output', help="The directory in which to store the output files", type=is_dir)
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--exclude',
        type=str,
        default=[],
        nargs='+',
        help='Taxa to exclude')
    group.add_argument('--include',
        type=str,
        default=[],
        nargs='+',
        help='Taxa to include')
    return parser.parse_args()

def get_samples_to_run(args, all_names):
    """docstring for get_samples_to_run"""
    if args.exclude:
        return set([name for name in all_names if name not in args.exclude])
    elif args.include:
        return set([name for name in all_names if name in args.include])
    else:
        return all_names
    
def get_all_taxon_names(nexus_files):
    taxa = set()
    for align_file in nexus_files:
        for align in AlignIO.parse(align_file, "nexus"):
            for taxon in list(align):
                #pdb.set_trace()
                taxa.add(taxon.name)
    return taxa

def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    args = get_args()
    nexus_files = get_files(args.input)
    taxa = get_all_taxon_names(nexus_files)
    taxa_to_keep = get_samples_to_run(args, taxa)
    for count, align_file in enumerate(nexus_files):
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        for align in AlignIO.parse(align_file, "nexus"):
            for taxon in list(align):
                if taxon.name in taxa_to_keep:
                    new_align.add_sequence(taxon.name, str(taxon.seq))
        outf = os.path.join(args.output, os.path.basename(align_file))
        AlignIO.write(new_align, open(outf, 'w'), 'nexus')
        print count


if __name__ == '__main__':
    main()
