#!/usr/bin/env python
# encoding: utf-8

"""
change_taxa_names.py

Created by Brant Faircloth on 22 September 2010 12:48 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  Remove the UCE locus name from nexus alignments.

USAGE:  python remove_locus_name_from_nexus_lines.py \
    --input my/input/folder/nexus \
    --output my/input/folder/nexus-renamed
"""


import os
import re
import glob
import argparse
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from phyluce.helpers import FullPaths, is_dir

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Remove the UCE locus name from nexus alignments.""")
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""The input directory containing nexus files to filter"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""The output directory to hold the converted nexus files""",
        )
    parser.add_argument(
            "--taxa",
            type=int,
            default=None,
            help="""The expected number of taxa in all alignments""",
        )
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex*'))


def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.input)
    all_taxa = set([])
    for count, f in enumerate(files):
        new_align = MultipleSeqAlignment([], generic_dna)
        for align in AlignIO.parse(f, 'nexus'):
            for seq in list(align):
                fname = os.path.splitext(os.path.basename(f))[0]
                new_seq_name = re.sub("^(_R_)*{}_*".format(fname), "", seq.name)
                all_taxa.add(new_seq_name)
                seq.id = new_seq_name
                seq.name = new_seq_name
                new_align.append(seq)
        if args.taxa is not None:
            assert len(all_taxa) == args.taxa, "Taxon names are not identical"
        outf = os.path.join(args.output, os.path.split(f)[1])
        try:
            AlignIO.write(new_align, open(outf, 'w'), 'nexus')
            print count
        except ValueError:
            raise IOError("Cannot write output file.")
    print "Taxon names in alignments: {0}".format(','.join(list(all_taxa)))

if __name__ == '__main__':
    main()
