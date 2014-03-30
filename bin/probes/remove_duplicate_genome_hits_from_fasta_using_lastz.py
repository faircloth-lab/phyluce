#!/usr/bin/env python
# encoding: utf-8
"""
File: remove_duplicate_genome_hits_from_fasta_using_lastz.py
Author: Brant Faircloth

Created by Brant Faircloth on 23 September 2013 15:09 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import re
import argparse
from collections import defaultdict

from Bio import SeqIO

from phyluce import lastz
from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Remove duplicate loci from a fasta given self-to-self lastz"""
    )
    parser.add_argument(
        "--fasta",
        required=True,
        action=FullPaths,
        help="""The fasta file to screen"""
    )
    parser.add_argument(
        "--lastz",
        required=True,
        action=FullPaths,
        help="""The lastz file to use"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        help="""The output fasta to create"""
    )
    parser.add_argument(
        "--long",
        action="store_true",
        default=False,
        help="""If the lastz output is the longfield format""",
    )
    return parser.parse_args()


def get_dupes(lastz_file, format):
    """Given a lastz_file of probes aligned to themselves, get duplicates"""
    matches = defaultdict(list)
    dupes = set()
    # get names and strip probe designation since loci are the same
    print "Parsing lastz file..."
    for lz in lastz.Reader(lastz_file, long_format=format):
        target_name = "{}:{}-{}".format(lz.name1, lz.zstart1, lz.end1)
        query_name = lz.name2
        matches[query_name].append(target_name)
    # see if one probe matches any other probes
    # other than the children of the locus
    print "Screening results..."
    for k, v in matches.iteritems():
        name, pos = k.split('|')
        name = name.strip()
        # if the probe doesn't match itself, we have
        # problems
        if len(v) == 1:
            assert pos == v[0], "no match"
        elif v > 1:
            dupes.add(name)
    # make sure all names are lowercase
    return set([d.lower() for d in dupes])


def main():
    args = get_args()
    dupes = get_dupes(args.lastz, args.long)
    kept = 0
    with open(args.fasta, 'rU') as infile:
        with open(args.output, 'w') as outf:
            for cnt, seq in enumerate(SeqIO.parse(infile, 'fasta')):
                if not seq.id in dupes:
                    outf.write(seq.format('fasta'))
                    kept += 1
                else:
                    pass
    print "Screened {} sequences.  Filtered {} duplicates. Kept {}.".format(
        cnt + 1,
        len(dupes),
        kept
    )


if __name__ == '__main__':
    main()
