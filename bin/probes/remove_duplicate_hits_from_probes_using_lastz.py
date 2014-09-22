#!/usr/bin/env python
# encoding: utf-8
"""
File: remove_duplicate_hits_from_fasta_using_lastz.py
Author: Brant Faircloth

Created by Brant Faircloth on 09 June 2013 12:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
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
        '--probe-prefix',
        required=True,
        type=str,
        default=None,
        help='The prefix (e.g. "uce-") added to all probes designed'
    )
    parser.add_argument(
        '--probe-regex',
        type=str,
        default='^({}\d+)(?:_p\d+.*)',
        help='The regular expression to use for matching probes'
    )
    parser.add_argument(
        '--probe-bed',
        type=str,
        default=None,
        help='The path to a file contaning the probe coordinates in BED format'
    )
    parser.add_argument(
        '--locus-bed',
        type=str,
        default=None,
        help='The path to a file contaning the locus coordinates in BED format'
    )
    parser.add_argument(
        "--long",
        action="store_true",
        default=False,
        help="""If the lastz output is the longfield format""",
    )
    return parser.parse_args()


def new_get_probe_name(header, regex):
    match = re.search(regex, header)
    return match.groups()[0]


def get_dupes(lastz_file, regex, format):
    """Given a lastz_file of probes aligned to themselves, get duplicates"""
    matches = defaultdict(list)
    dupes = set()
    # get names and strip probe designation since loci are the same
    print "Parsing lastz file..."
    for lz in lastz.Reader(lastz_file, long_format=format):
        target_name = new_get_probe_name(lz.name1, regex)
        query_name = new_get_probe_name(lz.name2, regex)
        matches[target_name].append(query_name)
    # see if one probe matches any other probes
    # other than the children of the locus
    print "Screening results..."
    for k, v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    # make sure all names are lowercase
    return set([d.lower() for d in dupes])


def return_screened_filename(s):
    pth, old_name = os.path.split(s)
    old_name, ext = os.path.splitext(old_name)
    new_name = "{}-DUPE-SCREENED{}".format(
        old_name,
        ext
    )
    return os.path.join(pth, new_name)


def main():
    args = get_args()
    regex = args.probe_regex.format(args.probe_prefix)
    dupes = get_dupes(args.lastz, regex, args.long)
    fasta_output = return_screened_filename(args.fasta)
    fasta_kept = 0
    with open(args.fasta, 'rU') as infile:
        with open(fasta_output, 'w') as outf:
            for fasta_cnt, seq in enumerate(SeqIO.parse(infile, 'fasta')):
                match = re.search(regex, seq.id)
                locus = match.groups()[0]
                if not locus in dupes:
                    outf.write(seq.format('fasta'))
                    fasta_kept += 1
                else:
                    pass
    print "Screened {} fasta sequences.  Filtered {} duplicates. Kept {}.".format(
        fasta_cnt,
        len(dupes),
        fasta_kept
    )
    if args.probe_bed:
        probe_kept = 0
        probe_bed_output = return_screened_filename(args.probe_bed)
        with open(args.probe_bed, 'rU') as infile:
            with open(probe_bed_output, 'w') as outf:
                for probe_cnt, line in enumerate(infile):
                    if line.startswith("track"):
                        header = True
                        outf.write(line)
                    else:
                        ls = line.strip().split("\t")
                        match = re.search(regex, ls[3])
                        locus = match.groups()[0]
                        if not locus in dupes:
                            outf.write(line)
                            probe_kept += 1
                        else:
                            pass
        if header:
            probe_cnt -= 1
        print "Screened {} BED probe sequences.  Filtered {} duplicates. Kept {}.".format(
            probe_cnt,
            len(dupes),
            probe_kept
        )
    if args.locus_bed:
        locus_kept = 0
        locus_bed_output = return_screened_filename(args.locus_bed)
        with open(args.locus_bed, 'rU') as infile:
            with open(locus_bed_output, 'w') as outf:
                for locus_cnt, line in enumerate(infile):
                    if line.startswith("track"):
                        header = True
                        outf.write(line)
                    else:
                        ls = line.strip().split("\t")
                        locus = ls[3]
                        if not locus in dupes:
                            outf.write(line)
                            locus_kept += 1
                        else:
                            pass
        if header:
            locus_cnt -= 1
        print "Screened {} BED locus sequences.  Filtered {} duplicates. Kept {}.".format(
            locus_cnt,
            len(dupes),
            locus_kept
        )




if __name__ == '__main__':
    main()
