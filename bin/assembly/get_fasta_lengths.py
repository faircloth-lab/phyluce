#!/usr/bin/env python
# encoding: utf-8
"""
File: get_fastq_lengths.py
Author: Brant Faircloth

Created by Brant Faircloth on 21 July 2012 12:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""
import os
import gzip
import glob
import numpy
import tempfile
import argparse
import subprocess
from itertools import groupby, imap
from phyluce.helpers import is_dir, FullPaths
from Bio import SeqIO

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Get summary (length) data from fastq""")
    parser.add_argument(
            "input",
            action=FullPaths,
            help="""The directory of fastq files to summarize"""
        )
    parser.add_argument(
            "--csv",
            action="store_true",
            default=False,
            help="""Give output in CSV"""
        )
    return parser.parse_args()

def fasta_iter(fasta):
    """modified from @brent_p on stackoverflow.  yield tuple of header, sequence"""
    if fasta.endswith('.gz'):
        with gzip.open(fasta) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                yield sum(len(s.strip()) for s in faiter.next())
    else:
        with open(fasta) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                yield sum(len(s.strip()) for s in faiter.next())


def main():
    args = get_args()
    lengths = numpy.array([int(record) for record in fasta_iter(args.input)])
    std_error = numpy.std(lengths, ddof=1) / numpy.sqrt(len(lengths))
    if not args.csv:
        print "Reads:\t\t{:,}".format(len(lengths))
        print "Bp:\t\t{:,}".format(sum(lengths))
        print "Avg. len:\t{:,}".format(numpy.average(lengths))
        print "STDERR len:\t{:,}".format(std_error)
        print "Min. len:\t{:,}".format(min(lengths))
        print "Max. len:\t{:,}".format(max(lengths))
        print "Median len:\t{:,}".format(numpy.median(lengths))
        print "Contigs > 1kb:\t{:,}".format(sum(lengths >= 1000))
    else:
        try:
            print "{},{},{},{},{},{},{},{},{}".format(os.path.basename(args.input), len(lengths), sum(lengths), numpy.average(lengths), std_error, min(lengths), max(lengths), numpy.median(lengths), sum(lengths >= 1000))
        except:
            print "{},{},{},{},{},{},{}".format(os.path.basename(args.input), len(lengths), "Div/0", "Div/0", "Div/0","Div/0", sum(lengths >= 1000))


if __name__ == '__main__':
    main()
