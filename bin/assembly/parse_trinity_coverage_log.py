#!/usr/bin/env python
# encoding: utf-8
"""
File: parse_trinity_coverage_log.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 October 2013 10:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import re
import sys
import argparse
from collections import defaultdict
from phyluce.helpers import FullPaths, is_file

def get_args():
    """Parse the log file from get_trinity_coverage.py to output a nice table of results"""
    parser = argparse.ArgumentParser(
            description="""Parse the log file from get_trinity_coverage.py to output a nice table of results""")
    parser.add_argument(
            "--log",
            required=True,
            action=FullPaths,
            type=is_file,
            help="""The log file to parse"""
        )
    parser.add_argument(
            "--output",
            action=FullPaths,
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="""The output CSV file to create"""
        )
    return parser.parse_args()


def main():
    # -------------- Processing xxx-yyy --------------
    regex1 = re.compile(".*-+\sProcessing\s(.*)\s-+$")
    # xxx contigs, mean coverage = yy.y, mean length = zzz.z
    regex2 = re.compile(".*\s\s\s\s(\d+)\scontigs,\smean\scoverage\s=\s(\d+.\d+),\smean\slength\s=\s(\d+.\d+)$")
    args = get_args()
    args.output.write("taxon,Total contigs (after trimming),Total contigs coverage (x),Total contigs mean length\n")
    d = defaultdict(str)
    with open(args.log) as infile:
        for line in infile:
            match1 = regex1.search(line)
            match2 = regex2.search(line)
            if match1:
                taxon = match1.groups()[0]
            if match2:
                contigs, coverage, mean_length = match2.groups()
                s = "{},{},{}".format(
                    contigs,
                    coverage,
                    mean_length
                )
                d[taxon] = s
                # reset variables to None
                taxon, contigs, coverage, mean_length = None, None, None, None
    for taxon in sorted(d.keys()):
        args.output.write("{},{}\n".format(taxon, d[taxon]))

if __name__ == '__main__':
    main()


