#!/usr/bin/env python
# encoding: utf-8
"""
File: parse_trinity_coverage_for_uce_loci_log.py
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
    """Parse the log file from get_trinity_coverage_for_uce_loci.py to output a nice table of results"""
    parser = argparse.ArgumentParser(
            description="""Parse the log file from get_trinity_coverage_for_uce_loci.py to output a nice table of results""")
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
    args = get_args()
    # -------------- Processing xxx-yyy --------------
    regex1 = re.compile(".*-+\sProcessing\s(.*)\s-+$")
    #    xxx contigs, mean trimmed length = xx.x, mean trimmed coverage = xx.x, on-target bases (uce contigs) = xx.x%, unique reads aligned (all contigs) = xx.x%
    rcontigs = "\t(\d+)\scontigs,"
    rlength = "\smean\strimmed\slength\s=\s(\d+.\d+),"
    rcoverage = "\smean\strimmed\scoverage\s=\s(\d+.\d+)x,"
    rontarget = "\son-target\sbases\s\(uce\scontigs\)\s=\s(\d+.\d+%),"
    runique = "\sunique\sreads\saligned\s\(all\scontigs\)\s=\s(\d+.\d+%)"
    regex2 = re.compile(".*{}{}{}{}{}$".format(
        rcontigs,
        rlength,
        rcoverage,
        rontarget,
        runique
    ))
    d = defaultdict(str)
    args.output.write("taxon,UCE contigs,UCE contigs mean length,UCE contigs coverage (x),UCE contigs reads on target,UCE contigs unique reads aligned,\n")
    with open(args.log) as infile:
        for line in infile:
            match1 = regex1.search(line)
            match2 = regex2.search(line)
            if match1:
                taxon = match1.groups()[0]
            if match2:
                contigs, length, coverage, on_target, unique = match2.groups()
                s = "{},{},{},{},{}".format(
                    contigs,
                    length,
                    coverage,
                    on_target,
                    unique,
                )
                d[taxon] = s
                # reset variables to None
                contigs, length, coverage, on_target, unique = None, None, None, None, None
    for taxon in sorted(d.keys()):
        args.output.write("{},{}\n".format(taxon, d[taxon]))

if __name__ == '__main__':
    main()


