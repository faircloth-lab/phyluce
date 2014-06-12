#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 12 June 2014 09:58 PDT (-0700)
"""


import re
import sys
import argparse
from phyluce.helpers import FullPaths, is_file


def get_args():
    parser = argparse.ArgumentParser(
            description="""Parse the log file from match_contigs_to_barcodes to output a table of results""")
    parser.add_argument(
            "--log",
            required=True,
            action=FullPaths,
            type=is_file,
            help="""The log file to parse"""
        )
    parser.add_argument(
            "--output",
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="""The output CSV file to create"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    # -------------- Processing xxx-yyy --------------
    regex1 = re.compile(".*Parsing\sFASTA\s(.*)$")
    regex2 = re.compile(".*\tBest\sBOLD\ssystems\smatch\sfor\slocus\s(.*):\s.*$")
    old_match = None
    args.output.write("[assemblies]\n")
    cnt = 0
    with open(args.log) as infile:
        for line in infile:
            match1 = regex1.search(line)
            match2 = regex2.search(line)
            if match1 is not None and old_match is None:
                old_match = match1.groups()[0]
            elif match1 is not None and match1 != old_match:
                old_match = match1.groups()[0]
            if match2:
                # add cnt to key to ensure unique
                args.output.write("{}|{}:{}\n".format(old_match, cnt, match2.groups()[0]))
                cnt += 1




if __name__ == '__main__':
    main()
