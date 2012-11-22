#!/usr/bin/env python
# encoding: utf-8
"""
File: concatenate_raxml_bootreps.py
Author: Brant Faircloth

Created by Brant Faircloth on 30 September 2012 10:09 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import glob
import argparse
from phyluce.helpers import is_dir
from collections import defaultdict

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Given an input directory containing RAxML bootreps, concatenate those into a file""")
    parser.add_argument(
            "input",
            type=is_dir,
            help="""The input directory, containing locus-specific output directories from RAxML"""
        )
    parser.add_argument(
            "output",
            help="""The output file in which to store the bootreps."""
        )
    return parser.parse_args()


def main():
    args = get_args()
    all_bootreps = defaultdict(lambda: defaultdict(str))
    for root, dirs, files in os.walk(args.input):
        for d in dirs:
            files = glob.glob(os.path.join(root, os.path.join(d, '*.ml')))
            bootreps = [f for f in files if "RAxML_bootstrap" in os.path.basename(f)]
            assert len(bootreps) == 1, "There appear to be >1 bootstrap files"
            bootrep = bootreps[0]
            with open(bootrep, 'rU') as file:
                for lineno, line in enumerate(file):
                    all_bootreps[lineno][d] = line.strip()
    #pdb.set_trace()
    outp = open(args.output, 'w')
    loci = sorted(all_bootreps[0].keys())
    for bootrep, trees in all_bootreps.iteritems():
        for locus in loci:
            outp.write('''{0}\t"{1}"\n'''.format(bootrep, trees[locus]))
    outp.close()


if __name__ == '__main__':
    main()
