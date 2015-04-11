#!/usr/bin/env python
# encoding: utf-8
"""
File: check_raxml_zcluster_runs.py
Author: Brant Faircloth

Created by Brant Faircloth on 19 March 2013 11:03 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import glob
import argparse
from phyluce.helpers import FullPaths, is_dir

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Check zcluster results for completion""")
    parser.add_argument(
            "alignments",
            action=FullPaths,
            type=is_dir,
            help="""The directory holding the alignments"""
        )
    parser.add_argument(
            "genetrees",
            action=FullPaths,
            type=is_dir,
            help="""The directory holding the genetrees"""
        )
    parser.add_argument(
            "bootreps",
            action=FullPaths,
            type=is_dir,
            help="""The directory holding the boostrap replications"""
        )
    parser.add_argument(
            "--replicates",
            type=int,
            default=200,
            help="""The number of bootreps to check for""",
        )
    return parser.parse_args()


def get_locus_names(alignments):
    pth = os.path.join(alignments , "*.phylip")
    loci = []
    for align in glob.glob(pth):
        loci.append(os.path.splitext(os.path.basename(align))[0])
    print "There are {} alignment files".format(len(loci))
    return loci


def check_genetree_files(genetrees, loci):
    """docstring for check_genetree_info_files"""
    regex = re.compile("Best-scoring\sML\stree\swritten\sto:\s(.*)")
    for locus in loci:
        filename = "RAxML_info.{}.best".format(locus)
        pth = os.path.join(genetrees, filename)
        data = open(pth, 'rU').read()
        result = regex.search(data)
        treefile = os.path.basename(result.groups()[0])
        assert os.path.exists(os.path.join(genetrees, treefile)), "File {} does not exist".format(treefile)
    print "All {} genetree files exist".format(len(loci))


def get_file_len(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def check_bootrep_files(bootreps, replicates, loci):
    """docstring for check_bootrep_files"""
    regex = re.compile("All\s200\sbootstrapped\strees\swritten\sto:\s(.*)")
    for locus in loci:
        filename = "RAxML_info.{}.bootrep".format(locus)
        pth = os.path.join(bootreps, filename)
        data = open(pth, 'rU').read()
        result = regex.search(data)
        try:
            bootrepfile = os.path.basename(result.groups()[0])
            lines = get_file_len(os.path.join(bootreps, bootrepfile))
            assert lines == replicates, "Bootrep file {} only has {} lines".format(bootrepfile, lines)
        except:
            print "{} does not have bootstrapped trees".format(filename)
        
def main():
    """docstring for main"""
    args = get_args()
    loci = get_locus_names(args.alignments)
    check_genetree_files(args.genetrees, loci)
    check_bootrep_files(args.bootreps, args.replicates, loci)

if __name__ == '__main__':
    main()