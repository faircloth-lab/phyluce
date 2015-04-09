#!/usr/bin/env python
# encoding: utf-8
"""
File: new_gt_bootrep_sorter.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 October 2013 17:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import glob
import gzip
import argparse
import tempfile
import subprocess

from itertools import groupby
from operator import itemgetter
from phyluce.helpers import FullPaths, CreateDir, is_dir


import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Program description"""
    )
    parser.add_argument(
        "--input",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="""The path to the directory containing bootstrap replicates."""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The path to the output file to create."""
    )
    return parser.parse_args()

args = get_args()
handle1, pth1 = tempfile.mkstemp(suffix=".tre")
for file in glob.glob(os.path.join(args.input, "*.bootrep")):
    name = os.path.basename(file)
    with open(file, 'rU') as infile:
        for cnt, line in enumerate(infile):
            os.write(handle1, "{}\t{}\t{}".format(
                name,
                cnt,
                line
            ))
os.close(handle1)
handle2, pth2 = tempfile.mkstemp(suffix=".sorted.tre")
os.close(handle2)
# sort the file rapidly
print "sorting..."
cmd = [
    "sort",
    "-k",
    "2,2",
    "-k",
    "1n",
    "-o",
    pth2,
    pth1
]
proc = subprocess.Popen(cmd)
proc.communicate()
# remove the unsorted file
os.remove(pth1)
print "removed {}".format(pth1)
# split the sorted file into subfiles
split = (line.split("\t") for line in open(pth2))
for key, rows in groupby(split, itemgetter(1)):
    with gzip.open(os.path.join(args.output, "{}.tre.gz".format(key)), "w") as output:
        for row in rows:
            output.write("\t".join(row))
# remove the sorted file
os.remove(pth2)
print "removed {}".format(pth2)

