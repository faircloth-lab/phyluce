#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 15:02 PDT (-0700)
"""

import sys
import argparse
import subprocess
from phyluce.helpers import is_file

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given an input BED file and BAM file, output an interval.list file"""
    )
    parser.add_argument(
        "--bed",
        required=True,
        type=is_file,
        help="""The BED file to filter."""
    )
    parser.add_argument(
        "--bam",
        required=True,
        type=is_file,
        help="""The BED file to filter."""
    )
    parser.add_argument(
        "--output",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="""The output BED file"""
    )
    return parser.parse_args()


def get_bam_header(bam):
    cmd = ["samtools", "view", "-H", bam]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    return stdout


def main():
    args = get_args()
    header = get_bam_header(args.bam)
    args.output.write(header)
    with open(args.bed, "rU") as infile:
        for line in infile:
            if line.startswith("track"):
                pass
            else:
                ls = line.strip().split("\t")
                args.output.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(ls[0],ls[1],ls[2],ls[5],ls[3])
                )
    args.output.close()

if __name__ == '__main__':
    main()
