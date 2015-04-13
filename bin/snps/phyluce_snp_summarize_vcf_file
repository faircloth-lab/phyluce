#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 29 August 2014 13:28 CDT (-0500)
"""


import os
import vcf
import argparse
from collections import defaultdict, Counter

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument(
            '--input',
            required=True,
            help="The vcf file to process"
        )
    return parser.parse_args()

def main():
    args = get_args()
    # setup counter for locus
    locus_count = []
    with open(args.input, 'r') as infile:
        vcf_reader = vcf.Reader(infile)
        for record in vcf_reader:
            # make sure that we get only passing (just in case)
            locus = record.CHROM.split("_")[0]
            locus_count.append(locus)
            if record.FILTER == []:
                complete = []
                for sample in record.samples:
                    # check for min completeness
                    if sample.data.GT == None:
                        complete.append(0)
                    else:
                        complete.append(1)
    c = Counter(locus_count)
    print "{} variable loci".format(len(c.keys()))
    print "{} variable SNPs per locus".format(float(sum(c.values()))/len(c.values()))
    print "{} min SNPs".format(min(c.values()))
    print "{} max SNPs".format(max(c.values()))

if __name__ == '__main__':
    main()
