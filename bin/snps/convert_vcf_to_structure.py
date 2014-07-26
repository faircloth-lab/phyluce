#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2014 12:57 PDT (-0700)
"""


import sys
import vcf
import argparse
import pandas as pd
from collections import Counter

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument(
        '--input',
        required=True,
        help="The VCF file to process"
    )
    parser.add_argument(
        '--output',
        required=True,
        help="The STRUCTURE file to write"
    )
    subgroup = parser.add_mutually_exclusive_group()
    subgroup.add_argument(
        "--filter-identical",
        action="store_true",
        default=False,
        help="""Remove sites with only a single genotype.""",
    )
    subgroup.add_argument(
        "--filter-informative",
        action="store_true",
        default=False,
        help="""Remove uninformative sites.""",
    )
    return parser.parse_args()


def identical(items):
    return all(x == items[0] for x in items if x != "-9\t-9")


def informative(items):
    cnt = Counter([x for x in items if x != "-9\t-9"])
    #pdb.set_trace()
    if len(cnt) <= 1:
        return False
    elif len(cnt) >= 3:
        return True
    elif min(cnt.values()) >= 2:
        return True
    else:
        return False


def progress(cnt):
    if cnt % 100 == 0:
        sys.stdout.write(".")
        sys.stdout.flush()


def main():
    args = get_args()
    alleles = {"A":1, "T":2, "G":3, "C":4}
    sys.stdout.write("Working")
    sys.stdout.flush()
    with open(args.input, 'r') as infile:
        vcf_reader = vcf.Reader(infile)
        # create data frame with indiv as row id
        data = pd.DataFrame()
        for cnt, record in enumerate(vcf_reader):
            # prog
            progress(cnt)
            # get ref/alt alleles and make those into a list - integer values
            # in genotypes refer to 0-index positions of this list.
            reference = [record.REF] + [str(i) for i in record.ALT]
            # we'll keep a list of all genotypes
            all_genotypes = []
            for sample in record.samples:
                # deal with missing data
                if sample.data.GT == None:
                    gt = "\t".join(["-9","-9"])
                    # bung those into list
                    all_genotypes.append(gt)
                else:
                    # split the integer code
                    coded_gt = [int(i) for i in sample.data.GT.split("/")]
                    # get the ref/alt allele based on integer code and lookup
                    # new integer value in dict/hash table
                    gt = sorted([alleles[reference[i]] for i in coded_gt])
                    # join those to a string
                    gt = "\t".join([str(i) for i in gt])
                    # bung those into list
                    all_genotypes.append(gt)
            # filter monomorphic loci
            if (args.filter_identical and not identical(all_genotypes)):
                temp_df = pd.DataFrame({cnt:pd.Series(all_genotypes,index=vcf_reader.samples)})
                data = pd.concat([data, temp_df], axis=1)
            # filter monomorphic loci
            elif (args.filter_informative and informative(all_genotypes)):
                temp_df = pd.DataFrame({cnt:pd.Series(all_genotypes,index=vcf_reader.samples)})
                data = pd.concat([data, temp_df], axis=1)
            elif (not args.filter_identical and not args.filter_informative):
                temp_df = pd.DataFrame({cnt:pd.Series(all_genotypes,index=vcf_reader.samples)})
                data = pd.concat([data, temp_df], axis=1)
    with open(args.output, 'w') as outf:
        loci_header = "\t".join(["{0}\t{0}".format(i) for i in data.columns.values.tolist()])
        outf.write("\t{}\n".format(loci_header))
        for row in data.iterrows():
            #pdb.set_trace()
            outf.write("{0}\t{1}\n".format(
                    row[0],
                    "\t".join(row[1])
                ))
    # flush
    print ""
    print "Data contain {} individuals and {} loci output to {}.".format(
        len(data),
        len(data.columns.values.tolist()),
        args.output
    )


if __name__ == '__main__':
    main()
