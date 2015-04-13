#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 14 April 2014 12:57 PDT (-0700)
"""

import re
import sys
import vcf
import argparse
import pandas as pd
from collections import Counter
from Bio import AlignIO
from Bio.Alphabet import generic_dna, Gapped
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Nexus import Nexus

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

'''
class SnappNumericAlphabet(Alphabet):
    """Generic alphabet with letters of size one."""
    size = 1
    letters = "012"   # string of all letters in the alphabet

class DNAAlphabet(SnappNumericAlphabet):
    """Generic alphabet with letters of size one."""
    pass
'''

def identical(items):
    return all(x == items[0] for x in items if x != "-9\t-9")


def informative(items):
    cnt = Counter([x for x in items if x != "?"])
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


def geno_match(gt):
    # homozygote
    if re.match("^0/0$", gt):
        return "0"
    # heterozygote, may have mult. non-ref alleles as alt. starting w/ "1"
    elif re.match("^0/[1-9]+$", gt):
        return "1"
    # reverse of above, to be safe
    elif re.match("^[1-9]+/0$", gt):
        return "1"
    # homozygote, alt alleles
    elif re.match("^[1-9]+/[1-9]+$", gt):
        return "2"
    else:
        return "-"


def main():
    args = get_args()
    sys.stdout.write("Working")
    sys.stdout.flush()
    with open(args.input, 'r') as infile:
        vcf_reader = vcf.Reader(infile)
        # create data frame with indiv as row id
        data = pd.DataFrame()
        for cnt, record in enumerate(vcf_reader):
            # prog
            progress(cnt)
            # we'll keep a list of all genotypes
            all_genotypes = []
            for sample in record.samples:
                # deal with missing data
                if sample.data.GT == None:
                    #gt = "\t".join(["-9","-9"])
                    # bung those into list
                    all_genotypes.append("?")
                else:
                    gt = geno_match(sample.data.GT)
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
    # setup nexus file - cheat here by making this look like a DNA alphabet
    minimal_record = '''#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=dna missing=? gap=-; end;'''
    nex = Nexus.Nexus(minimal_record)
    nex.alphabet = generic_dna
    for row in data.iterrows():
        nex.add_sequence(row[0].replace("-", "_"), "".join(row[1]))
    nex.write_nexus_data(args.output, interleave=False)
    # flush
    print ""
    print "Data contain {} individuals and {} loci output to {}.".format(
        len(data),
        len(data.columns.values.tolist()),
        args.output
    )


if __name__ == '__main__':
    main()
