#!/usr/bin/env python
# encoding: utf-8

"""
get_dbsnp_positions_along_uces.py

Created by Brant Faircloth on 02 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import numpy
import argparse
from phyluce.helpers import get_dupes, get_xml_data

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get the average contig length, by species, within a combined fasta')
    parser.add_argument('dbsnp', help='CSV input from dbSNP giving SNP positions within UCE', type=argparse.FileType('rU'))
    parser.add_argument('xml', help='The XML file holiding locus data from dbSNP')
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--dupefile', help='The path to a lastz file of lastz-against-self results')
    return parser.parse_args()

def main():
    args = get_args()
    if args.dupefile:
        dupes = get_dupes(args.dupefile)
    else:
        dupes = None
    #pdb.set_trace()
    # get dbSNP data
    all_snps = get_xml_data(args.xml)
    used = set()
    # iterate over intersections
    args.output.write('rsid,pos,maf,1000g\n')
    for row in args.dbsnp:
        if not row.startswith('UCE'):
            uce, chromo, start, end, snp, snps, snpe = row.strip('\n').split(',')
            start, end, snps, snpe = map(int, [start, end, snps, snpe])
            # get relative position
            if not snpe - snps > 1 and snp not in used and not uce in dupes:
                middle = int(round((start + end)/2, 0))
                rel_snp_pos = snps - middle
                # lookup data for snps
                if all_snps[snp.strip('rs')].val_1000G and all_snps[snp.strip('rs')].val_1000G.lower() == 'true':
                    thousandg = True
                else:
                    thousandg = False
                if not all_snps[snp.strip('rs')].freq_freq:
                    freq = 0.0
                else:
                    freq = float(all_snps[snp.strip('rs')].freq_freq)
                args.output.write("{0},{1},{2},{3}\n".format(
                    snp, 
                    rel_snp_pos,
                    freq, 
                    thousandg
                    )
                )
                # make sure we skip any duplicates
                used.add(snp)
                #pdb.set_trace()
        
    # check to see if UCE is in dupes

if __name__ == '__main__':
    main()
