#!/usr/bin/env python
# encoding: utf-8

"""
get_dbsnp_freq_stats.py

Created by Brant Faircloth on 02 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import copy
import numpy
import argparse
from phyluce.helpers import get_dupes, get_xml_data

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get the average contig length, by species, within a combined fasta')
    parser.add_argument('dbsnp', help='CSV input from dbSNP giving SNP positions within UCE')
    parser.add_argument('xml', help='The XML file holiding locus data from dbSNP')
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--output2', help = 'The output file',
        type=argparse.FileType('w'))
    parser.add_argument('--dupefile', help='The path to a lastz file of lastz-against-self results')
    return parser.parse_args()

def main():
    args = get_args()
    if args.dupefile:
        dupes = get_dupes(args.dupefile)
    else:
        dupes = None
    used = set()
    mx = max([int(row.strip('\n').split(',')[3]) \
            - int(row.strip('\n').split(',')[2]) \
            for row in open(args.dbsnp,'rU') if not row.startswith('UCE')])
    # get the SNP metadata
    all_snps = get_xml_data(args.xml)
    # find the middle
    overall_middle = int(round(mx/2, 0))
    # list to hold results 
    l = numpy.zeros(mx + 1)
    positions = copy.deepcopy(l)
    # create a dict to hold the results by position in longest array
    #differences = dict((d,numpy.array([])) for d in range(-middle, middle + 1))
    # iterate over intersections
    d = {}
    if args.output2:
        args.output2.write('UCE,chromo,uce-start,uce-end,snp-name,snp-start,snp-end,1000gvalidated,freq\n')
    for row in open(args.dbsnp, 'rU'):
        if not row.startswith('UCE'):
            uce, chromo, start, end, snp, snps, snpe = row.strip('\n').split(',')
            start, end, snps, snpe = map(int, [start, end, snps, snpe])
            # get middle of this UCE
            middle = int(round((start + end)/2, 0))
            #pdb.set_trace()
            if snp not in used:
                if not snpe - snps > 1 \
                    and (uce not in dupes) \
                    and all_snps[snp.strip('rs')].val_1000G == 'true' \
                    and all_snps[snp.strip('rs')].freq_freq is not None:
                    if not uce in d.keys():
                        d[uce] = numpy.zeros(mx + 1)
                    rel_snp_pos = snps - middle
                    d[uce][overall_middle + rel_snp_pos] = all_snps[snp.strip('rs')].freq_freq
                if args.output2 and not snpe - snps > 1 and (snp not in used) and (uce not in dupes):
                    args.output2.write("{},{},{},{},{},{},{},{},{}\n".format(
                        uce, chromo, start, end, snp, snps, 
                        snpe, all_snps[snp.strip('rs')].val_1000G, 
                        all_snps[snp.strip('rs')].freq_freq))
                used.add(snp)
    stack = numpy.array([d[uce] for uce in d.keys()])
    #pdb.set_trace()
    # compute the running average
    win = 25
    data = sum(stack > 0)
    weightings = numpy.repeat(1.0, win) / win
    running = numpy.convolve(data, weightings)[win-1:-(win-1)]
    args.output.write("pos,avg,ci,datatype\n")
    for base in range(len(running)):
        pos = base - overall_middle
        args.output.write("{},{},,running\n".format(pos,running[base]))
    # also output the average heterozygosity of 1000 Genome validated, hetero SNPs.
    for base in range(len(stack[0])):
        pos = base - overall_middle
        values = numpy.where(stack[:,base] != 0)[0]
        # reindex
        avg = numpy.mean(stack[:,base][values])
        ci = 1.96 * (numpy.std(stack[:,base][values], ddof = 1)/numpy.sqrt(len(stack[:,base][values])))
        args.output.write("{},{},{},mean_hetero\n".format(pos, avg, ci))
    win = 25
    data = numpy.mean(stack, axis = 1)
    weightings = numpy.repeat(1.0, win) / win
    running = numpy.convolve(data, weightings)[win-1:-(win-1)]
    for base in range(len(stack[0])):
        pos = base - overall_middle
        args.output.write("{},{},,running_hetero\n".format(pos,running[base]))

if __name__ == '__main__':
    main()
