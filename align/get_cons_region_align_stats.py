#!/usr/bin/env python
# encoding: utf-8

"""
get_cons_region_align_stats.py

Created by Brant Faircloth on 10 July 2010 10:09 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  this program iterates through a folder of aligned reads (either
trimmed or untrimmed) and tallies the count of bases at which each species
differs.

"""

import pdb
import os
import sys
import glob
import numpy as np
import optparse
from Bio import AlignIO


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the input directory.', 
metavar='FILE')

    p.add_option('--output', dest = 'output', action='store', 
type='string', default = None, help='The path to the output file.', 
metavar='FILE')

#    p.add_option('--chicken', dest = 'chicken', action='store_true', default=False, 
#help='Include chicken sequence reads')

    (options,arg) = p.parse_args()
    if not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isdir(os.path.expanduser(options.input)):
        print "You must provide a valid path to an input directory."
        p.print_help()
        sys.exit(2)
    return options, arg 


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    options, args = interface()
    # iterate through all the files to determine the longest alignment
    files = get_files(options.input)
    align_lengths = [AlignIO.read(f, 'nexus').get_alignment_length() \
                            for f in files]
    max_align_length = max(align_lengths)
    # find the middle
    middle = int(round(max_align_length/2, 0))
    # create a dict to hold the results by position in longest array
    differences = dict((d,np.array([])) for d in range(-middle, middle + 1))
    # iterate through all the files again
    for f in files:
        align = AlignIO.read(f, 'nexus')
        align_length = align.get_alignment_length()
        align_diff = int((round((max_align_length - align_length)/2.,0) - middle))
        # determine relative start of this alignment to longest
        for col in xrange(align_length):
            bases = align.get_column(col)
            b_counts = len(set(bases))
            #pdb.set_trace()
            differences[align_diff + col] = np.append(differences[align_diff + col],b_counts)
    #pdb.set_trace()
    position = differences.keys()
    # create bins in groups of 100
    #bins = np.array(range(0,500,100))
    #print bins
    #pdb.set_trace()
    if options.output:
        outp = open(options.output, 'w')
    else:
        outp = sys.stdout
    outp.write('bp, mean, ci, onediff, greaterthanonediff, fourdiff, count\n')
    ignore = []
    for p in sorted(position):
        #pdb.set_trace()
        # how many only have 0-1 difference
        try:
            one_diff = sum(differences[p] <= 1)/float(len(differences[p]))
            four_diff = sum(differences[p] >= 4)/float(len(differences[p]))
            greater_than_one_diff = sum(differences[p] > 1)/float(len(differences[p]))
            total = len(differences[p])
        except ZeroDivisionError:
            ignore.append(p)
        if p not in ignore:
            outp.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(p, np.mean(differences[p]), 1.96 * np.std(differences[p]), one_diff, four_diff, greater_than_one_diff, total))
    if options.output:
        outp.close()

if __name__ == '__main__':
    main()