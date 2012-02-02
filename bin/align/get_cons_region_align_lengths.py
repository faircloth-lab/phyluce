#!/usr/bin/env python
# encoding: utf-8

"""
get_cons_region_lengths.py

Created by Brant Faircloth on 30 July 2010 10:09 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  this program iterates through a folder of aligned reads (either
trimmed or untrimmed) and tallies the count of bases at which each species
differs.

"""

import pdb
import os
import sys
import glob
import numpy
import optparse
from Bio import AlignIO


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the input directory.', 
metavar='FILE')

    p.add_option('--verbose', dest = 'verbose', action='store_true', \
default = False, help='Verbose output')

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
    files = glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))
    files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nexus')))
    return files

def main():
    options, args = interface()
    # iterate through all the files to determine the longest alignment
    files = get_files(options.input)
    align_lengths = [[AlignIO.read(f, 'nexus').get_alignment_length(),os.path.split(f)[1]] \
                            for f in files]
    #pdb.set_trace()
    l = numpy.array([l[0] for l in align_lengths])
    #pdb.set_trace()
    print "Average length: ", numpy.mean(l)
    print "95 CI length: ", 1.96 * (numpy.std(l, ddof=1)/numpy.sqrt(len(l)))
    if options.verbose:
        for a in align_lengths:
            print a[1].strip('.nex'), a[0]
    #pdb.set_trace()

if __name__ == '__main__':
    main()
