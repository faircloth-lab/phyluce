#!/usr/bin/env python
# encoding: utf-8

"""
change_taxa_names.py

Created by Brant Faircloth on 22 September 2010 12:48 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  Add the filename to the taxa name in each of the output nexus files - for use
with *BEAST.

USAGE:  python change_taxa_names.py \
    --input=Bats_177_Loci_12Species_alignments_nexus \
    --output=Bats_177_Loci_12Species_alignments_nexus_renamed
"""

import pdb
import os
import sys
import glob
import optparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment

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
    p.add_option('--shorten-name', dest = 'shorten_name', action='store_true', default=False, \
help='Shorten long names.')
    p.add_option('--positions', dest = 'positions', action='store', default='0,1')

    (options,arg) = p.parse_args()
    options.input = os.path.abspath(os.path.expanduser(options.input))
    options.output = os.path.abspath(os.path.expanduser(options.output))
    if not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isdir(options.input) or not os.path.isdir(options.output):
        print "You must provide a valid path to an input directory."
        p.print_help()
        sys.exit(2)
    return options, arg

def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def rename(align, first, second):
    for a in align:
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        for seq in a:
            split_name = seq.id.split('_')
            #pdb.set_trace()
            if first and second:
                new_seq_name = '_'.join([split_name[first][0:3], split_name[second][0:3]])
            elif not second:
                new_seq_name = split_name[first]
            new_align.add_sequence(new_seq_name, str(seq.seq))
        yield new_align

def main():
    options, args = interface()
    # iterate through all the files to determine the longest alignment
    files = get_files(options.input)
    if options.positions:
        pos1, pos2 = [eval(i) for i in options.positions.strip().split(',')]
    for count, f in enumerate(files):
        align = AlignIO.parse(f, "nexus")
        #pdb.set_trace()
        new_name = os.path.splitext(os.path.split(f)[1])[0] + '.phylip'
        #pdb.set_trace()
        outf = os.path.join(options.output, new_name)
        #try:
        if options.shorten_name:
            for item in rename(align, pos1, pos2):
                AlignIO.write(item, open(outf, 'w'), 'phylip')
        else:
            AlignIO.write(align, open(outf, 'w'), 'phylip')
        #except ValueError:
            #pdb.set_trace()
            
        print count

if __name__ == '__main__':
    main()