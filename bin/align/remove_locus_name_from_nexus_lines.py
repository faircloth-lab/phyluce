#!/usr/bin/env python
# encoding: utf-8

"""
change_taxa_names.py

Created by Brant Faircloth on 22 September 2010 12:48 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  Add the filename to the taxa name in each of the output nexus files - for use
with *BEAST.

USAGE:  python remove_locus_name_from_nexus_lines.py \
    --input my/input/folder/nexus \
    --output my/input/folder/nexus-renamed
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
    p.add_option('--position', dest = 'position', action='store', 
type='int', default = -2, help='The position of the name to keep')

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
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex*'))

def main():
    options, args = interface()
    # iterate through all the files to determine the longest alignment
    files = get_files(options.input)
    for count, f in enumerate(files):
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        #filename = os.path.basename(f)
        #chromo_name = filename.split('.')[0]
        for align in AlignIO.parse(f, 'nexus'):
            for seq in list(align):
                if '.copy' in seq.name:
                    pass
                else:
                #pdb.set_trace()
                #new_seq_name = seq.name.split('|')[0]
                    new_seq_name = '_'.join(seq.name.split('_')[options.position:])
                    new_align.add_sequence(new_seq_name, str(seq.seq))
        #pdb.set_trace()
        outf = os.path.join(options.output, os.path.split(f)[1])
        try:
            AlignIO.write(new_align, open(outf, 'w'), 'nexus')
        except ValueError:
            pdb.set_trace()
        print count

if __name__ == '__main__':
    main()
