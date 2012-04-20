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
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment as Alignment

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
    p.add_option('--format', dest='format', action='store', 
type='string', default = 'phylip', help='The output format',)
    p.add_option('--splitchar', dest="splitchar", action='store', 
type='string', default = '_', help='The character to put between name components',)

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

def rename(align, first, second, splitchar="_"):
    for a in align:
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        #pdb.set_trace()
        for seq in a:
            split_name = seq.id.split('_')
            if second:
                if splitchar == "_":
                    new_seq_name = splitchar.join([split_name[first][0:3], split_name[second][0:3]])
                else:
                    new_seq_name = splitchar.join([split_name[first][0:3], split_name[second][0:3].title()])
            else:
                new_seq_name = split_name[first]
            seq.id, seq.name = new_seq_name, new_seq_name
            new_align.append(seq)
        yield new_align

def main():
    options, args = interface()
    # iterate through all the files to determine the longest alignment
    files = get_files(options.input)
    pos1, pos2 = [eval(i) for i in options.positions.strip().split(',')]
    for count, f in enumerate(files):
        align = AlignIO.parse(f, "nexus")
        new_name = os.path.splitext(os.path.split(f)[1])[0] + '.' + options.format
        outf = os.path.join(options.output, new_name)
        if options.shorten_name:
            for item in rename(align, pos1, pos2, options.splitchar):
                AlignIO.write(item, open(outf, 'w'), options.format)
        else:
            AlignIO.write(align, open(outf, 'w'), options.format)
            
        print count

if __name__ == '__main__':
    main()