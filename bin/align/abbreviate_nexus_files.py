"""
File: abbreviate_nexus_files.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 February 2012 12:02 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:  Abbreviate taxon names used in a Nexus file

"""

import os
import sys
import glob
import numpy
import shutil
import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment

from phyluce.helpers import is_dir

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('nexus', help='The directory containing the nexus files', type=is_dir)
    parser.add_argument('output', help='''The directory that will contain the
    abbreviated Nexus files''', type=is_dir)
    return parser.parse_args()

def abbreviator(taxon_list):
    abbr = {}
    for taxa in taxon_list:
        cap_taxa = taxa.split('_')
        cap_taxa[1] = cap_taxa[1].capitalize()
        cap_taxa[0] = cap_taxa[0].lower()
        nn = "{}{}".format(cap_taxa[0][0:3], cap_taxa[1][0:3])
        for i in xrange(1, 10):
            if not "{}{}".format(nn, i) in abbr.keys():
                abbr["{}{}".format(nn, i)] = taxa
                break
    # reverse dict
    return {v:k for k,v in abbr.iteritems()}

def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.nexus)
    old_names = set()
    for f in files:
        for align in AlignIO.parse(f, 'nexus'):
            for seq in list(align):
                old_names.update([seq.name])
    #pdb.set_trace()
    name_map = abbreviator(old_names)
    for count, f in enumerate(files):
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        #filename = os.path.basename(f)
        #chromo_name = filename.split('.')[0]
        for align in AlignIO.parse(f, 'nexus'):
            for seq in list(align):
                new_seq_name = name_map[seq.name]
                new_align.add_sequence(new_seq_name, str(seq.seq))
        #pdb.set_trace()
        outf = os.path.join(args.output, os.path.split(f)[1])
        try:
            AlignIO.write(new_align, open(outf, 'w'), 'nexus')
        except ValueError:
            pdb.set_trace()
        print count

if __name__ == '__main__':
    main()
