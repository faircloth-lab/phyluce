#!/usr/bin/env python
# encoding: utf-8

"""
add_gaps_for_missing_taxa.py

Created by Brant Faircloth on 05 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import copy
import glob
import argparse
import ConfigParser
from Bio import AlignIO
from collections import defaultdict
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument('input', help = "Alignment files to process")
    parser.add_argument('output', nargs='?', default=sys.stdout, help = "Output dir for alignment files")
    parser.add_argument('config', help='The config file containing match information')
    parser.add_argument('notstrict', help = "Name of conf file containing notstrict species and data")
    parser.add_argument('--genera', dest='genera', action='append', help = "Name of conf file containing notstrict species and data")
    parser.add_argument('--min-taxa', dest='min_taxa', help='''The minimum number of taxa to keep
    (default = 3)''', default = 3, type = int)
    return parser.parse_args()

def get_names_from_config(config, group):
    try:
        return [i[0].rstrip('*') for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None

def add_gaps_to_align(organisms, missing, align, genera = False, min_taxa = 3):
    local_organisms = copy.deepcopy(organisms)
    for a in align:
        if len(a) < min_taxa:
            new_align = None
            break
        elif len(a) >= min_taxa:
            #pdb.set_trace()
            new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
            overall_length = len(a[0])
            for seq in a:
                if genera and any(sp for sp in genera if sp in seq.name):
                    new_seq_name = '_'.join(seq.name.split('_')[-1:])
                else:
                    new_seq_name = '_'.join(seq.name.split('_')[-2:])
                new_align.add_sequence(new_seq_name, str(seq.seq))
                local_organisms.remove(new_seq_name)
            for org in local_organisms:
                if genera and any(sp for sp in genera if sp in seq.name):
                    loc = '_'.join(seq.name.split('_')[:-1])
                else:
                    loc = '_'.join(seq.name.split('_')[:-2])
                assert loc in missing[org], "Locus missing"
                new_align.add_sequence(org, '-' * overall_length)
    return new_align

def get_missing_loci_from_conf_file(config):
    missing = defaultdict(list)
    for sec in config.sections():
    #    pdb.set_trace()
        for item in config.items(sec):
            missing[sec].append(item[0])
    return missing

def main():
    args = get_args()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.config)
    notstrict = ConfigParser.RawConfigParser(allow_no_value=True)
    notstrict.read(args.notstrict)
    missing = get_missing_loci_from_conf_file(notstrict)
    organisms = get_names_from_config(config, 'Organisms')
    for count, nex in enumerate(glob.glob(os.path.join(args.input, '*.nex'))):
        align = AlignIO.parse(nex, "nexus")
        new_align = add_gaps_to_align(organisms, missing, align, args.genera,
                args.min_taxa)
        if new_align is not None:
            outf = os.path.join(args.output, os.path.basename(nex))
            AlignIO.write(new_align, open(outf, 'w'), 'nexus')
            print count
        else:
            print "{0} Skipping {1} - < {2} taxon".format(count, nex,
                    args.min_taxa)

if __name__ == '__main__':
    main()
        
