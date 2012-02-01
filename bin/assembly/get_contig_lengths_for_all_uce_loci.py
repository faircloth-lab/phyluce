#!/usr/bin/env python
# encoding: utf-8

"""
get_contig_lengths_for_all_uce_loci.py

Created by Brant Faircloth on 02 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import sqlite3
import argparse
import ConfigParser
from collections import defaultdict
from phyluce.helpers import is_dir
from phyluce.helpers import get_names_from_config
from seqtools.sequence import fasta

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get the average contig length, by species, of contigs matching UCEs')
    parser.add_argument('fasta', help='The directory holding contig files (fasta format)', type=is_dir)
    parser.add_argument('db', help='The database containing match information')
    parser.add_argument('config', help='The config file holding the organismal group whose fastas we want')
    parser.add_argument('group', help='The config group whose results you want', type=str)
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def get_matching_node_names(c, org):
    query = "SELECT {0} FROM match_map WHERE {0} IS NOT NULL".format(org)
    data = c.execute(query)
    return set([node[0].upper().split('(')[0] for node in c.execute(query)])

def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.config)
    organisms = get_names_from_config(config, args.group)
    excludes = get_names_from_config(config, 'Excludes')
    if excludes:
        organisms = [org for org in organisms if org not in excludes]
    args.output.write("org\tcontigs\tavg len\n")
    for org in organisms:
        # skip extended data, which are typically from genome-enabled orgs,
        # not capture data
        if not org.endswith('*'):
            # get the uce-matching node names from the db
            matching_nodes = get_matching_node_names(c, org)
            # parse the contig file for the organism, and return contig
            # lengths
            f = os.path.join(args.fasta, "{0}.{1}".format(org.replace('_','-'),'contigs.fasta'))
            records = fasta.FastaReader(f)
            contig_lens = [len(seq) for seq in records 
                if '_'.join(seq.identifier.strip('>').split('_')[0:2]) in matching_nodes]
            # write the average contig length of contigs matching UCEs
            args.output.write("{0}\t{1}\t{2}\n".format(org, len(contig_lens), float(sum(contig_lens))/len(contig_lens)))
            

if __name__ == '__main__':
    main()
    
