#!/usr/bin/env python
# encoding: utf-8

"""
get_fastas_from_match_counts.py

Created by Brant Faircloth on 04 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import sqlite3
import argparse
import ConfigParser
from tools.sequence import fasta
from tools.sequence import transform
from match_contigs_to_probes import is_dir
from match_contigs_to_probes import get_name
from get_match_counts import get_names_from_config

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('contigs', help='The directory containing the contigs to match against probes', type=is_dir)
    parser.add_argument('db', help='The database holding the match and match_map tables from match_contigs_to_probes.py')
    parser.add_argument('config', help='The config file holding the organismal group whose fastas we want')
    parser.add_argument('--output', help='The output file')
    parser.add_argument('--extend-db', dest = 'extend_db', help='The match database to add as an extension')
    parser.add_argument('--extend-dir', dest = 'extend_dir', help='The directory holding extension fastas/contigs')
    parser.add_argument('--notstrict', help = 'The outfile for notstrict data', type=argparse.FileType('w'), default=False)
    return parser.parse_args()

def get_nodes_for_uces(c, organism, uces, extend = False, notstrict = False):
    # get only those UCEs we know are in the set
    uces = [("\'{0}\'").format(u) for u in uces]
    if not extend:
        query = "SELECT lower({0}), uce FROM match_map where uce in ({1})".format(organism, ','.join(uces))
    else:
        query = "SELECT lower({0}), uce FROM extended.match_map where uce in ({1})".format(organism, ','.join(uces))
    c.execute(query)
    rows = c.fetchall()
    node_dict = {node[0].split('(')[0]:[node[1], node[0].split('(')[1].strip(')')] for node in rows if node[0] is not None}
    if notstrict:
        missing = [node[1] for node in rows if node[0] is None]
    else:
        missing = None
    return node_dict, missing
    
def get_coverage(header):
    return '_'.join(header.split('_')[-2:])

def main():
    args = get_args()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.config)
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    if args.extend_db:
        query = "ATTACH DATABASE '{0}' AS extended".format(args.extend_db)
        c.execute(query)
    organisms = get_names_from_config(config, 'Organisms')
    uces = get_names_from_config(config, 'Loci')
    #pdb.set_trace()
    uce_fasta_out = fasta.FastaWriter(args.output)
    for organism in organisms:
        print "\tGetting {0} reads...".format(organism)
        written = []
        # going to need to do something more generic w/ suffixes
        if args.notstrict:
            if not organism.endswith('*'):
                reads = os.path.join(args.contigs, organism.replace('_', '-') + '.contigs.fasta')
                node_dict, missing = get_nodes_for_uces(c, organism, uces, extend = False, notstrict = True)
            elif args.extend_dir:
                # remove the asterisk
                organism = organism.rstrip('*')
                try:
                    reads = os.path.join(args.extend_dir, organism.replace('_', '-') + '.contigs.fasta')
                    assert os.path.exists(reads)
                except AssertionError:
                    reads = os.path.join(args.extend_dir, organism.replace('_', '-') + '.fasta')
                    assert os.path.exists(reads)
                node_dict, missing = get_nodes_for_uces(c, organism, uces, extend = True, notstrict = True)
        else:
            if not organism.endswith('*'):
                reads = os.path.join(args.contigs, organism.replace('_', '-') + '.contigs.fasta')
                node_dict, missing = get_nodes_for_uces(c, organism, uces)
            elif args.extend_dir:
                # remove the asterisk
                organism = organism.rstrip('*')
                try:
                    reads = os.path.join(args.extend_dir, organism.replace('_', '-') + '.contigs.fasta')
                    assert os.path.exists(reads)
                except AssertionError:
                    reads = os.path.join(args.extend_dir, organism.replace('_', '-') + '.fasta')
                    assert os.path.exists(reads)
                node_dict, missing = get_nodes_for_uces(c, organism, uces, extend = True)
            else:
                organism = organism.rstrip('*')
                reads = os.path.join(args.contigs, organism.replace('_', '-') + '.fasta')
                node_dict, missing = get_nodes_for_uces(c, organism, uces)
        for read in fasta.FastaReader(reads):
            name = get_name(read.identifier)
            coverage = get_coverage(read.identifier)
            if name in node_dict.keys():
                uce_seq = fasta.FastaSequence()
                uce_seq.identifier = ">{0}_{1} |{0}|{2}".format(node_dict[name][0], organism, coverage)
                # deal with strandedness because aligners dont, which
                # is annoying
                if node_dict[name][1] == '-':
                    uce_seq.sequence = transform.DNA_reverse_complement(read.sequence)
                else:
                    uce_seq.sequence = read.sequence
                uce_fasta_out.write(uce_seq)
                written.append(str(node_dict[name][0]))
            else:
                pass
        #pdb.set_trace()
        if args.notstrict and missing:
            args.notstrict.write("[{0}]\n".format(organism))
            for name in missing:
                args.notstrict.write("{0}\n".format(name))
                written.append(name)
        assert set(written) == set(uces), "UCE names do not match"
        #assert set(written) == set(uces), pdb.set_trace()
        
    uce_fasta_out.close()

if __name__ == '__main__':
    main()
