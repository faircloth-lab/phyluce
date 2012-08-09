#!/usr/bin/env python
# encoding: utf-8

"""
get_fastas_from_match_counts.py

Created by Brant Faircloth on 04 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import re
import sys
import sqlite3
import argparse
import ConfigParser
from seqtools.sequence import fasta
from seqtools.sequence import transform
from phyluce.helpers import is_dir
from phyluce.helpers import get_name
from phyluce.helpers import get_names_from_config

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('contigs',
            help='The directory containing the contigs to match against probes',
            type=is_dir
        )
    parser.add_argument('db',
            help='The database holding the match and match_map tables from match_contigs_to_probes.py'
        )
    parser.add_argument('config',
            help='The config file holding the organismal group whose fastas we want'
        )
    parser.add_argument('--output',
            help='The output file'
        )
    parser.add_argument('--extend-db',
            dest='extend_db',
            help='The match database to add as an extension'
        )
    parser.add_argument('--extend-dir',
            dest='extend_dir',
            help='The directory holding extension fastas/contigs'
        )
    parser.add_argument('--notstrict',
            help='The outfile for notstrict data',
            type=argparse.FileType('w'),
            default=False
        )
    return parser.parse_args()


def get_nodes_for_uces(c, organism, uces, extend=False, notstrict=False):
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


def find_file(contigs, name):
    extensions = ['.fa', '.fasta', '.contigs.fasta', '.contigs.fa', '.gz', '.fasta.gz', '.fa.gz']
    for ext in extensions:
            reads = os.path.join(contigs, name) + ext
            if os.path.isfile(reads):
                break
            elif os.path.isfile(reads.lower()):
                reads = reads.lower()
                break
            else:
                reads = None
    if reads is None:
        raise ValueError("Cannot find the a fasta file for {} with any of the extensions ({}) ".format(
                name,
                ', '.join(extensions)
            )
        )
    return reads


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
    regex = re.compile("[N,n]{1,21}")
    for organism in organisms:
        print "Getting {0} reads...".format(organism)
        written = []
        # going to need to do something more generic w/ suffixes
        #pdb.set_trace()
        name = organism.replace('_', '-')
        if args.notstrict:
            if not organism.endswith('*'):
                reads = find_file(args.contigs, name)
                node_dict, missing = get_nodes_for_uces(c, organism, uces, extend=False, notstrict=True)
            elif args.extend_dir:
                # remove the asterisk
                name = name.rstrip('*')
                reads = find_file(args.extend_dir, name)
                node_dict, missing = get_nodes_for_uces(c, organism.rstrip('*'), uces, extend=True, notstrict=True)
        else:
            if not name.endswith('*'):
                reads = find_file(args.contigs, name)
                node_dict, missing = get_nodes_for_uces(c, organism, uces)
            elif name.endswith('*') and args.extend_dir:
                # remove the asterisk
                name = name.rstrip('*')
                reads = find_file(args.extend_dir, name)
                node_dict, missing = get_nodes_for_uces(c, organism.rstrip('*'), uces, extend=True)
        for read in fasta.FastaReader(reads):
            name = get_name(read.identifier).lower()
            coverage = get_coverage(read.identifier)
            if name in node_dict.keys():
                uce_seq = fasta.FastaSequence()
                uce_seq.identifier = ">{0}_{1} |{0}|{2}".format(node_dict[name][0], organism.rstrip('*'), coverage)
                # deal with strandedness because aligners dont, which
                # is annoying
                if node_dict[name][1] == '-':
                    uce_seq.sequence = transform.DNA_reverse_complement(read.sequence)
                else:
                    uce_seq.sequence = read.sequence
                # replace any occurrences of <21 Ns
                if regex.search(uce_seq.sequence):
                    uce_seq.sequence = re.sub(regex, "", uce_seq.sequence)
                    print "\tReplaced < 20 ambiguous bases in {0}".format(uce_seq.identifier.split(' ')[0])
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
