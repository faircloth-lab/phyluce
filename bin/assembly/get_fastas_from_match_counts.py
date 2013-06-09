#!/usr/bin/env python
# encoding: utf-8

"""
get_fastas_from_match_counts.py

Created by Brant Faircloth on 04 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import re
import sqlite3
import argparse
import ConfigParser
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phyluce.helpers import is_dir
from phyluce.helpers import get_names_from_config

#import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument(
        'contigs',
        help='The directory containing the contigs to match against probes',
        type=is_dir
    )
    parser.add_argument(
        'db',
        help='The database holding the match and match_map tables from match_contigs_to_probes.py'
    )
    parser.add_argument(
        'config',
        help='The config file holding the organismal group whose fastas we want'
    )
    parser.add_argument(
        '--output',
        help='The output file'
    )
    parser.add_argument(
        '--extend-db',
        dest='extend_db',
        help='The match database to add as an extension'
    )
    parser.add_argument(
        '--extend-dir',
        dest='extend_dir',
        help='The directory holding extension fastas/contigs'
    )
    parser.add_argument(
        '--incomplete-matrix',
        dest='notstrict',
        help='The outfile for incomplete-matrix data',
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
    node_dict = defaultdict()
    missing = []
    for node in rows:
        if node[0] is not None:
            match = re.search('^(node_\d+|comp\d+_c\d+_seq\d+)\(([+-])\)', node[0])
            node_dict[match.groups()[0]] = (node[1], match.groups()[1])
        elif notstrict:
            missing.append(node[1])
        else:
            raise IOError("Complete matrices should have no missing data")
    return node_dict, missing


def find_file(contigs, name):
    extensions = ['.fa', '.fasta', '.contigs.fasta', '.contigs.fa', '.gz', '.fasta.gz', '.fa.gz']
    for ext in extensions:
        reads1 = os.path.join(contigs, name) + ext
        reads2 = os.path.join(contigs, name.replace('-', '_')) + ext
        for reads in [reads1, reads2]:
            if os.path.isfile(reads):
                break
            elif os.path.isfile(reads.lower()):
                reads = reads.lower()
                break
            else:
                reads = None
        if reads is not None:
            break
    if reads is None:
        raise ValueError("Cannot find the a fasta file for {} with any of the extensions ({}) ".format(
            name,
            ', '.join(extensions)
        ))
    return reads


def get_contig_name(header):
    """parse the contig name from the header of either velvet/trinity assembled contigs"""
    match = re.search("^(Node_\d+|comp\d+_c\d+_seq\d+).*", header)
    return match.groups()[0]


def replace_and_remove_bases(regex, seq):
    new_seq_string = str(seq.seq)
    if regex.search(new_seq_string):
        new_seq_string = re.sub(regex, "", new_seq_string)
        print "\tReplaced < 20 ambiguous bases in {0}".format(seq.id)
    new_seq_string = re.sub("^[acgtn]+", "", new_seq_string)
    new_seq_string = re.sub("[acgtn]+$", "", new_seq_string)
    new_seq = Seq(new_seq_string)
    new_seq_record = SeqRecord(new_seq, id=seq.id, name='', description='')
    return new_seq_record


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
    uce_fasta_out = open(args.output, 'w')
    regex = re.compile("[N,n]{1,21}")
    for organism in organisms:
        print "Getting {0} reads...".format(organism)
        written = []
        # going to need to do something more generic w/ suffixes
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
        for seq in SeqIO.parse(open(reads, 'rU'), 'fasta'):
            name = get_contig_name(seq.id).lower()
            if name in node_dict.keys():
                seq.id = "{0}_{1} |{0}".format(node_dict[name][0], organism.rstrip('*'))
                seq.name = ''
                seq.description = ''
                # deal with strandedness because aligners sometimes dont, which
                # is annoying
                if node_dict[name][1] == '-':
                    seq.seq = seq.seq.reverse_complement()
                # Replace any occurrences of <21 Ns in a given sequence with
                # blanks.  These should gap out during alignment. Also, replace
                # leading/trailing lowercase bases from velvet assemblies.
                # Lowercase bases indicate low coverage, and these
                # have been problematic in downstream alignments).
                seq = replace_and_remove_bases(regex, seq)
                uce_fasta_out.write(seq.format('fasta'))
                written.append(str(node_dict[name][0]))
            else:
                pass
        if args.notstrict and missing:
            args.notstrict.write("[{0}]\n".format(organism))
            for name in missing:
                args.notstrict.write("{0}\n".format(name))
                written.append(name)
        assert set(written) == set(uces), "UCE names do not match"
    uce_fasta_out.close()

if __name__ == '__main__':
    main()
