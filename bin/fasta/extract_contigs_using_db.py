#!/usr/bin/env python
# encoding: utf-8
"""
File: extract_contigs_using_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 September 2012 14:09 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sqlite3
import argparse
from seqtools.sequence.fasta import FastaReader, FastaWriter
from phyluce.helpers import is_file, is_dir, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            type=is_file,
            action=FullPaths,
            help="""Help text"""
        )
    parser.add_argument(
            "contigs",
            type=is_file,
            action=FullPaths,
            help="""Help text"""
        )
    parser.add_argument(
            "output",
            type=is_dir,
            action=FullPaths,
            help="""Help text"""
        )
    parser.add_argument(
            "--flag",
            action="store_true",
            default=False,
            help="""Help text""",
        )
    return parser.parse_args()


def get_taxa_and_loci(cur):
    results = {}
    cur.execute("PRAGMA table_info(match_map)")
    columns = [i[1] for i in cur.fetchall() if i[1] != 'uce']
    for taxon in columns:
        query = "SELECT {0} FROM match_map WHERE {0} IS NOT NULL".format(taxon)
        cur.execute(query)
        rows = [i[0].split('(')[0] for i in cur.fetchall()]
        results[taxon] = rows
    return results


def parse_fasta_and_write_new_file(results, contigs, output):
    #pdb.set_trace()
    for taxon, rows in results.iteritems():
        outp = FastaWriter(os.path.join(output, "{}.fasta".format(taxon)))
        inp = "{}.contigs.fasta".format(taxon.replace('_', '-'))
        fasta_file = FastaReader(os.path.join(contigs, inp))
        for fasta in fasta_file:
            name = '_'.join(fasta.identifier.lstrip('>').split('_')[:2]).lower()
            if name in rows:
                outp.write(fasta)
        outp.close()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    # get the taxa and loci we're after
    results = get_taxa_and_loci(cur)
    # parse the fasta file and write the UCE contigs to a new file
    parse_fasta_and_write_new_file(results, args.contigs, args.output)

if __name__ == '__main__':
    main()
