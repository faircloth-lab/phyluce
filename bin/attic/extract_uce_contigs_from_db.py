#!/usr/bin/env python
# encoding: utf-8
"""
File: extract_uce_contigs_from_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2012 16:10 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Given a fasta file of contigs and a database file of matches,
output matches of contigs to UCEs

"""

import os
import argparse
import sqlite3
from phyluce.helpers import is_file, FullPaths
from Bio import SeqIO

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "fasta",
            action=FullPaths,
            type=is_file,
            help="""The fasta file to screen"""
        )
    parser.add_argument(
            "db",
            action=FullPaths,
            type=is_file,
            help="""The database matches to use"""
        )
    parser.add_argument(
            "column",
            type=str,
            help="""The name of the column to use"""
        )
    parser.add_argument(
            "output",
            type=str,
            help="""The name of the output fasta file"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    # get all the UCE records from the db
    query = "SELECT uce, {0} FROM match_map WHERE {0} IS NOT NULL".format(args.column)
    cur.execute(query)
    data = {row[1].split("(")[0]:row[0] for row in cur.fetchall()}
    nodenames = set(data.keys())
    # make sure we don't lose any dupes
    assert len(data) == len(nodenames), "There were duplicate contigs."
    outp = open(args.output, 'w')
    for record in SeqIO.parse(open(args.fasta), 'fasta'):
        name = '_'.join(record.id.split('_')[:2])
        if name.lower() in nodenames:
            record.id = "{0}|{1}".format(data[name.lower()], record.id)
            outp.write(record.format('fasta'))
    outp.close()

if __name__ == '__main__':
    main()