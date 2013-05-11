#!/usr/bin/env python
# encoding: utf-8
"""
File: get_probe_sets.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 August 2012 15:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sqlite3
import argparse
from seqtools.sequence import fasta
from phyluce.helpers import is_file

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            help="""The database containing tables to convert to BED"""
        )
    parser.add_argument(
            "--add",
            type=is_file,
            help="""The path to store the output fasta"""
        )
    parser.add_argument(
            "--add-name",
            dest='name',
            type=is_file,
            help="""The path to store the output fasta"""
        )
    parser.add_argument(
            "--build",
            action="store_true",
            default=False,
            help="""Help text""",
        )
    parser.add_argument(
            "--rebuild",
            action="store_true",
            default=False,
            help="""Help text""",
        )
    return parser.parse_args()


def create_probeset_table(args, conn, c):
    try:
        # create main probe column
        c.execute("""PRAGMA foreign_keys = ON""")
        c.execute("""CREATE TABLE probeset (id int REFERENCES probes, probes12k int DEFAULT 0)""")
        c.execute("""INSERT INTO probeset SELECT id, 1 from probes""")
    except sqlite3.OperationalError, msg:
        if 'table probeset already exists' in msg.message:
            if args.rebuild:
                c.execute("""DROP TABLE probeset""")
                create_probeset_table(args, conn, c)
            else:
                print "Looks like the probeset table exists. " \
                + "Use --rebuild to DROP and start fresh."


def add_additional_columns(args, conn, c):
    assert args.name is not None, "You need to include --add-name to add a table"
    query = """ALTER TABLE probeset ADD COLUMN {0} int DEFAULT 0""".format(args.name)
    c.execute(query)
    for seq in fasta.FastaReader(args.add):
        locus = seq.identifier.lstrip('>').split('|')[0]
        query = """SELECT id, locus, probe, source, sequence, oldprobe
                FROM probes
                WHERE oldprobe LIKE '%{0}%'""".format(locus)
        c.execute(query)
        rows = c.fetchall()
        hit = False
        for row in rows:
            idx, locus, probe, source, sequence, oldlocus = row
            if seq.sequence in sequence:
                hit = True
                query = """UPDATE probeset set {0} = 1 WHERE id = {1}""".format(
                        args.name,
                        idx
                        )
                c.execute(query)
        if not hit:
            print "Miss: {0}".format(seq.identifier)
            #pdb.set_trace()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    if args.build or args.rebuild:
        create_probeset_table(args, conn, c)
    add_additional_columns(args, conn, c)
    conn.commit()
    c.close()
    conn.close()


if __name__ == '__main__':
    main()
