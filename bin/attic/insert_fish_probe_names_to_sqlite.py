#!/usr/bin/env python
# encoding: utf-8
"""
File: insert_probe_names_to_sqlite.py
Author: Brant Faircloth

Created by Brant Faircloth on 19 July 2012 09:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import re
import sqlite3
import argparse
from collections import defaultdict
from Bio import SeqIO

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            help="""The database to which to add the probe names"""
        )
    parser.add_argument(
            "probes",
            help="""The probe fasta file to enter to the database"""
        )
    return parser.parse_args()

def get_all_probes(probes):
    all_probes = defaultdict(lambda: defaultdict(list))
    for record in SeqIO.parse(open(probes, 'rU'), 'fasta'):
        rs = record.description.split('|')
        name, span = rs[0], rs[1]
        all_probes[name][span].append(record)
    return all_probes



def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    all_probes = get_all_probes(args.probes)
    cur.execute("PRAGMA foreign_keys = ON")
    query = '''CREATE TABLE probes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            locus int,
            probe int,
            source text,
            oldlocus text,
            oldprobe text,
            sequence text
        )'''
    cur.execute(query)
    query = '''CREATE TABLE probeset (
            id int REFERENCES probes,
            probes500 int DEFAULT 0,
            probes1k int DEFAULT 0
        )'''
    cur.execute(query)
    for lpos, locus in enumerate(all_probes):
        spans = sorted(all_probes[locus].keys())
        for spos, span in enumerate(spans):
            probe = all_probes[locus][span][0]
            ns = probe.description.split('|')
            query = '''INSERT INTO probes (locus, probe, source, oldlocus, oldprobe, sequence)
                    values ({0},{1},'{2}','{3}','{4}','{5}')'''.format(
                    lpos + 1,
                    spos + 1,
                    ns[2],
                    ns[0],
                    probe.description,
                    str(probe.seq)
                )
            cur.execute(query)
            if ns[2] == '500':
                cur.execute('''INSERT INTO probeset (id, probes500, probes1k) VALUES (?, 1, 0)''', (cur.lastrowid,))
            elif ns[2] == '1000':
                cur.execute('''INSERT INTO probeset (id, probes500, probes1k) VALUES (?, 0, 1)''', (cur.lastrowid,))
            elif ns[2] == 'both':
                cur.execute('''INSERT INTO probeset (id, probes500, probes1k) VALUES (?, 1, 1)''', (cur.lastrowid,))
    conn.commit()
    cur.close()
    conn.close()

if __name__ == '__main__':
    main()
