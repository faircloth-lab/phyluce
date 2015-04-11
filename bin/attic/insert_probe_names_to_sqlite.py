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
from collections import Counter, defaultdict

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
            "tablename",
            help="""The table to create"""
        )
    parser.add_argument(
            "probes1",
            help="""The probe fasta file to enter to the database"""
        )
    parser.add_argument(
            "probes2",
            help="""Another probe fasta from which to get original source"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    # get the original source of each probe
    if args.probes2:
        probes2 = defaultdict(list)
        loci2 = []
        probe2_file = open(args.probes2, 'rU')
        for line in probe2_file:
            if line.startswith('>'):
                locus_name = line.split('|')
                loci2.append(locus_name[0].strip('>'))
            seq = probe2_file.next().strip()
            probes2[seq].extend([locus_name[0].strip('>'), locus_name[1]])
        probe2_file.close()
    #pdb.set_trace()
    # make oldname dict
    loci2 = list(set(loci2))
    loci2_dict = {v:k for k,v in enumerate(loci2)}
    query = '''CREATE TABLE {0} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            locus int,
            probe int,
            source text,
            oldlocus text,
            oldprobe text,
            sequence text
        )'''.format(args.tablename)
    #pdb.set_trace()
    cur.execute(query)
    probe1_file = open(args.probes1, 'rU')
    regex = re.compile("_(\d+)$")
    for line in probe1_file:
        if line.startswith('>'):
            full_name = line.strip().split('|')[0]
            locus_name = re.sub("_\d+$", "", full_name.lstrip('>'))
            probe_number = int(regex.search(full_name).groups()[0])
            seq = probe1_file.next().strip()
            # get original metadata
            oldname, source = probes2[seq]
            # make sure meta and locus name match
            assert oldname == locus_name, "metadata and locus name do not match"
            query = '''INSERT INTO {0} (locus, probe, source, oldlocus, oldprobe, sequence)
                    values ({1},{2},'{3}','{4}','{5}','{6}')'''.format(
                    args.tablename,
                    loci2_dict[locus_name],
                    probe_number,
                    source,
                    locus_name,
                    full_name,
                    seq
                )
            #pdb.set_trace()
            cur.execute(query)
    conn.commit()
    cur.close()
    conn.close()

if __name__ == '__main__':
    main()
