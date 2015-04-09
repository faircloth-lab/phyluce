#!/usr/bin/env python
# encoding: utf-8
"""
File: get_binary_match_data.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 September 2012 11:09 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import argparse
import sqlite3
from phyluce.helpers import is_file, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Get binary alignment data from a database of UCE matches""")
    parser.add_argument(
            "db",
            type=is_file,
            action=FullPaths,
            help="""The path to the database"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""The path to the output file"""
        )
    return parser.parse_args()


def get_taxa_from_db(cur):
    cur.execute("PRAGMA table_info(matches);")
    columns = cur.fetchall()
    return [i[1] for i in columns if i[1] != 'uce']


def get_match_data_for_taxon(cur, taxon):
    query = "SELECT uce, {} FROM matches ORDER BY uce".format(taxon)
    cur.execute(query)
    return cur.fetchall()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    taxa = get_taxa_from_db(cur)
    matches = {}
    for taxon in taxa:
        binary = []
        rows = get_match_data_for_taxon(cur, taxon)
        for result in rows:
            locus, present = result
            if present == '1':
                binary.append('1')
            elif present is None:
                binary.append('0')
        matches[taxon] = binary
    outp = open(args.output, 'w')
    outp.write("{0} {1}\n".format(
        len(matches),
        len(matches[matches.keys()[0]])
        ))
    for taxon, binary in matches.iteritems():
        outp.write("{0}\t{1}\n".format(taxon, ''.join(binary)))
    outp.close()


if __name__ == '__main__':
    main()
