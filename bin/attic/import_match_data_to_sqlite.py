#!/usr/bin/env python
# encoding: utf-8
"""
File: import_match_data_to_sqlite.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 August 2012 13:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import sqlite3
import argparse
from phyluce.helpers import is_dir, FullPaths

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
            "exports",
            type=is_dir,
            action=FullPaths,
            help="""The exported data to import"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    c.execute("""PRAGMA foreign_keys = ON""")
    for file in glob.glob(os.path.join(args.exports, '*.csv')):
        table_name = os.path.basename(file).split('-')[0]
        print "Inserting {0}".format(table_name)
        query = """CREATE TABLE {0} (
                id INTEGER REFERENCES probes,
                score int,
                name1 text,
                strand1 text,
                zstart1 int,
                end1 int,
                length1 int,
                name2 text,
                strand2 text,
                zstart2 int,
                end2 int,
                length2 int,
                diff text,
                cigar text,
                identity text,
                percent_identity float,
                continuity text,
                percent_continuity float,
                coverage text,
                percent_coverage float,
                duplicate int)""".format(table_name)
        c.execute(query)
        for row in open(file, 'rU'):
            rs = row.strip().split(',')
            c.execute("""SELECT id FROM probes WHERE oldprobe = ?""", (rs[7],))
            idx = c.fetchall()
            if len(idx) == 1:
                query = """INSERT INTO {} VALUES ({},{},'{}','{}',{},{},{},'{}','{}',{},{},{},'{}','{}','{}',{},'{}',{},'{}',{},{})""".format(
                    table_name,
                    idx[0][0],
                    rs[1],
                    rs[2],
                    rs[3],
                    rs[4],
                    rs[5],
                    rs[6],
                    rs[7],
                    rs[8],
                    rs[9],
                    rs[10],
                    rs[11],
                    rs[12],
                    rs[13],
                    rs[14],
                    rs[15],
                    rs[16],
                    rs[17],
                    rs[18],
                    rs[19],
                    rs[20],
                    )
                c.execute(query)
            elif len(idx) == 0:
                print "Dropped {}".format(rs[7])
            elif len(idx) > 1:
                print "Probe {} has more than one index".format(rs[7])
    conn.commit()
    conn.close()


if __name__ == '__main__':
    main()
