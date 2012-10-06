#!/usr/bin/env python
# encoding: utf-8
"""
File: get_bed_from_sqlite.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 July 2012 15:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sqlite3
import argparse
from collections import Counter
from phyluce.helpers import FullPaths

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
            "output",
            help="""The path to store the output (defaults to tablename.bed)""",
            action=FullPaths
        )
    parser.add_argument(
            "--tables",
            nargs='+',
            default=None,
            help="""The table or table of lastz results to convert to BED""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    cur.execute("SELECT * FROM sqlite_master WHERE type='table'")
    if args.tables:
        tables = [t[2] for t in cur.fetchall() if t[2] in args.tables]
        assert len(tables) == len(args.tables), "You have specified tables that are not in the DB"
    else:
        tables = [t[2] for t in cur.fetchall() if t[2]]
    for table in tables:
        outf = open(os.path.join(args.output, "{0}.bed".format(table)), 'w')
        outf.write('''track name=probes description="{0} All UCE probes" useScore=1\n'''.format(table))
        query = '''SELECT name1, strand1, zstart1, end1, name2, probes.id, locus, probe, source
                FROM {0}, probes
                WHERE duplicate = 0
                AND {0}.name2 = probes.oldprobe
                ORDER BY name1, zstart1'''.format(table)
        #pdb.set_trace()
        cur.execute(query)
        matches = cur.fetchall()
        for match in matches:
            chromo, strand, start, end, oldname, dbid, dblocus, dbprobe, dbsource = match
            # parse good stuff from name
            outf.write("{0}\t{1}\t{2}\tprobes-id:{5}|probes-locus:{3}|probes-probe:{4}|probes-oldprobe:{6}\t1000\t{7}\n".format(
                    chromo,
                    start,
                    end,
                    dblocus,
                    dbprobe,
                    dbid,
                    oldname,
                    strand
                ))
        outf.close()



if __name__ == '__main__':
    main()
