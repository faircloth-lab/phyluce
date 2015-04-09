#!/usr/bin/env python
# encoding: utf-8
"""
File: format_similarity_data_for_R.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 April 2012 13:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import sqlite3
import argparse

from phyluce.helpers import FullPaths, is_file

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Given a config file of similarity results, enter those to database""")
    parser.add_argument(
            "db",
            action=FullPaths,
            type=is_file,
            help="The path to the config file containing directory mappings"
        )
    parser.add_argument(
            "table",
            type=str,
            choices=['diverge_gaps', 'diverge_no_gaps', 'distance', 'length'],
            help="The table name to summarize"
        )
    parser.add_argument(
            "matrix",
            type=str,
            choices=['complete', 'incomplete'],
            help="Pull values assuming complete/incomplete matrix"
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="The name of the output database",
        )
    parser.add_argument(
            "--exclude",
            type=str,
            help="Taxa to exclude",
        )
    parser.add_argument(
            "--limit",
            type=int,
            default=None,
            help="Return only --limit records",
        )
    return parser.parse_args()


def get_complete_matrix(cur, taxa, null=False, limit=None):
    if not null:
        s = ["{} IS NOT NULL".format(taxon) for taxon in taxa]
    elif null:
        s = ["{} IS NULL".format(taxon) for taxon in taxa]
    s = ' AND '.join(s)
    if not limit:
        query = "SELECT uce FROM distance WHERE {0}".format(s)
    else:
        query = "SELECT uce FROM distance WHERE {0} LIMIT {1}".format(s, limit)
    cur.execute(query)
    #loci = [locus[0] for locus in cur.fetchall()]
    loci = ["'{0}'".format(locus[0]) for locus in cur.fetchall()]
    return loci


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    cur.execute("""pragma table_info(distance)""")
    if args.exclude:
        exclude = [e for e in args.exclude.split(',')]
        exclude.extend(['uce'])
    else:
        exclude = ['uce']
    taxa = [element[1] for element in cur.fetchall() if element[1] not in exclude]
    outfile = open(args.output, 'w')
    outfile.write("uce,length,{0},taxon\n".format(args.table))
    # get only loci having values at all loci
    if args.matrix == 'complete':
        complete = get_complete_matrix(cur, taxa, limit=args.limit)
    else:
        complete = get_complete_matrix(cur, taxa, null=True, limit=args.limit)
    for taxon in taxa:
        if args.matrix == 'complete':
            query = """SELECT {1}.uce, length.{0}, round({1}.{0}, 3)
                FROM length, {1}
                WHERE length.uce = {1}.uce
                AND {1}.uce IN ({2})""".format(taxon, args.table, ', '.join(complete))
        elif args.matrix == 'incomplete':
            query = """SELECT {1}.uce, length.{0}, round({1}.{0}, 3)
                FROM length, {1}
                WHERE length.uce = {1}.uce
                AND {1}.uce NOT IN ({2})""".format(taxon, args.table, ', '.join(complete))
        cur.execute(query)
        data = cur.fetchall()
        for result in data:
            r = [str(i) for i in result]
            s = ','.join(r).replace('None', '')
            outfile.write("{0},{1}\n".format(s, taxon))
    outfile.close()

if __name__ == '__main__':
    main()
