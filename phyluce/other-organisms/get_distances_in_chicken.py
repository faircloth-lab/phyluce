#!/usr/bin/env python
# encoding: utf-8

"""
get_distances_in_chicken.py

Created by Brant Faircloth on 03 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import numpy
import sqlite3

import pdb

def get_chromos(c):
    c.execute('''select distinct(chromo) from cons, blast 
        where blast.id = cons.id 
        and blast.matches = 1 
        and duplicate = 0''')
    return c.fetchall()

def get_data_from_db(c, chromo):
    c.execute('''select chromo, cons_start, cons_end from cons, blast 
        where blast.id = cons.id 
        and blast.matches = 1 
        and duplicate = 0 
        and chromo = ? 
        order by cons_start''', chromo)
    return c.fetchall()
    
def main():
    conn = sqlite3.connect('/Users/bcf/git/brant/seqcap/probe.sqlite')
    c = conn.cursor()
    chromos = get_chromos(c)
    lengths = []
    #pdb.set_trace()
    for chromo in chromos:
        rows = get_data_from_db(c, chromo)
        for k,v in enumerate(rows):
            if k == 0:
                old = v[2]
            else:
                new = v[1]
                lengths.append(new - old)
                old = v[2]
    n = numpy.array(lengths)
    print numpy.mean(n)
    print 1.96 * (numpy.std(n, ddof=1)/numpy.sqrt(len(n)))

if __name__ == '__main__':
    main()