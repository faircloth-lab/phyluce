#!/usr/bin/env python
# encoding: utf-8

"""
get_probes_matches_from_lastz.py

Created by Brant Faircloth on 23 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import glob
import sqlite3
import argparse
import textwrap
import bx.seq.twobit
from seqcap.lib import lastz
from operator import itemgetter
from collections import defaultdict
from tools.sequence import transform
from uce_helpers import get_name
from uce_helpers import get_dupes
from get_fake_velvet_contigs_from_genomes import get_matches
from get_fake_velvet_contigs_from_genomes import run_checks

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('lastz', help='The lastz output')
    parser.add_argument('--name-components', dest = 'components', help = 'number of parts in the name', default = 2, type = int)
    parser.add_argument('--splitchar', help = 'The name character on which to split', default = "_", type = str)
    parser.add_argument('--dupefile', help='The path to a lastz file of lastz-against-self results')
    return parser.parse_args()

def main():
    args = get_args()
    if args.dupefile:
        dupes = get_dupes(args.dupefile)
    else:
        dupes = None
    matches, probes = get_matches(args.lastz, args.splitchar, args.components)
    pdb.set_trace()

if __name__ == '__main__':
    main()