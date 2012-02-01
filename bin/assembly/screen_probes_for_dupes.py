#!/usr/bin/env python
# encoding: utf-8

"""
screen_probes_for_dupes.py

Created by Brant Faircloth on 13 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import argparse
from phyluce import lastz
from operator import itemgetter
from collections import defaultdict
from phyluce.helpers import get_name
from phyluce.helpers import get_dupes

def get_args():
    parser = argparse.ArgumentParser(description='Screen a lastz file of self-to-self matches for dupes')
    parser.add_argument('lastz', help='The lastz input')
    return parser.parse_args()

def main():
    args = get_args()
    print get_dupes(args.lastz)

if __name__ == '__main__':
    main()
