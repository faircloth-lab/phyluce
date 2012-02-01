#!/usr/bin/env python
# encoding: utf-8

"""
parse_dbsnp_xml.py

Created by Brant Faircloth on 21 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import argparse
from lxml import etree
from collections import namedtuple
from phyluce.helpers import get_xml_data

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Parse data from dbSNP')
    parser.add_argument('xml', help='The xml file to parse')
    return parser.parse_args()

def main():
    args = get_args()
    snps = get_xml_data(args.xml, True)
    #pdb.set_trace()

if __name__ == '__main__':
    main()
