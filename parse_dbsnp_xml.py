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

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Parse data from dbSNP')
    parser.add_argument('xml', help='The xml file to parse')
    return parser.parse_args()

def main():
    args = get_args()
    xml = etree.parse(args.xml)
    validity_terms = set(['byHapMap', 'byOtherPop', 'suspect', 'byFrequency', 
        'by1000G', 'by2Hit2Allele', 'byCluster'])
    print "rsid,type,genotype,het-type,het-value,het-std-error,freq-allele,freq-freq,"+ \
        "freq-sample-size,val-hapmap,val-other-pop,val-freq,val-2hit,val-cluster,"+ \
        "val-1000G,val-suspect"
    for cnt, t in enumerate(xml.getiterator( '{http://www.ncbi.nlm.nih.gov/SNP/docsum}Rs' )):
        rsid = t.get('rsId')
        typ  = t.get('snpType')
        geno = t.get('genotype')
        h = t.find('{http://www.ncbi.nlm.nih.gov/SNP/docsum}Het')
        if h is not None:
            het = h.attrib
        else:
            het = {'type': None, 'value': None, 'stdError': None}
        f = t.find('{http://www.ncbi.nlm.nih.gov/SNP/docsum}Frequency')
        if f is not None:
            freq = f.attrib
        else:
            freq = {'allele': None, 'freq': None, 'sampleSize': None}
        validity = dict(t.find('{http://www.ncbi.nlm.nih.gov/SNP/docsum}Validation').attrib)
        for missing in validity_terms.difference(set(validity.keys())):
            validity[missing] = None
        print "rs{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(rsid, typ, geno, 
            het['type'], het['value'],het['stdError'], freq['allele'],freq['freq'],
            freq['sampleSize'],validity['byHapMap'],validity['byOtherPop'],
            validity['byFrequency'],validity['by2Hit2Allele'],validity['byCluster'],
            validity['by1000G'],validity['suspect'])
        #pdb.set_trace()

if __name__ == '__main__':
    main()