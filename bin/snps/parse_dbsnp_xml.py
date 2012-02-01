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

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Parse data from dbSNP')
    parser.add_argument('xml', help='The xml file to parse')
    return parser.parse_args()

def get_xml_data(xml, prnt = False):
    xml = etree.parse(xml)
    dbsnp = namedtuple('dbsnp', "rsid,type,genotype,het_type,het_value,het_std_error,freq_allele,freq_freq,"+ \
        "freq_sample_size,val_hapmap,val_other_pop,val_freq,val_2hit,val_cluster,"+ \
        "val_1000G,val_suspect")
    validity_terms = set(['byHapMap', 'byOtherPop', 'suspect', 'byFrequency', 
        'by1000G', 'by2Hit2Allele', 'byCluster'])
    if prnt:
        print "rsid,type,genotype,het-type,het-value,het-std-error,freq-allele,freq-freq,"+ \
            "freq-sample-size,val-hapmap,val-other-pop,val-freq,val-2hit,val-cluster,"+ \
            "val-1000G,val-suspect"
    snps = {}
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
        metadata = [
                rsid, typ, geno, het['type'],het['value'],het['stdError'], freq['allele'],freq['freq'],
                freq['sampleSize'],validity['byHapMap'],validity['byOtherPop'],
                validity['byFrequency'],validity['by2Hit2Allele'],validity['byCluster'],
                validity['by1000G'],validity['suspect']
                ]
        if prnt:
            print "rs{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(rsid, typ, geno, 
                het['type'],het['value'],het['stdError'], freq['allele'],freq['freq'],
                freq['sampleSize'],validity['byHapMap'],validity['byOtherPop'],
                validity['byFrequency'],validity['by2Hit2Allele'],validity['byCluster'],
                validity['by1000G'],validity['suspect'])
        else:
            snps[rsid] = dbsnp._make(metadata)
    return snps

def main():
    args = get_args()
    snps = get_xml_data(args.xml, True)
    #pdb.set_trace()

if __name__ == '__main__':
    main()