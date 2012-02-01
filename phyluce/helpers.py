#!/usr/bin/env python
# encoding: utf-8

"""
helpers.py

Created by Brant Faircloth on 15 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import argparse
from seqcap.lib import lastz
from operator import itemgetter
from collections import defaultdict

import pdb

def get_name(header, splitchar = "_", items = 2):
    """use own function vs. import from match_contigs_to_probes - we don't want lowercase"""
    if splitchar:
        return "_".join(header.split(splitchar)[:items]).lstrip('>')
    else:
        return header.lstrip('>')
        
def get_matches(lastz_file, splitchar = "|", pos = 1):
    matches = defaultdict(list)
    for lz in lastz.Reader(lastz_file):
        target_name = get_name(lz.name1, splitchar, pos)
        query_name = get_name(lz.name2, splitchar, pos)
        matches[target_name].append(query_name)
    return matches

def get_dupes(lastz_file, splitchar = "|", pos = 1):
    dupes = set()
    matches = get_matches(lastz_file, splitchar, pos)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k,v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    return dupes

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_names_from_config(config, group):
    try:
        return [i[0] for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None

def run_checks(k, v, probes, verbose = True):
    try:
        # make sure we have as many matches as probes
        assert len(v) == probes[k]              # more matches than expected
        assert len(set([i[0] for i in v])) == 1 # matches to more than one chromosome
        assert len(set([i[1] for i in v])) == 1 # multiple match directions
        return True
    except AssertionError:
        if verbose:print "Probe: ", k
        result_string = ["{0}:{1}-{2}|{3}".format(i[0], i[2], i[3], i[1]) for i in v]
        if verbose:print "\t\tExpected hits: {0}\n\t\tObserved Hits: {1}\n\t\t\t{2}".format(probes[k], len(v), '\n\t\t\t'.join(result_string))
        return False

def get_matches(lastz_file, splitchar, components, fish = False):
    matches = defaultdict(list)
    probes = defaultdict(int)
    for lz in lastz.Reader(lastz_file, long_format = True):
        # skip silly hg19 mhc haplotypes
        if "hap" in lz.name1:
            print "Skipping: ", lz.name1
        else:
            if not fish:
                uce_name = get_name(lz.name2, "|", 1)
                probe_number = int(lz.name2.split(':')[-1])
            else:
                uce_name = get_name(lz.name2, "_", 1)
                # add 1 because fish probe indexing starts @ 0
                probe_number = int(lz.name2.split('|')[1].split('_')[1]) + 1
            #pdb.set_trace()
            if probe_number > probes[uce_name]:
                probes[uce_name] = probe_number
            matches[uce_name].append([get_name(lz.name1, splitchar = splitchar, items = components), lz.strand2, lz.zstart1, lz.end1])
    return matches, probes

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
