#!/usr/bin/env python
# encoding: utf-8

"""
helpers.py

Created by Brant Faircloth on 15 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import re
import sys
import glob
import argparse
import shutil
import ConfigParser
#from operator import itemgetter
from collections import defaultdict

from phyluce import lastz
from phyluce.pth import get_all_user_params

#import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_name(header, splitchar = "_", items = 2):
    """use own function vs. import from match_contigs_to_probes - we don't want lowercase"""
    if splitchar:
        return "_".join(header.split(splitchar)[:items]).lstrip('>')
    else:
        return header.lstrip('>')

def get_dupe_matches(lastz_file, splitchar = "|", pos = 1, longfile = False):
    matches = defaultdict(list)
    for lz in lastz.Reader(lastz_file, longfile):
        target_name = get_name(lz.name1, splitchar, pos)
        query_name = get_name(lz.name2, splitchar, pos)
        matches[target_name].append(query_name)
    return matches

def get_dupes(lastz_file, splitchar = "|", pos = 1, longfile = False):
    dupes = set()
    matches = get_dupe_matches(lastz_file, splitchar, pos, longfile)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k, v in matches.iteritems():
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

def get_names_from_config(config, group):
    try:
        return [i[0] for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None

def run_checks(k, v, probes, verbose = True):
    try:
        assert probes[k] >= len(v)               # don't allow more matches than expected
        assert len(set([i[0] for i in v])) == 1  # matches to more than one chromosome
        assert len(set([i[1] for i in v])) == 1  # multiple match directions
        return True
    except AssertionError:
        if verbose:
            print "Probe: ", k
            result_string = ["{0}:{1}-{2}|{3}".format(i[0], i[2], i[3], i[1]) for i in v]
            print "\t\tExpected hits: {0}\n\t\tObserved Hits: {1}\n\t\t\t{2}".format(
                    probes[k],
                    len(v),
                    '\n\t\t\t'.join(result_string)
                )
        return False

def get_matches(lastz_file, splitchar, components, fish = False):
    matches = defaultdict(list)
    probes = defaultdict(int)
    for lz in lastz.Reader(lastz_file, long_format = True):
        # skip silly hg19 mhc haplotypes
        if "hap" in lz.name1:
            print "Skipping: ", lz.name1
        else:
            if fish:
                uce_name = get_name(lz.name2, "_", 1)
                # add 1 because fish probe indexing starts @ 0
                probe_number = int(lz.name2.split('|')[1].split('_')[1]) + 1
            else:
                uce_name = get_name(lz.name2, "|", 1)
                probe_number = int(lz.name2.split(':')[-1])

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


def which(program):
    """ from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python"""
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def snip_if_many_N_bases(regex, chromo, seq, uce, verbose = True):
    """Some genome builds contain long runs of Ns.  Since we're
    slicing reads from these genomes, sometimes these slices contains
    giant runs of Ns.  Remove these by finding the UCE and trimming out
    from the middle to retain the UCE while removing the Ns"""
    # find uce in seq
    uce_start = seq.find(uce)
    uce_end = uce_start + len(uce)
    # slice front
    seq_slice = seq[:uce_start]
    # reverse it - we want first occurence moving 5'
    # from uce start (so first going backwards)
    seq_slice = seq_slice[::-1]
    # search for Ns
    r = regex.search(seq_slice)
    if r:
        new_start = len(seq_slice) - r.start()
    else:
        new_start = 0
    # slice rear
    seq_slice = seq[uce_end:]
    r = regex.search(seq_slice)
    if r:
        new_end = uce_end + r.start()
    else:
        new_end = len(seq)
    seq = seq[new_start:new_end]
    if verbose:
        print "{0} trimmed for > 20 N bases".format(chromo)
    return seq


def get_uce_names_from_probes(probes, regex=None, repl=None, lower=False):
    loci = []
    for line in open(probes, 'rU'):
        if line.startswith('>'):
            ls = line.strip().lstrip('>').split('|')
            if regex and repl:
                locus = re.sub(regex, repl, ls[0])
            else:
                ls[0]
            if not lower:
                loci.append(locus)
            else:
                loci.append(locus.lower())
    # return unique set
    return set(loci)


def get_file_extensions(ftype):
    ext = {
        'fasta': ('.fasta', '.fsa', '.aln', '.fa'),
        'nexus': ('.nexus', '.nex'),
        'phylip': ('.phylip', '.phy'),
        'phylip-relaxed': ('.phylip', '.phy', '.phylip-relaxed'),
        'phylip-sequential': ('.phylip', '.phy', '.phylip-sequential'),
        'clustal': ('.clustal', '.clw'),
        'emboss': ('.emboss',),
        'stockholm': ('.stockholm',)
    }
    return ext[ftype]


def get_alignment_files(log, input_dir, input_format):
    log.info("Getting alignment files")
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def write_alignments_to_outdir(log, outdir, alignments, format):
    log.info('Writing output files')
    for tup in alignments:
        locus, aln = tup
        if aln.trimmed is not None:
            outname = "{}{}".format(
                os.path.join(outdir, locus),
                get_file_extensions(format)[0]
            )
            with open(outname, 'w') as outf:
                outf.write(aln.trimmed.format(format))
        else:
            log.warn("DROPPED {0} from output".format(locus))

def get_contig_header_string():
    return "|".join(get_all_user_params("headers"))
