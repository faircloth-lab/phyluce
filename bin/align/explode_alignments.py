#!/usr/bin/env python
# encoding: utf-8
"""
File: explode_alignments.py
Author: Brant Faircloth

Created by Brant Faircloth on 29 January 2013 10:01 PST (-0800)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import glob
import argparse
import ConfigParser
from Bio import AlignIO
from phyluce.helpers import get_file_extensions, is_dir, is_file, FullPaths


import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""Input folder of alignments"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""Output folder of fasta files"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format"""
        )
    parser.add_argument(
            "--conf",
            action=FullPaths,
            type=is_file,
            help="""Config file for renaming""",
        )
    parser.add_argument(
            "--section",
            type=str,
            help="""Section of config file to use""",
        )
    parser.add_argument(
            "--exclude",
            type=str,
            nargs='+',
            default = [],
            help="""Taxa/taxon to exclude""",
        )
    parser.add_argument(
            "--by-taxon",
            action="store_true",
            default=False,
            help="""Explode file by taxon instead of by-locus""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    extensions = get_file_extensions(input_format)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    # ensure we collapse duplicate filenames
    return list(set(files))


def main():
    args = get_args()
    files = get_files(args.input, args.input_format)
    if args.conf:
        conf = ConfigParser.ConfigParser()
        conf.read(args.conf)
        names = {item[0].replace(' ', '_'):item[1] for item in conf.items(args.section)}
        print "Original taxon count = ", len(names.keys())
        for taxon in args.exclude:
            del names[taxon]
    #pdb.set_trace()
    if args.by_taxon:
        d = {}
        for file in files:
            sys.stdout.write('.')
            sys.stdout.flush()
            basename = os.path.basename(file)
            locus = os.path.splitext(basename)[0]
            aln = AlignIO.read(file, args.input_format)
            for taxon in aln:
                name = taxon.id.replace(locus, '').lstrip('_')
                if name not in args.exclude:
                    try:
                        shortname = names[name]
                    except:
                        shortname = name
                if shortname not in d.keys():
                    new_file = shortname + ".fasta"
                    d[shortname] = open(os.path.join(args.output, new_file), 'w')
                seq = str(taxon.seq).replace('-', '')
                seq = str(taxon.seq).replace('?', '')
                if not len(seq) == 0:
                    d[shortname].write(">{0}\n{1}\n".format(taxon.id, seq))
        for k, v in d.iteritems():
            v.close()
    else:
        for file in files:
            sys.stdout.write('.')
            sys.stdout.flush()
            basename = os.path.basename(file)
            locus = os.path.splitext(basename)[0]
            new_file = locus + ".fasta"
            taxon_count = []
            #pdb.set_trace()
            outp = open(os.path.join(args.output, new_file), 'w')
            aln = AlignIO.read(file, args.input_format)
            count = 0
            for taxon in aln:
                name = taxon.id.replace(locus, '').lstrip('_')
                if name not in args.exclude:
                    try:
                        shortname = names[name]
                    except:
                        shortname = name
                    seq = str(taxon.seq).replace('-', '')
                    seq = str(taxon.seq).replace('?', '')
                    if not len(seq) == 0:
                        outp.write(">{0}\n{1}\n".format(shortname, seq))
                        count += 1
                    else:
                        print locus
            taxon_count.append(count)
            outp.close()
        print "\n"
        print "Final taxon count = ", set(taxon_count)

if __name__ == '__main__':
    main()
