#!/usr/bin/env python
# encoding: utf-8
"""
File: get_locus_bed_from_lastz_files.py
Author: Brant Faircloth

Created by Brant Faircloth on 31 January 2013 11:01 PST (-0800)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description: Given an input directory of lastz files, output a
directory of BED files giving the locus locations

"""

import os
import re
import glob
import logging
import argparse
from collections import defaultdict
from phyluce.helpers import FullPaths, is_dir
from phyluce import lastz

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""get_locus_bed_from_lastz""")
    parser.add_argument(
            "input",
            type=is_dir,
            action=FullPaths,
            help="""The input directory containing lastz files"""
        )
    parser.add_argument(
            "output",
            type=is_dir,
            action=FullPaths,
            help="""The output directory to hold BED-formatted files""",
        )
    return parser.parse_args()


def setup_logger():
    # get a logger and give it a name
    logger = logging.getLogger("PHYLUCE")
    formatter = logging.Formatter('%(name)-5s: %(levelname)-8s %(message)s')
    # setup handler
    hdlr = logging.StreamHandler()
    hdlr.setFormatter(formatter)
    # add handlder and set logging level
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    return logger


def write_bed_file(outf, chromo, start, end, probe):
    outf.write("{0}\t{1}\t{2}\t{3}\t1000\t+\t{1}\t{2}\t100,149,237\n".format(chromo, start, end, probe))


def main():
    args = get_args()
    log = setup_logger()
    for file in glob.glob(os.path.join(args.input, "*lastz*")):
        lz = lastz.Reader(file, long_format=True)
        probes = defaultdict(list)
        # get output file name
        #outname = os.path.basename(file).split('.')[1].split('_')[-1]
        search_result = re.search('_v_([A-Za-z0-9]+).lastz', os.path.basename(file))
        outname = search_result.groups()[0]
        log.info("Working on {}".format(outname))
        outf = open(os.path.join(args.output, "{}.bed".format(outname)), 'w')
        outf.write('''track name="uce-v-{0}" description="UCE probe matches to {0}" visibility=2 itemRgb="On"\n'''.format(outname))
        written = set([])
        for match in lz:
            probe = match.name2.split('|')[0].strip()
            probes[probe].append([match.name1, match.zstart1, match.end1])
        #pdb.set_trace()
        for probe in sorted(probes.keys()):
            for match in probes[probe]:
                chromo, start, end = match
                if probe in written:
                    log.warn("{0} may have >1 hit".format(probe))
                else:
                    written.add(probe)
                write_bed_file(outf, chromo, start, end, probe)
        outf.close()
        #pdb.set_trace()

if __name__ == '__main__':
    main()
