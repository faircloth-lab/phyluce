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
    parser.add_argument(
            "--regex",
            type=str,
            default=None,
            help="""A regular expression to apply to the probe sequences for replacement""",
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


def write_bed_file(outf, chromo, locs, locus):
    outf.write("{0}\t{1}\t{2}\t{3}\t1000\t+\t{1}\t{2}\t100,149,237\n".format(chromo, min(locs), max(locs), locus))


def main():
    args = get_args()
    regex = re.compile(args.regex)
    log = setup_logger()
    for file in glob.glob(os.path.join(args.input, "*lastz*")):
        lz = lastz.Reader(file, long_format=True)
        loci = defaultdict(lambda: defaultdict(list))
        # get output file name
        #outname = os.path.basename(file).split('.')[1].split('_')[-1]
        search_result = re.search('_v_([A-Za-z0-9]+).lastz', os.path.basename(file))
        outname = search_result.groups()[0]
        log.info("Working on {}".format(outname))
        for match in lz:
            locus = re.sub(regex, '', match.name2.split('|')[0].strip())
            for pos in [match.zstart1, match.end1]:
                loci[locus][match.name1].append(pos)
        outf = open(os.path.join(args.output, "{}.bed".format(outname)), 'w')
        outf.write('''track name="uce-v-{0}" description="UCE locus matches to {0}" visibility=2 itemRgb="On"\n'''.format(outname))
        written = set([])
        for locus in sorted(loci.keys()):
            matches = loci[locus]
            for chromo, locs in matches.iteritems():
                if locus in written:
                    log.warn("{0} may have >1 hit".format(locus))
                else:
                    written.add(locus)
                if len(locs) == 2:
                    write_bed_file(outf, chromo, locs, locus)
                else:
                    mn = min(locs)
                    mx = max(locs)
                    dist = mx - mn
                    if dist > 1000:
                        log.warn("Region ({0} bp) is large for {1} at {2}:{3}-{4}".format(dist, locus, match.name1, mn, mx))
                    write_bed_file(outf, chromo, locs, locus)
        outf.close()
        #pdb.set_trace()




if __name__ == '__main__':
    main()
