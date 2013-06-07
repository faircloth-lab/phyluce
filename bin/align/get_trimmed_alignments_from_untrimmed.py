#!/usr/bin/env python
# encoding: utf-8
"""
File: trim_align_only.py
Author: Brant Faircloth

Created by Brant Faircloth on 06 May 2012 14:05 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Trim the edges of alignments files to remove nasties.

Usage:

    python ~/Git/brant/phyluce/bin/align/trim_align_only.py \
        sate-fasta-untrimmed \
        sate-nexus-trimmed \
        --output-format nexus \
        --multiprocessing

"""

import os
import sys
import glob
import argparse
import multiprocessing
from phyluce.helpers import is_dir, FullPaths, get_file_extensions
from phyluce.generic_align import GenericAlign

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Use the PHYLUCE trimming algorithm to trim existing alignments""")
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""The directory containing alignments"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""The directory to contain trimmed alignments"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument(
            "--output-format",
            dest="output_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The output alignment format""",
        )
    parser.add_argument(
            "--window",
            type=int,
            default=20,
            help="Sliding window size for trimming"
        )
    parser.add_argument(
            "--threshold",
            type=float,
            default=0.75,
            help="Threshold cutoff for trimming"
        )
    parser.add_argument(
            "--proportion",
            type=float,
            default=0.65,
            help="Proportional removal of gaps"
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""Process alignments with X number of --cores""",
            )
    return parser.parse_args()


def get_and_trim_alignments(params):
    trimming_params, align_file = params
    input_format, window, threshold, proportion = trimming_params
    #pdb.set_trace()
    name = os.path.basename(os.path.splitext(align_file)[0])
    aln = GenericAlign(align_file)
    # call private method to read alignment into alignment object
    try:
        aln._read(input_format)
        # dont return consensus
        aln.trim_alignment(
                method='running',
                window_size=window,
                threshold=threshold,
                proportion=proportion
            )
        sys.stdout.write(".")
        sys.stdout.flush()
        return (name, aln)
    except ValueError, e:
        if e.message == 'No records found in handle':
            return (name, False)
        else:
            raise ValueError('Something is wrong with alignment {0}'.format(name))


def write_alignments_to_outdir(outdir, alignments, format):
    print '\nWriting output files...'
    for tup in alignments:
        locus, aln = tup
        if aln.trimmed_alignment is not None:
            outname = "{}{}".format(
                    os.path.join(outdir, locus),
                    get_file_extensions(format)[0]
                )
            outf = open(outname, 'w')
            outf.write(aln.trimmed_alignment.format(format))
            outf.close()
        else:
            print "\tSkipped writing {0}, there was no record".format(locus)


def main():
    args = get_args()
    alignments = []
    for ftype in get_file_extensions(args.input_format):
        alignments.extend(glob.glob(os.path.join(args.input, "*{}".format(ftype))))
    # package up needed arguments for map()
    package = [args.input_format, args.window, args.threshold, args.proportion]
    params = zip([package] * len(alignments), alignments)
    # print some output for user
    sys.stdout.write('Trimming')
    sys.stdout.flush()
    # if --multprocessing, use Pool.map(), else use map()
    # can also extend to MPI map, but not really needed on multicore
    # machine
    if args.cores > 1:
        pool = multiprocessing.Pool(args.cores - 1)
        alignments = pool.map(get_and_trim_alignments, params)
    else:
        alignments = map(get_and_trim_alignments, params)
    write_alignments_to_outdir(args.output, alignments, args.output_format)


if __name__ == '__main__':
    main()
