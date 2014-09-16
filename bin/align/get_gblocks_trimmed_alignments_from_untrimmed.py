#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 June 2014 13:28 PDT (-0700)
"""

import os
import sys
import glob
import argparse
import subprocess
import multiprocessing

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

from phyluce.log import setup_logging
from phyluce.helpers import FullPaths, CreateDir, is_dir, get_file_extensions

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Use GBLOCKS to trim existing alignments in parallel""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--alignments",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="""The directory containing alignments to be trimmed."""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The directory in which to store the resulting alignments."""
    )
    parser.add_argument(
        "--input-format",
        choices=["fasta", "nexus", "phylip", "clustal", "emboss", "stockholm"],
        default="fasta",
        help="""The input alignment format.""",
    )
    parser.add_argument(
        "--output-format",
        choices=["fasta", "nexus", "phylip", "clustal", "emboss", "stockholm"],
        default="nexus",
        help="""The output alignment format.""",
    )
    parser.add_argument(
        "--b1",
        type=float,
        default=0.5,
        help="""The GBLOCKS -b1 proportion.""",
    )
    parser.add_argument(
        "--b2",
        type=float,
        default=0.85,
        help="""The GBLOCKS -b2 proportion.""",
    )
    parser.add_argument(
        "--b3",
        type=int,
        default=8,
        help="""The GBLOCKS -b3 integer value.""",
    )
    parser.add_argument(
        "--b4",
        type=int,
        default=10,
        help="""The GBLOCKS -b4 integer value.""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""Process alignments in parallel using --cores for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )
    args = parser.parse_args()
    assert args.b4 >= 2, "--b4 must be ≥ 2"
    return args


def get_and_trim_alignments(params):
    args, align_file = params
    name = os.path.splitext(os.path.basename(align_file))[0]
    try:
        # determine the number of sequences in align
        with open(align_file, 'rU') as input_aln:
            taxa = len(AlignIO.read(input_aln, args.input_format))
        b1 = int(round(args.b1 * taxa)) + 1
        b2 = int(round(args.b2 * taxa))
        if b2 < b1:
            b2 = b1
        cmd = [
            "gblocks",
            align_file,
            "-t=DNA",
            "-b1={}".format(b1),
            "-b2={}".format(b2),
            "-b3={}".format(args.b3),
            "-b4={}".format(args.b4),
            "-b5=h",
            "-p=n"
        ]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        # read new alignment
        trimmed_aln_file = "{}-gb".format(align_file)
        with open(trimmed_aln_file, 'rU') as trimmed_aln:
            aln = AlignIO.read(trimmed_aln, "fasta", alphabet=Gapped(IUPAC.ExtendedIUPACDNA()))
        # remove the qblocks file
        os.remove(trimmed_aln_file)
        sys.stdout.write(".")
        sys.stdout.flush()
        return name, aln
    except ValueError, e:
        if e.message == "No records found in handle":
            return name, False
    except:
        return name, False

def write_gblocks_alignments_to_outdir(log, outdir, alignments, format):
    log.info('Writing output files')
    for tup in alignments:
        locus, aln = tup
        if aln:
            outname = "{}{}".format(
                os.path.join(outdir, locus),
                get_file_extensions(format)[0]
            )
            #pdb.set_trace()
            try:
                with open(outname, 'w') as outf:
                    AlignIO.write(aln, outf, format)
            except ValueError:
                log.warn("Unable to write {} - alignment too short".format(locus))
        else:
            log.warn("Missing information for locus {}".format(locus))


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    alignments = []
    log.info("Getting aligned sequences for trimming")
    for ftype in get_file_extensions(args.input_format):
        alignments.extend(glob.glob(os.path.join(args.alignments, "*{}".format(ftype))))
    # package up needed arguments for map()
    params = [[args, aln] for aln in alignments]
    log.info("Alignment trimming begins.")
    # if --multprocessing, use Pool.map(), else use map()
    # can also extend to MPI map, but not really needed on multicore
    # machine
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores - 1)
        alignments = pool.map(get_and_trim_alignments, params)
    else:
        alignments = map(get_and_trim_alignments, params)
    print("")
    # drop back into logging
    log.info("Alignment trimming ends")
    # write the output files
    write_gblocks_alignments_to_outdir(log, args.output, alignments, args.output_format)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))

if __name__ == '__main__':
    main()
