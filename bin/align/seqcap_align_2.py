#!/usr/bin/env python
# encoding: utf-8
"""
File: seqcap_align_2.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 March 2012 11:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Parallel aligner for UCE fasta files generated with
assembly/get_fastas_from_match_counts.py

"""

import os
import sys
import copy
import shutil
import argparse
import tempfile
import multiprocessing
from collections import defaultdict

from seqtools.sequence import fasta

import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Align records in a file of UCE fastas"""
        )
    parser.add_argument('infile',
            help='The file containing fasta reads associated with UCE loci'
        )
    parser.add_argument('outdir',
            help='A directory for the output.'
        )
    parser.add_argument('species',
            type=int,
            default=None, \
            help='Number of species expected in each alignment.'
        )
    parser.add_argument('--aligner',
            choices=['dialign', 'muscle', 'mafft'],
            default='muscle',
            help='The aligner to use.'
        )
    parser.add_argument('--faircloth',
            action='store_true',
            default=False,
            help='Take faircloth+stephens probe names'
        )
    parser.add_argument('--notstrict',
            action='store_true',
            default=False,
            help='Allow alignments containing not all species'
        )
    parser.add_argument('--notrim',
            action='store_true',
            default=False,
            help='Do not trim alignments'
        )
    parser.add_argument('--window',
            type=int,
            default=20,
            help='Sliding window size for trimming'
        )
    parser.add_argument('--threshold',
            type=float,
            default=0.5,
            help='Threshold cutoff for trimming'
        )
    parser.add_argument('--proportion',
            type=float,
            default=0.3,
            help='Proportional removal of gaps'
        )
    parser.add_argument('--ambiguous',
            action='store_true',
            default=False,
            help='Allow reads in alignments containing N-bases'
        )
    parser.add_argument('--cores',
            type=int,
            default=1,
            help='Use multiple cores for alignment'
        )
    return parser.parse_args()


def build_locus_dict(loci, locus, record, ambiguous=False):
    if not ambiguous:
        if not 'N' in record.sequence:
            loci[locus].append(record)
        else:
            print 'Skipping {0} because it contains ambiguous bases'.format(record.identifier)
    else:
        loci[locus].append(record)
    return loci


def create_locus_specific_fasta(sequences):
    fd, fasta_file = tempfile.mkstemp(suffix='.fasta')
    os.close(fd)
    fasta_writer = fasta.FastaWriter(fasta_file)
    for seq in sequences:
        fasta_writer.write(seq)
    fasta_writer.close()
    return fasta_file


def align(params):
    locus, opts = params
    name, sequences = locus
    # get additional params from params tuple
    window, threshold, notrim, proportion = opts
    fasta = create_locus_specific_fasta(sequences)
    aln = Align(fasta)
    aln.run_alignment(consensus=False)
    if notrim:
        aln.trimmed_alignment = aln.alignment
    else:
        aln.trim_alignment(
                method='running',
                window_size=window,
                threshold=threshold,
                proportion=proportion
            )
    sys.stdout.write(".")
    sys.stdout.flush()
    return (name, aln)


def get_fasta_dict(args):
    print 'Building the locus dictionary...'
    if args.ambiguous:
        print 'NOT removing sequences with ambiguous bases...'
    else:
        print 'Removing ALL sequences with ambiguous bases...'
    loci = defaultdict(list)
    for record in fasta.FastaReader(args.infile):
        if not args.faircloth:
            locus = record.identifier.split('|')[1]
        else:
            locus = '_'.join([record.identifier.split('|')[0], \
                record.identifier.split('|')[1].split('_')[0]])
        loci = build_locus_dict(loci, locus, record, args.ambiguous)
    # workon a copy so we can iterate and delete
    snapshot = copy.deepcopy(loci)
    # iterate over loci to check for all species at a locus
    for locus, data in snapshot.iteritems():
        if args.notstrict:
            if len(data) < 3:
                t = "\tDropping Locus {0} because it has fewer " + \
                        "than the minimum number " + \
                        "of taxa for alignment (N < 2)"
                print(t).format(locus)
                del loci[locus]
        else:
            if len(data) < args.species:
                del loci[locus]
                t = "\tDropping Locus {0} because it has fewer " + \
                        "than the minimum number " + \
                        "of taxa for alignment (N < 2)"
                print(t).format(locus)
    return loci


def create_output_dir(outdir):
    print 'Creating output directory...'
    if os.path.exists(outdir):
        answer = raw_input("Output directory exists, remove [Y/n]? ")
        if answer == "Y":
            shutil.rmtree(outdir)
        else:
            sys.exit()
    os.makedirs(outdir)


def write_alignments_to_outdir(outdir, alignments, format='nexus'):
    formats = {
            'clustal': '.clw',
            'emboss': '.emboss',
            'fasta': '.fa',
            'nexus': '.nex',
            'phylip': '.phylip',
            'stockholm': '.stockholm'
        }
    print '\nWriting output files...'
    for tup in alignments:
        locus, aln = tup
        if aln.trimmed_alignment is not None:
            outname = "{}{}".format(os.path.join(outdir, locus), formats[format])
            outf = open(outname, 'w')
            outf.write(aln.trimmed_alignment.format(format))
            outf.close()
        else:
            print "Dropped {0} from output".format(locus)


def main(args):
    create_output_dir(args.outdir)
    loci = get_fasta_dict(args)
    sys.stdout.write("\nAligning with {}".format(str(args.aligner).upper()))
    sys.stdout.flush()
    opts = [[args.window, args.threshold, args.notrim, args.proportion] \
            for i in range(len(loci))]
    params = zip(loci.items(), opts)
    if args.cores:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        alignments = pool.map(align, params)
    else:
        alignments = map(align, params)
    write_alignments_to_outdir(args.outdir, alignments)


if __name__ == '__main__':
    args = get_args()
    # globally import Align method
    if args.aligner == 'muscle':
        from phyluce.muscle import Align
    elif args.aligner == 'mafft':
        from phyluce.mafft import Align
    elif args.aligner == 'dialign':
        from phyluce.dialign import Align
    main(args)
