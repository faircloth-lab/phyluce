    #!/usr/bin/env python
# encoding: utf-8
"""
File: mpi_sate.py
Author: Brant Faircloth

Created by Brant Faircloth on 04 May 2012 15:05 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""


import os
import sys
import copy
import glob
import shutil
import argparse
import tempfile
import subprocess
from collections import defaultdict
from phyluce.helpers import is_dir, is_file, FullPaths
from seqtools.sequence import fasta

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "infile",
            action=FullPaths,
            help="""Help text"""
        )
    parser.add_argument(
            "output",
            type=is_dir,
            action=FullPaths,
            help="""Help text""",
        )
    parser.add_argument('species',
            type=int,
            default=None, \
            help='Number of species expected in each alignment.'
        )
    parser.add_argument(
            "sate",
            action=FullPaths,
            help="""Help text""",
        )
    parser.add_argument(
            "cfg",
            action=FullPaths,
            help="""Help text""",
        )
    parser.add_argument('--faircloth',
            action='store_true',
            default=False,
            help='Take faircloth+stephens probe names'
        )
    parser.add_argument('--ambiguous',
            action='store_true',
            default=False,
            help='Allow reads in alignments containing N-bases'
        )
    parser.add_argument('--notstrict',
            action='store_true',
            default=False,
            help='Allow alignments containing not all species'
        )
    parser.add_argument('--verbose',
            action='store_true',
            default=False,
            help='Give verbose output'
        )
    parser.add_argument(
            "--parallelism",
            choices=['mpi', 'multiprocessing', 'single'],
            default='single',
            help="""The type of parallelism to use.""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default='8',
            help="""The number of compute cores to use in multiprocessing.""",
        )
    return parser.parse_args()


def build_locus_dict(loci, locus, record, ambiguous=False):
    if not ambiguous:
        if not 'N' in record.sequence:
            loci[locus].append(record)
        else:
            sys.stdout.write('Skipping {0} because it contains ambiguous bases\n'.format(record.identifier))
            sys.stdout.flush()
    else:
        loci[locus].append(record)
    return loci


def get_fasta_dict(args):
    if args.verbose:
        sys.stdout.write('Building the locus dictionary...\n')
        if args.ambiguous:
            sys.stdout.write('NOT removing sequences with ambiguous bases...\n')
        else:
            sys.stdout.write('Removing ALL sequences with ambiguous bases...\n')
    sys.stdout.flush()
    loci = defaultdict(list)
    if os.path.isfile(args.infile):
        for record in fasta.FastaReader(args.infile):
            if not args.faircloth:
                locus = record.identifier.split('|')[1]
            else:
                locus = '_'.join([record.identifier.split('|')[0], \
                    record.identifier.split('|')[1].split('_')[0]])
            loci = build_locus_dict(loci, locus, record, args.ambiguous)
    # work with a directory of fastas if we have those - get locus name from
    # filename
    elif os.path.isdir(args.infile):
        for ff in glob.glob(os.path.join(args.infile, '*.fa*')):
            locus = os.path.splitext(os.path.basename(ff))[0]
            for record in fasta.FastaReader(ff):
                loci = build_locus_dict(loci, locus, record, args.ambiguous)
    # workon a copy so we can iterate and delete
    snapshot = copy.deepcopy(loci)
    # iterate over loci to check for all species at a locus
    for locus, data in snapshot.iteritems():
        if args.notstrict:
            if len(data) < 3:
                t = "\tDropping Locus {0} because it has fewer " + \
                        "than the minimum number " + \
                        "of taxa for alignment (N < 2)\n"
                sys.stdout.write((t).format(locus))
                sys.stdout.flush()
                del loci[locus]
        else:
            if len(data) < args.species:
                del loci[locus]
                t = "\tDropping Locus {0} because it has fewer " + \
                        "than the minimum number " + \
                        "of taxa for alignment (N < 2)\n"
                sys.stdout.write((t).format(locus))
                sys.stdout.flush()
    return loci


def worker(params):
    locus, opts = params
    name, sequences = locus
    sate, cfg = opts
    # create a tempdir to hold all our stuff
    working = tempfile.mkdtemp()
    # write content to outfile
    descriptor, path = tempfile.mkstemp(dir=working, suffix='.mpi.fasta')
    os.close(descriptor)
    tf = fasta.FastaWriter(path)
    [tf.write(seq) for seq in sequences]
    tf.close()
    # run SATe
    cli = [
            'python',
            sate,
            '--input',
            path,
            '--output-directory',
            working,
            '--temporaries',
            working,
            cfg
        ]
    stderr, stdout = subprocess.Popen(
            cli,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE
        ).communicate()
    # get contents of output file(s)
    aln_name = "satejob.marker001.{0}.aln".format(os.path.splitext(os.path.basename(path))[0])
    aln_file = os.path.join(working, aln_name)
    aln = open(aln_file, 'rU').read()
    # zap working tempdir
    shutil.rmtree(working)
    sys.stdout.write('.')
    sys.stdout.flush()
    # return filename and align so we can store resulting alignments reasonably
    return (name, aln)


def main():
    args = get_args()
    # iterate over files reading contents into a list that we'll pass to sate
    loci = get_fasta_dict(args)
    opts = [[args.sate, args.cfg] for i in range(len(loci))]
    params = zip(loci.items(), opts)
    sys.stdout.write('Aligning')
    # pass to map or Pool.map or mpimap
    alignments = mmap(worker, params)
    for data in alignments:
        name, aln = data
        out_file = os.path.join(args.output, name) + '.aln'
        out = open(out_file, 'w')
        out.write(aln)
        out.close()

if __name__ == '__main__':
    args = get_args()
    if args.parallelism == 'mpi':
        from deap.dtm import map as mmap
        from deap.dtm import start
        start(main)
    elif args.parallelism == 'multiprocessing':
        from multiprocessing import Pool
        pool = Pool(args.cores)
        mmap = pool.map
        main()
    elif args.parallelism == 'single':
        mmap = map
        main()
