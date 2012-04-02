#!/usr/bin/env python
# encoding: utf-8
"""
File: run_many_lastzs.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 April 2012 15:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import time
import argparse
import tempfile
import subprocess
import ConfigParser

from bx.seq import twobit
from multiprocessing import Pool, cpu_count
from phyluce.helpers import is_file, is_dir, FullPaths

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Run many lastz queries against the genomes in conf file""")
    parser.add_argument(
            "conf",
            action=FullPaths,
            type=is_file,
            help="""Path to the configuration file"""
        )
    parser.add_argument(
            "query",
            action=FullPaths,
            type=is_file,
            help="""Path to the query file""",
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""Path to the output directory to hold results""",
        )
    parser.add_argument(
            '--nprocs',
            type=int,
            default=1,
            help='The number of processors to use'
        )
    parser.add_argument(
            '--coverage',
            type=float,
            default=83.,
            help="""The fraction of bases in the entire input sequence
            (target or query, whichever is shorter) that are
            included in the alignment block, expressed as a percentage"""
        )
    parser.add_argument(
            '--identity',
            type=float,
            default=92.5,
            help="""The fraction of aligned bases (excluding columns containing
                gaps or non-ACGT characters) that are matches, expressed as a
                percentage"""
        )
    args = parser.parse_args()
    if args.nprocs > cpu_count:
        args.nprocs = cpu_count - 1
    return args


def start():
    start = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y %H:%M:%S", time.localtime(start))
    return start


def stop(start):
    end = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y %H:%M:%S", time.localtime(end))
    print 'Time for execution: ', (end - start) / 60, 'minutes'


def align_query_to_genomes(work):
    chromo, params = work
    query, coverage, identity = params
    format = "general-:score,name1,strand1,zstart1,end1,length1," + \
            "name2,strand2,zstart2,end2,length2,diff,cigar,identity," + \
            "continuity,coverage"
    qry = "{}[nameparse=full]".format(query)
    fd, out = tempfile.mkstemp(suffix='.lastz')
    os.close(fd)
    cmd = [
            'lastz',
            chromo,
            qry,
            "--strand=both",
            "--seed=12of19",
            "--transition",
            "--nogfextend",
            "--nochain",
            "--gap=400,30",
            "--xdrop=910",
            "--ydrop=8370",
            "--hspthresh=3000",
            "--gappedthresh=3000",
            "--noentropy",
            "--coverage={}".format(coverage),
            "--identity={}".format(identity),
            "--output={}".format(out),
            "--format={}".format(format)
        ]
    process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    stdout, stderr = process.communicate()
    print "\t\t{}".format(out)
    return out, stderr


def chunk_scaffolds(genome, size=10000000):
    """Given a genome in many scaffolds, build temp files of `size` Mbp
    for easier querying"""
    print '\tChunking files into {0} bp...'.format(size)
    chromos = []
    tb = twobit.TwoBitFile(file(genome))
    # split target file into `options.size` (~10 Mbp) chunks
    fd, out = tempfile.mkstemp(suffix='.fasta')
    os.close(fd)
    temp = open(out, 'w')
    length = 0
    for seq in tb.keys():
        sequence = tb[seq][0:]
        length += len(sequence)
        # write it to the outfile
        temp.write('>{0}\n{1}\n'.format(seq, sequence))
        if length > size:
            temp.close()
            # put tempfile name on stack
            chromos.append(out + '[multiple]')
            # open a new temp file
            fd, out = tempfile.mkstemp(suffix='.fasta')
            os.close(fd)
            temp = open(out, 'w')
            # reset sequence length
            length = 0
    return chromos


def save_results_and_cleanup(outdir, name, results, chunks=None):
    outp = os.path.join(outdir, name) + ".lastz"
    out = open(outp, 'wb')
    print "\tWriting the results file {}".format(outp)
    #pdb.set_trace()
    for temp, stderr in results:
        print '\t%s' % temp
        # read the file
        out.write(open(temp, 'rb').read())
        # cleanup the lastz output files
        os.remove(temp)
    out.close()
    print '\tCleaning up the chunked files...'
    if chunks:
        # cleanup the chunked files
        [os.remove(f.strip('[multiple]')) for f in chunks]


def main():
    args = get_args()
    if args.nprocs > 1:
        pool = Pool(args.nprocs)
    # get and print start time
    begin_run = start()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    params = (args.query, args.coverage, args.identity)
    # get align types ("Chromos"/"Scaffolds")
    if conf.has_section("chromos"):
        for genome in conf.items("chromos"):
            name, f = genome
            print "{0}\nWorking on {1}\n{0}\n".format("=" * 30, name)
            chromos = [os.path.join(f, chromo)
                    for chromo in twobit.TwoBitFile(file(f)).keys()]
            work = [(chromo, params) for chromo in chromos]
            if args.nprocs > 1:
                results = pool.map(align_query_to_genomes, work)
            else:
                results = map(align_query_to_genomes, work)
            save_results_and_cleanup(args.output, name, results)
    if conf.has_section("scaffolds"):
        for genome in conf.items("scaffolds"):
            name, f = genome
            print "{0}\nWorking on {1}\n{0}\n".format("=" * 30, name)
            chunks = chunk_scaffolds(f)
            work = [(chunk, params) for chunk in chunks]
            if args.nprocs > 1:
                results = pool.map(align_query_to_genomes, work)
            else:
                results = map(align_query_to_genomes, work)
            save_results_and_cleanup(args.output, name, results, chunks)
    # get and print end time
    stop(begin_run)


if __name__ == '__main__':
    main()