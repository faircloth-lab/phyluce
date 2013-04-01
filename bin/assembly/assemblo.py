"""
File: assemblo.py
Author: Brant Faircloth

Created by Brant Faircloth on 7 March 2012 09:31 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.
"""

import os
import re
import glob
import errno
import argparse
import subprocess
from phyluce.helpers import is_dir, FullPaths

#import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Batch assembly folders of reads
            using velvet and VelvetOptimiser""", )
    parser.add_argument('input',
        type=is_dir,
        action=FullPaths,
        help='The directory containing species-specific data')
    parser.add_argument('s',
        type=str,
        help='The starting kmer value (-s in VelvetOptimiser)')
    parser.add_argument('e',
        type=str,
        help='The ending kmer value (-e in VelvetOptimiser')
    parser.add_argument('t',
        type=str,
        help='The number of parallel processes to run (-t in VelvetOptimiser)')
    parser.add_argument('--single-end',
        action="store_true",
        default=False,
        help='Choose this flag if you are assembling single-end reads'
        )
    parser.add_argument('--exclude',
        type=str,
        default=[],
        nargs='+',
        help='Directories in <input> to exclude')
    return parser.parse_args()


def get_values(regex, stdout):
    return regex.search(stdout).groups()[0]


def mkdir_p(path):
    """borrowed from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python"""
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise
    return path


def assemble_paired_end_reads(args, indiv, read, interleaved_dir):
    singletons = os.path.join(
            interleaved_dir,
            "{0}{1}".format(indiv, '-READ-singleton.fastq.gz')
        )
    interleaved = os.path.join(
            interleaved_dir,
            "{0}{1}".format(indiv, '-READ1and2-interleaved.fastq.gz')
        )
    assert is_dir(interleaved_dir), \
            "NO interleaved directory"
    for f in [singletons, interleaved]:
        assert os.path.isfile(f), "Missing sequence file(s): {0}".format(f)
    assembly_dir = os.path.join(read, 'assembly')
    mkdir_p(assembly_dir)
    os.chdir(assembly_dir)
    velveth = "-fastq.gz -short {} -shortPaired {}".format(singletons, interleaved)
    cmd = ["VelvetOptimiser",
            "-s", args.s,
            "-e", args.e,
            "-c", "'ncon'",
            "-t", args.t,
            "-a",
            "-f", velveth
        ]
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout, stderr, assembly_dir


def assemble_single_end_reads(args, indiv, read, interleaved_dir):
    interleaved = os.path.join(
            interleaved_dir,
            "{0}{1}".format(indiv, '-READS-interleaved.fastq.gz')
        )
    assert is_dir(interleaved_dir), \
            "NO interleaved directory"
    for f in [interleaved]:
        assert os.path.isfile(f), "Missing sequence file(s): {0}".format(f)
    assembly_dir = os.path.join(read, 'assembly')
    mkdir_p(assembly_dir)
    os.chdir(assembly_dir)
    velveth = "-fastq.gz -short {}".format(interleaved)
    cmd = ["VelvetOptimiser",
            "-s", args.s,
            "-e", args.e,
            "-c", "'ncon'",
            "-t", args.t,
            "-a",
            "-f", velveth
        ]
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout, stderr, assembly_dir

def main():
    args = get_args()
    basedir = args.input
    contig_dir = mkdir_p(os.path.join(basedir, 'contigs'))
    # Total number of contigs: 4882
    contig = re.compile('Total\snumber\sof\scontigs:\s(\d+)')
    kmer = re.compile('Velvet\shash\svalue:\s(\d+)')
    #n50: 335
    n50 = re.compile('n50:\s(\d+)')
    print "Running...\n"
    print "taxon,opt-kmer,contigs,n50,comments"
    # skip contigs dir that we just added
    reads = [f for f in glob.glob(os.path.join(args.input, '*')) \
            if os.path.basename(f) != 'contigs']
    for read in reads:
        indiv = os.path.basename(read)
        if indiv not in args.exclude:
            interleaved_dir = os.path.join(read, 'interleaved-adapter-quality-trimmed')
            if not args.single_end:
                stdout, stderr, assembly_dir = assemble_paired_end_reads(args, indiv, read, interleaved_dir)
            else:
                stdout, stderr, assembly_dir = assemble_single_end_reads(args, indiv, read, interleaved_dir)
            results = [indiv]
            for search in [kmer, contig, n50]:
                results.append(get_values(search, stderr))
            if results[1] == args.s or results[1] == args.e:
                results.extend(["kmer at edge of range"])
            print "{}".format(','.join(results))
            # symlink assembly contig file
            #pdb.set_trace()
            fasta = glob.glob(os.path.join(assembly_dir, 'auto_data_*/contigs.fa'))[0]
            fastalink = os.path.join(contig_dir, "{}.contigs.fasta".format(indiv))
            os.symlink(fasta, fastalink)
            # move up to basedir
            os.chdir(basedir)

if __name__ == '__main__':
    main()
