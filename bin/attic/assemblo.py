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
from phyluce.third_party import which

import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Batch assembly folders of reads
            using velvet and VelvetOptimiser""", )
    parser.add_argument('input',
        type=is_dir,
        action=FullPaths,
        help='The directory containing species-specific data')
    parser.add_argument('--output',
        action=FullPaths,
        help='The directory to hold the assembled contigs')
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
    parser.add_argument('--separate-reads',
        action="store_true",
        default=False,
        help='Choose this flag if you are assembling paired-end data in separate read files'
        )
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--exclude',
        type=str,
        default=[],
        nargs='+',
        help='Directories in <input> to exclude')
    group.add_argument('--include',
        type=str,
        default=[],
        nargs='+',
        help='Directories in <input> to include')
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


def assemble_paired_end_reads(args, velvet_optimiser, indiv, read, interleaved_dir):
    assert is_dir(interleaved_dir), IOError("Directory {} containing interleaved reads does not exist".format(interleaved_dir))
    singletons = os.path.join(
            interleaved_dir,
            "{0}{1}".format(indiv, '-READ-singleton.fastq.gz')
        )
    interleaved = os.path.join(
            interleaved_dir,
            "{0}{1}".format(indiv, '-READ1and2-interleaved.fastq.gz')
        )
    for f in [singletons, interleaved]:
        assert os.path.isfile(f), IOError("Missing sequence file(s): {0}".format(f))
    assembly_dir = os.path.join(args.output, indiv, 'assembly')
    mkdir_p(assembly_dir)
    os.chdir(assembly_dir)
    velveth = "-fastq.gz -shortPaired {} -short {}".format(singletons, interleaved)
    cmd = [velvet_optimiser,
            "-s", args.s,
            "-e", args.e,
            "-c", "'ncon'",
            "-t", args.t,
            "-a",
            "-f", velveth
        ]
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout, stderr, assembly_dir


def assemble_separate_paired_end_reads(args, velvet_optimiser, indiv, read, paired_dir):
    assert is_dir(paired_dir), IOError("Directory {} containing separate reads does not exist".format(paired_dir))
    singletons = os.path.join(
            paired_dir,
            "{0}{1}".format(indiv, '-READ-singleton.fastq.gz')
        )
    read1 = os.path.join(
            paired_dir,
            "{0}{1}".format(indiv, '-READ1.fastq.gz')
        )
    read2 = os.path.join(
            paired_dir,
            "{0}{1}".format(indiv, '-READ2.fastq.gz')
        )
    for f in [singletons, read1, read2]:
        assert os.path.isfile(f), IOError("Missing sequence file(s): {0}".format(f))
    assembly_dir = os.path.join(args.output, indiv, 'assembly')
    mkdir_p(assembly_dir)
    os.chdir(assembly_dir)
    velveth = "-fastq.gz -separate -shortPaired {} {} -short {}".format(read1, read2, singletons)
    cmd = [velvet_optimiser,
            "-s", args.s,
            "-e", args.e,
            "-c", "'ncon'",
            "-t", args.t,
            "-a",
            "-f", velveth
        ]
    #pdb.set_trace()
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout, stderr, assembly_dir


def assemble_single_end_reads(args, velvet_optimiser, indiv, read, singles_dir):
    assert is_dir(singles_dir), IOError("Directory {} containing single reads does not exist".format(singles_dir))
    singles = os.path.join(
            singles_dir,
            "{0}{1}".format(indiv, '-READS-interleaved.fastq.gz')
        )
    for f in [singles]:
        assert os.path.isfile(f), IOError("Missing sequence file(s): {0}".format(f))
    assembly_dir = os.path.join(args.output, indiv, 'assembly')
    mkdir_p(assembly_dir)
    os.chdir(assembly_dir)
    velveth = "-fastq.gz -short {}".format(singles)
    cmd = [velvet_optimiser,
            "-s", args.s,
            "-e", args.e,
            "-c", "'ncon'",
            "-t", args.t,
            "-a",
            "-f", velveth
        ]
    stdout, stderr = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return stdout, stderr, assembly_dir


def get_samples_to_run(args, reads):
    """docstring for get_samples_to_run"""
    all_names = [os.path.basename(read) for read in reads]
    if args.exclude:
        return set([name for name in all_names if name not in args.exclude])
    elif args.include:
        return set([name for name in all_names if name in args.include])
    else:
        return all_names
        
def get_velvet_optimiser(name='VelvetOptimiser'):
    """ensure that velvetg, velveth, VelvetOpt, and VelvetOptimiser are in $PATH"""
    # ensure velvetg and velveth are in $PATH
    velvetg = which("velvetg")
    velveth = which("velveth")
    velvet_opt = which("VelvetOpt")
    # we need velvetoptimiser - ensure that is in $PATH and return
    try:
        velvet_optimiser = which("{}".format(name))[0]
        return velvet_optimiser
    except EnvironmentError, e:
        velvet_optimiser = which("{}.pl".format(name))[0]
        return velvet_optimiser

def main():
    args = get_args()
    basedir = args.input
    if not args.output:
        args.output = args.input
    contig_dir = mkdir_p(os.path.join(args.output, 'contigs'))
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
    taxa_to_run = get_samples_to_run(args, reads)
    # find VelvetOptimiser
    velvet_optimiser = get_velvet_optimiser()
    for read in reads:
        indiv = os.path.basename(read)
        if indiv in taxa_to_run:
            if not args.separate_reads:
                read_dir = os.path.join(read, 'interleaved-adapter-quality-trimmed')
            else:
                read_dir = os.path.join(read, 'split-adapter-quality-trimmed')
            if not args.single_end and not args.separate_reads:
                stdout, stderr, assembly_dir = assemble_paired_end_reads(args, velvet_optimiser, indiv, read, read_dir)
            elif not args.single_end and args.separate_reads:
                stdout, stderr, assembly_dir = assemble_separate_paired_end_reads(args, velvet_optimiser, indiv, read, read_dir)
            else:
                stdout, stderr, assembly_dir = assemble_single_end_reads(args, velvet_optimiser, indiv, read, read_dir)
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
