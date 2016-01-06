#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 17:16 PDT (-0700)
"""

import os
import subprocess

from phyluce.pth import get_user_path

import pdb


def index(log, sample, sample_dir, bam):
    log.info("Indexing BAM for {}".format(sample))
    cmd = [
        get_user_path("samtools", "samtools"),
        "index",
        bam
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-index-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def create_faidx(log, sample, sample_dir, fasta):
    log.info("Indexing fasta for {}".format(sample))
    cmd = [
        get_user_path("samtools", "samtools"),
        "faidx",
        fasta
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-faidx-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def sort(log, sample, sample_dir, bam):
    log.info("Sorting BAM for {}".format(sample))
    out_prefix = "{}.sorted".format(os.path.splitext(bam)[0])
    cmd = [
        get_user_path("samtools", "samtools"),
        "sort",
        bam,
        out_prefix
    ]
    samtools_out_fname = '{}.samtools-sort-out.log'.format(sample_dir)
    with open(samtools_out_fname, 'a') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return "{}.bam".format(out_prefix)


def call(log, sample, sample_dir, reference, bam, phase=None):
    if phase is None:
        log.info("Creating REF/ALT allele FASTQ file --Unphased--")
    else:
        log.info("Creating REF/ALT allele FASTQ file {}".format(phase))
    cmd1 = [
        get_user_path("samtools", "samtools"),
        "mpileup",
        "-u",
        "-f",
        reference,
        bam
    ]
    cmd2 = [
        get_user_path("samtools", "bcftools"),
        "view",
        "-cg",
        "-"
    ]
    cmd3 = [
        get_user_path("samtools", "vcfutils"),
        "vcf2fq"
    ]
    mpileup_out_fname = "{}.samtools-mpileup-out.log".format(sample_dir)
    bcftools_out_fname = "{}.samtools-bcftools-out.log".format(sample_dir)
    vcfutils_out_fname = "{}.samtools-vcfutils-out.log".format(sample_dir)
    if phase is None:
        vcfutils_fastq_fname = "{}.fq".format(sample_dir)
    else:
        vcfutils_fastq_fname = "{}.{}.fq".format(sample_dir, phase)
    with open(mpileup_out_fname, 'w') as mpileup_out:
        with open(bcftools_out_fname, 'w') as bcftools_out:
            with open(vcfutils_out_fname, 'w') as vcfutils_out:
                with open(vcfutils_fastq_fname, 'w') as vcfutils_fastq:
                    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=mpileup_out)
                    proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=bcftools_out)
                    proc3 = subprocess.Popen(cmd3, stdin=proc2.stdout, stdout=vcfutils_fastq, stderr=vcfutils_out)
                    proc1.stdout.close()
                    proc2.stdout.close()
                    proc3.communicate()
    return vcfutils_fastq_fname

def phase(log, sample, sample_dir, bam):
    log.info("Phasing BAM file for {}".format(sample))
    cmd = [
        get_user_path("samtools", "samtools"),
        "phase",
        "-A",
        "-F",
        "-Q",
        "20",
        "-b",
        sample_dir,
        bam
    ]
    samtools_out_fname = '{}.samtools-phase-out.log'.format(sample_dir)
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return "{}.0.bam".format(sample_dir), "{}.1.bam".format(sample_dir)
