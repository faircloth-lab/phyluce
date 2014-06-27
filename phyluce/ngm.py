#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 17:48 PDT (-0700)
"""


import os
import subprocess


def create_index_files(log, reference):
    log.info("Running NextGenMap indexing against {}".format(reference))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(reference))
    cmd = ["ng", "-r", reference]
    with open('ngm.index.log', 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)


def pe_align(log, sample, sample_dir, ref, cores, r1, r2):
    bam_out_fname = os.path.join(sample_dir, '{}.bam'.format(sample))
    cmd = [
        "ngm"
        "-r",
        ref,
        "-1",
        r1,
        "-2",
        r2,
        "-b",
        "-o",
        bam_out_fname,
        "-t",
        str(cores),
        "-p"
    ]
    ngmpe_out_fname = os.path.join(sample_dir, '{}.ngm.pe.log'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(ngmpe_out_fname, 'w') as sampe_out:
        proc = subprocess.Popen(cmd, stdout=sampe_out, stderr=sampe_out)
        proc.communicate()
    return bam_out_fname


def se_align(log, sample, sample_dir, ref, cores, r1, r2):
    bam_out_fname = os.path.join(sample_dir, '{}-se.bam'.format(sample))
    cmd = [
        "ngm"
        "-r",
        ref,
        "-1",
        r1,
        "-b",
        "-o",
        bam_out_fname,
        "-t",
        str(cores)
    ]
    ngmse_out_fname = os.path.join(sample_dir, '{}.ngm.se.log'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(ngmse_out_fname, 'w') as sampe_out:
        proc = subprocess.Popen(cmd, stdout=sampe_out, stderr=sampe_out)
        proc.communicate()
    return bam_out_fname
