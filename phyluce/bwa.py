#!/usr/bin/env python
# encoding: utf-8
"""
File: bwa.py
Author: Brant Faircloth

Created by Brant Faircloth on 30 September 2013 11:09 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""
import os
import subprocess

from phyluce.pth import get_user_param, get_user_path

# import pdb


def create_index_files(log, reference):
    log.info("Running bwa indexing against {}".format(reference))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(reference))
    cmd = [get_user_path("binaries", "bwa"), "index", reference]
    with open("bwa-index-file.log", "a") as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)


def create_sai(log, sample, sample_dir, ref, cores, reads, read):
    log.info("Creating read index file for {}".format(reads.file))
    cmd = [get_user_path("binaries", "bwa"), "aln", "-t", str(cores), ref, reads.pth]
    aln_out_fname = os.path.join(sample_dir, "{}-r{}.sai".format(sample, read))
    aln_err_fname = os.path.join(
        sample_dir, "{}-r{}.bwa-aln-out.log".format(sample, read)
    )
    with open(aln_out_fname, "w") as aln_out:
        with open(aln_err_fname, "w") as aln_err:
            proc = subprocess.Popen(cmd, stdout=aln_out, stderr=aln_err)
            proc.communicate()
    return aln_out_fname


# def bwa_create_sai(log, sample, sample_dir, ref, cores, reads, read):
#    if read == 1:
#        aln_out_fname = bwa_run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
#    elif read == 2:
#        aln_out_fname = bwa_run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
#    return aln_out_fname


def se_align(log, sample, sample_dir, ref, cores, rS):
    rSsai = create_sai(log, sample, sample_dir, ref, cores, rS, "S")
    cmd1 = [get_user_path("binaries", "bwa"), "samse", ref, rSsai, rS.pth]
    cmd2 = [get_user_path("binaries", "samtools"), "view", "-bS", "-"]
    sampe_out_fname = os.path.join(sample_dir, "{}.se.bwa-samse-out.log".format(sample))
    samtools_out_fname = os.path.join(
        sample_dir, "{}.se.samtools-se-out.log".format(sample)
    )
    bam_out_fname = os.path.join(sample_dir, "{}-se.bam".format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, "w") as sampe_out:
        with open(samtools_out_fname, "w") as samtools_out:
            with open(bam_out_fname, "w") as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(
                    cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out
                )
                proc1.stdout.close()
                proc2.communicate()
    # remove the sai files (they'll be stale soon)
    os.remove(rSsai)
    return bam_out_fname


def pe_align(log, sample, sample_dir, ref, cores, r1, r2):
    r1sai = create_sai(log, sample, sample_dir, ref, cores, r1, 1)
    r2sai = create_sai(log, sample, sample_dir, ref, cores, r2, 2)
    cmd1 = [
        get_user_path("binaries", "bwa"),
        "sampe",
        "-a",
        "700",
        ref,
        r1sai,
        r2sai,
        r1.pth,
        r2.pth,
    ]
    cmd2 = [get_user_path("binaries", "samtools"), "view", "-bS", "-"]
    sampe_out_fname = os.path.join(sample_dir, "{}.pe.bwa-sampe-out.log".format(sample))
    samtools_out_fname = os.path.join(
        sample_dir, "{}.pe.samtools-out.log".format(sample)
    )
    bam_out_fname = os.path.join(sample_dir, "{}.bam".format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, "w") as sampe_out:
        with open(samtools_out_fname, "w") as samtools_out:
            with open(bam_out_fname, "w") as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(
                    cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out
                )
                proc1.stdout.close()
                proc2.communicate()
    # remove the sai files (they'll be stale soon)
    os.remove(r1sai)
    os.remove(r2sai)
    return bam_out_fname


def mem_se_align(log, sample, sample_dir, ref, cores, rS):
    # pdb.set_trace()
    cmd1 = [
        get_user_path("binaries", "bwa"),
        "mem",
        "-t",
        str(cores),
        "-M",
        ref,
        rS.pth,
    ]
    cmd2 = [get_user_path("binaries", "samtools"), "view", "-bS", "-"]
    sampe_out_fname = os.path.join(sample_dir, "{}.se.bwa-sampe-out.log".format(sample))
    samtools_out_fname = os.path.join(
        sample_dir, "{}.se.samtools-view-out.log".format(sample)
    )
    bam_out_fname = os.path.join(sample_dir, "{}-se.bam".format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, "w") as sampe_out:
        with open(samtools_out_fname, "w") as samtools_out:
            with open(bam_out_fname, "w") as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(
                    cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out
                )
                proc1.stdout.close()
                proc2.communicate()
    return bam_out_fname


def mem_pe_align(log, sample, sample_dir, ref, cores, r1, r2):
    # pdb.set_trace()
    cmd1 = [
        get_user_path("binaries", "bwa"),
        "mem",
        "-t",
        str(cores),
        "-M",
        ref,
        r1.pth,
        r2.pth,
    ]
    cmd2 = [get_user_path("binaries", "samtools"), "view", "-bS", "-"]
    sampe_out_fname = os.path.join(sample_dir, "{}.pe.bwa-sampe-out.log".format(sample))
    samtools_out_fname = os.path.join(
        sample_dir, "{}.pe.samtools-view-out.log".format(sample)
    )
    bam_out_fname = os.path.join(sample_dir, "{}.bam".format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, "w") as sampe_out:
        with open(samtools_out_fname, "w") as samtools_out:
            with open(bam_out_fname, "w") as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(
                    cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out
                )
                proc1.stdout.close()
                proc2.communicate()
    return bam_out_fname
