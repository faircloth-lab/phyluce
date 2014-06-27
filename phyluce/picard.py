#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 17:13 PDT (-0700)
"""


import os
import subprocess


JAVA = "java"
JAVA_PARAMS = "-Xmx20g"
JAR_PATH = "/home/bcf/bin/"


def new_bam_name(bam, append):
    pth, bamfname = os.path.split(bam)
    bamfname = os.path.splitext(bamfname)[0]
    new_bamfname = "{}-{}.bam".format(bamfname, append)
    new_bam = os.path.join(pth, new_bamfname)
    return new_bam


def create_reference_dict(log, sample, sample_dir, reference):
    log.info("Creating FASTA dict for {}".format(sample))
    outf = os.path.splitext(reference)[0] + ".dict"
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "CreateSequenceDictionary.jar"),
        "R={}".format(reference),
        "O={}".format(outf)
    ]
    picard_ref_dict_fname = os.path.join(sample_dir, '{}.picard-reference-dict-out.log'.format(sample))
    with open(picard_ref_dict_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()


def clean_up_bam(log, sample, sample_dir, bam, type):
    log.info("Cleaning BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "CL")
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "CleanSam.jar"),
        "I={}".format(bam),
        "O={}".format(new_bam)
    ]
    picard_clean_out_fname = os.path.join(sample_dir, '{}.{}.picard-clean-out.log'.format(sample, type))
    with open(picard_clean_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


def fix_mate_information(log, sample, sample_dir, bam, type):
    log.info("Fixing mate information for {}".format(sample))
    new_bam = new_bam_name(bam, "CL")
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "FixMateInformation.jar"),
        "I={}".format(bam),
        "O={}".format(new_bam),
        "VALIDATION_STRINGENCY=SILENT"
    ]
    picard_clean_out_fname = os.path.join(sample_dir, '{}.{}.picard.fixmate.log'.format(sample, type))
    with open(picard_clean_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


def add_rg_header_info(log, sample, sample_dir, flowcell, bam, type):
    #pdb.set_trace()
    log.info("Adding RG header to BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "RG")
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "AddOrReplaceReadGroups.jar"),
        "I={}".format(bam),
        "O={}".format(new_bam),
        "SORT_ORDER=coordinate",
        "RGPL=illumina",
        "RGPU={}".format(flowcell),
        "RGLB=Lib1",
        "RGID={}".format(sample),
        "RGSM={}".format(sample),
        "VALIDATION_STRINGENCY=LENIENT"
    ]
    picard_rg_out_fname = os.path.join(sample_dir, '{}.{}.picard-RG-out.log'.format(sample, type))
    with open(picard_rg_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


def merge_two_bams(log, sample, sample_dir, bam, bam_se):
    log.info("Merging BAMs for {}".format(sample))
    new_bam = new_bam_name(bam, "M")
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "MergeSamFiles.jar"),
        "SO=coordinate",
        "AS=true",
        "I={}".format(bam),
        "I={}".format(bam_se),
        "O={}".format(new_bam),
        "VALIDATION_STRINGENCY=LENIENT",
    ]
    picard_merge_out_fname = os.path.join(sample_dir, '{}.picard-merge-out.log'.format(sample))
    with open(picard_merge_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    os.remove(bam_se)
    return new_bam


def mark_duplicates(log, sample, sample_dir, bam, type):
    log.info("Marking read duplicates from BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "MD")
    metricsfile = os.path.join(sample_dir, "{}.{}.picard-metricsfile.txt".format(sample, type))
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "MarkDuplicates.jar"),
        "I={}".format(bam),
        "O={}".format(new_bam),
        "METRICS_FILE={}".format(metricsfile),
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250",
        "ASSUME_SORTED=true",
        "VALIDATION_STRINGENCY=SILENT",
        "REMOVE_DUPLICATES=false",
    ]
    picard_dd_out_fname = os.path.join(sample_dir, '{}.{}.picard-MD-out.log'.format(sample, type))
    with open(picard_dd_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


def calculate_hs_metrics(log, sample, sample_dir, reference, bam, target, bait):
    log.info("Calculating coverage metrics for {}".format(sample))
    hs_metrics_file = os.path.join(sample_dir, "{}.reads-on-target.txt".format(sample))
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "CalculateHsMetrics.jar"),
        "I={}".format(bam),
        "O={}".format(hs_metrics_file),
        "REFERENCE_SEQUENCE={}".format(reference),
        "TARGET_INTERVALS={}".format(target),
        "BAIT_INTERVALS={}".format(bait)
    ]
    picard_hs_out_fname = os.path.join(sample_dir, '{}.picard-hs-metrics-out.log'.format(sample))
    with open(picard_hs_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    return hs_metrics_file


def get_percent_reads_on_target(log, hs_metrics_file, sample):
    lines = []
    with open(hs_metrics_file, "rU") as metrics:
        for line in metrics:
            if line.strip().startswith("#") or line.strip() == "":
                pass
            else:
                lines.append(line.strip("\n"))
    try:
        assert len(lines) == 2
    except:
        raise IOError("Picard metrics file for {} has more than two lines".format(sample))
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))
