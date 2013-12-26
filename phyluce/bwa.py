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
import re
import gzip
import glob
import numpy
import subprocess
from collections import OrderedDict
from Bio import SeqIO

import pdb

JAVA = "java"
JAVA_PARAMS = "-Xmx20g"
JAR_PATH = "/home/bcf/bin/"

def bwa_create_index_files(log, reference):
    log.info("Running bwa indexing against {}".format(reference))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(reference))
    cmd = ["bwa", "index", reference]
    with open('bwa-index-file.log', 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)

def bwa_create_sai(log, sample, sample_dir, ref, cores, reads, read):
    log.info("Creating read index file for {}".format(reads.file))
    cmd = [
        "bwa",
        "aln",
        "-t",
        str(cores),
        ref,
        reads.pth
    ]
    aln_out_fname = os.path.join(sample_dir, '{}-r{}.sai'.format(sample, read))
    aln_err_fname = os.path.join(sample_dir, '{}-r{}.bwa-aln-out.log'.format(sample, read))
    with open(aln_out_fname, 'w') as aln_out:
        with open(aln_err_fname, 'w') as aln_err:
            proc = subprocess.Popen(cmd, stdout=aln_out, stderr=aln_err)
            proc.communicate()
    return aln_out_fname


#def bwa_create_sai(log, sample, sample_dir, ref, cores, reads, read):
#    if read == 1:
#        aln_out_fname = bwa_run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
#    elif read == 2:
#        aln_out_fname = bwa_run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
#    return aln_out_fname


def bwa_se_align(log, sample, sample_dir, ref, cores, rS):
    rSsai = bwa_create_sai(log, sample, sample_dir, ref, cores, rS, 'S')
    cmd1 = [
        "bwa",
        "samse",
        ref,
        rSsai,
        rS.pth,
    ]
    cmd2 = [
        "samtools",
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, '{}.se.bwa-samse-out.log'.format(sample))
    samtools_out_fname = os.path.join(sample_dir, '{}.se.samtools-se-out.log'.format(sample))
    bam_out_fname = os.path.join(sample_dir, '{}-se.bam'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    # remove the sai files (they'll be stale soon)
    os.remove(rSsai)
    return bam_out_fname


def bwa_pe_align(log, sample, sample_dir, ref, cores, r1, r2):
    r1sai = bwa_create_sai(log, sample, sample_dir, ref, cores, r1, 1)
    r2sai = bwa_create_sai(log, sample, sample_dir, ref, cores, r2, 2)
    cmd1 = [
        "bwa",
        "sampe",
        "-a",
        "700",
        ref,
        r1sai,
        r2sai,
        r1.pth,
        r2.pth
    ]
    cmd2 = [
        "samtools",
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, '{}.pe.bwa-sampe-out.log'.format(sample))
    samtools_out_fname = os.path.join(sample_dir, '{}.pe.samtools-out.log'.format(sample))
    bam_out_fname = os.path.join(sample_dir, '{}.bam'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    # remove the sai files (they'll be stale soon)
    os.remove(r1sai)
    os.remove(r2sai)
    return bam_out_fname


def bwa_mem_se_align(log, sample, sample_dir, ref, cores, rS):
    #pdb.set_trace()
    cmd1 = [
        "bwa",
        "mem",
        "-t",
        str(cores),
        "-M",
        ref,
        rS.pth,
    ]
    cmd2 = [
        "samtools",
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, '{}.se.bwa-sampe-out.log'.format(sample))
    samtools_out_fname = os.path.join(sample_dir, '{}.se.samtools-view-out.log'.format(sample))
    bam_out_fname = os.path.join(sample_dir, '{}-se.bam'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    return bam_out_fname


def bwa_mem_pe_align(log, sample, sample_dir, ref, cores, r1, r2):
    #pdb.set_trace()
    cmd1 = [
        "bwa",
        "mem",
        "-t",
        str(cores),
        "-M",
        ref,
        r1.pth,
        r2.pth
    ]
    cmd2 = [
        "samtools",
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, '{}.pe.bwa-sampe-out.log'.format(sample))
    samtools_out_fname = os.path.join(sample_dir, '{}.pe.samtools-view-out.log'.format(sample))
    bam_out_fname = os.path.join(sample_dir, '{}.bam'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    return bam_out_fname



def new_bam_name(bam, append):
    pth, bamfname = os.path.split(bam)
    bamfname = os.path.splitext(bamfname)[0]
    new_bamfname = "{}-{}.bam".format(bamfname, append)
    new_bam = os.path.join(pth, new_bamfname)
    return new_bam


def picard_create_reference_dict(log, sample, sample_dir, reference):
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


def picard_clean_up_bam(log, sample, sample_dir, bam, type):
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


def picard_add_rg_header_info(log, sample, sample_dir, flowcell, bam, type):
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


def picard_merge_two_bams(log, sample, sample_dir, bam, bam_se):
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


def picard_mark_and_remove_dupes(log, sample, sample_dir, bam, type):
    log.info("Removing read duplicates from BAM for {}".format(sample))
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


def picard_calculate_hs_metrics(log, sample, sample_dir, reference, bam, target, bait):
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


def picard_get_percent_reads_on_target(log, hs_metrics_file, sample):
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


def samtools_index(log, sample, sample_dir, bam):
    log.info("Indexing BAM for {}".format(sample))
    cmd = [
        "samtools",
        "index",
        bam
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-index-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def samtools_create_faidx(log, sample, sample_dir, fasta):
    log.info("Indexing fasta for {}".format(sample))
    cmd = [
        "samtools",
        "faidx",
        fasta
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-faidx-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def gatk_coverage(log, sample, assembly_pth, assembly, cores, bam):
    log.info("Computing coverage with GATK for {}".format(sample))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(assembly_pth)
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, "GenomeAnalysisTK.jar"),
        "-T",
        "DepthOfCoverage",
        "-R",
        assembly,
        "-I",
        bam,
        "-o",
        "{}-coverage".format(sample),
        "-nt",
        str(cores),
        "--omitIntervalStatistics",
        "--omitLocusTable"
    ]
    gatk_coverage_fname = os.path.join(assembly_pth, '{}.GATK-coverage-out.log'.format(sample))
    with open(gatk_coverage_fname, 'w') as gatk_out:
        proc = subprocess.Popen(cmd, stdout=gatk_out, stderr=subprocess.STDOUT)
        proc.communicate()
    os.chdir(cwd)
    return os.path.join(assembly_pth, "{}-coverage".format(sample))


def compute_gatk_coverage_metrics(contig_depth):
    metadata = {
        "beginning-length":None,
        "beginning-mean-cov": None,
        "trim-start":None,
        "trim-end":None,
        "ending-length":None,
        "ending-mean-cov":None
    }
    depth_array = numpy.array(contig_depth)
    metadata["beginning-length"] = len(depth_array)
    metadata["beginning-mean-cov"] = numpy.round(numpy.mean(depth_array), 2)
    # find where coverage is >= 3
    min_cov = depth_array >= 3
    # get positions at edges where True
    true_positions = numpy.where(min_cov==True)
    try:
        metadata["trim-start"] = true_positions[0][0]
        metadata["trim-end"] = true_positions[0][-1] + 1
        good = depth_array[metadata["trim-start"]:metadata["trim-end"]]
    except IndexError as err:
        if err.message == 'index out of bounds':
            good = None
        else:
            raise err
    if good is not None:
        metadata["ending-length"] = len(good)
        metadata["ending-mean-cov"] = numpy.round(numpy.mean(good), 2)
    else:
        metadata["ending-length"] = None
        metadata["ending-mean-cov"] = None
    return metadata


def get_coverage_from_gatk(log, sample, assembly_pth, coverage, velvet):
    log.info("Screening and filtering contigs for coverage (3x ends, 5x avg.)")
    if not velvet:
        regex = re.compile("(comp\d+_c\d+_seq\d+).*:(\d+)")
    else:
        regex = re.compile("(NODE_\d+_length_\d+_cov_\d+.*):(\d+)")
    # setup starting values
    previous_match = None
    contig_depth = []
    contig_data = OrderedDict()
    overall_coverage = []
    overall_length = []
    overall_count = 0
    overall_contigs = {}
    pbc = os.path.join(
        assembly_pth,
        '{}-TRIMMED-per-base-coverage.txt.gz'.format(sample)
    )
    pcc = os.path.join(
        assembly_pth,
        '{}-TRIMMED-per-contig-coverage.txt'.format(sample)
    )
    with open(coverage, 'rU') as infile:
        with gzip.open(pbc, 'w') as per_base_cov:
            with open(pcc, 'w') as per_contig_cov:
                # read header line
                gatk_header = infile.readline()
                # write headers to outfiles
                per_contig_cov.write("name\tbeginning-length\tbeginning-mean-cov\ttrim-start\ttrim-end\tend-length\tend-mean-cov\n")
                per_base_cov.write(gatk_header)
                for line in infile:
                    ls = line.split()
                    search = regex.search(ls[0])
                    match_name, pos = search.groups()
                    if previous_match is None or match_name == previous_match:
                        # hold onto current match_name
                        previous_match = match_name
                        # compute metrics on current position
                        contig_data[int(pos)] = line
                        contig_depth.append(int(ls[1]))
                    elif match_name != previous_match:
                        metadata = compute_gatk_coverage_metrics(contig_depth)
                        if metadata["ending-mean-cov"] >= 5.0:
                            per_contig_cov.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                previous_match,
                                metadata["beginning-length"],
                                metadata["beginning-mean-cov"],
                                metadata["trim-start"],
                                metadata["trim-end"],
                                metadata["ending-length"],
                                metadata["ending-mean-cov"]
                            ))
                            for pos, line in contig_data.iteritems():
                                if pos-1 >= metadata["trim-start"] and pos-1 < metadata["trim-end"]:
                                    per_base_cov.write(line)
                            overall_contigs[previous_match] = metadata
                            overall_count += 1
                            overall_coverage.append(metadata["ending-mean-cov"])
                            overall_length.append(metadata["ending-length"])
                        # reset previous match to current
                        previous_match = match_name
                        # reset containers
                        contig_depth = []
                        contig_data = OrderedDict()
                        # compute metrics on current first position
                        contig_data[int(pos)] = line
                        contig_depth.append(int(ls[1]))
    log.info("\t{} contigs, mean coverage = {:.1f}, mean length = {:.1f}".format(
        overall_count,
        numpy.mean(overall_coverage),
        numpy.mean(overall_length)
    ))
    return overall_contigs


def remove_gatk_coverage_files(log, assembly_pth, coverage):
    log.info("Removing GATK coverage files we created (screened files remain)")
    for file in glob.glob("{}*".format(coverage)):
        # gzip coverage file - remove rest
        if file == coverage and not os.path.splitext(coverage)[1] == ".gz":
            log.info("[GZIP] {} (slow)".format(os.path.basename(file)))
            with open(file, "rb") as unzip:
                with gzip.open("{}.gz".format(file), "wb") as zip:
                    zip.writelines(unzip)
            # remove the unzipped file
            os.remove(file)
        else:
            # remove all other files
            log.info("[Remove] {}".format(os.path.basename(file)))
            os.remove(file)


def filter_screened_contigs_from_assembly(log, sample, assembly_pth, assembly, overall_contigs):
    log.info("Filtering assembly FASTAs for coverage")
    keep = set(overall_contigs.keys())
    outf_name = os.path.join(
        assembly_pth,
        "{}-trinity-FILTERED.fasta".format(sample)
    )
    with open(outf_name, 'w') as outf:
        with open(assembly, 'rU') as infile:
            for seq in SeqIO.parse(infile, 'fasta'):
                if seq.id in keep:
                    # get trimming coordinates
                    start = overall_contigs[seq.id]["trim-start"]
                    end = overall_contigs[seq.id]["trim-end"]
                    seq = seq[start:end]
                    seq.description = "len={}".format(len(seq))
                    outf.write(seq.format('fasta'))
    return outf_name



