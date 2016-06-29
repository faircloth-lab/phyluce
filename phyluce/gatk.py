#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 17:17 PDT (-0700)
"""


import os
import re
import gzip
import glob
import numpy
import subprocess
from collections import OrderedDict

from phyluce.pth import get_user_param, get_user_path

from Bio import SeqIO


JAVA = get_user_param("java", "executable")
JAVA_PARAMS = get_user_param("java", "mem")
JAR_PATH = get_user_path("java", "jar")
GATK = get_user_param("java", "gatk")


def coverage(log, sample, assembly_pth, assembly, cores, bam):
    log.info("Computing coverage with GATK for {}".format(sample))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(assembly_pth)
    cmd = [
        JAVA,
        JAVA_PARAMS,
        "-jar",
        os.path.join(JAR_PATH, GATK),
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


def compute_coverage_metrics(contig_depth, trim=False):
    metadata = {
        "beginning-length":None,
        "beginning-mean-cov": None,
        "contig-start": None,
        "contig-end": None,
        "trim-start":None,
        "trim-end":None,
        "ending-length":None,
        "ending-mean-cov":None
    }
    depth_array = numpy.array(contig_depth)
    metadata["beginning-length"] = len(depth_array)
    metadata["beginning-mean-cov"] = numpy.round(numpy.mean(depth_array), 2)
    metadata["contig-start"] = 0
    metadata["contig-end"] = len(depth_array)
    if trim:
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
    else:
        metadata["trim-start"] = None
        metadata["trim-end"] = None
        metadata["ending-length"] = len(depth_array)
        metadata["ending-mean-cov"] = numpy.round(numpy.mean(depth_array), 2)
    return metadata


def get_trimmed_coverage_from_output(log, sample, assembly_pth, coverage, assembler):
    log.info("Screening and filtering contigs for coverage (3x ends, 5x avg.)")
    if assembler == "trinity":
        regex = re.compile("({}).*:(\d+)".format(get_user_param("headers", "trinity")))
    elif assembler == "velvet":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "velvet")))
    elif assembler == "abyss":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "abyss")))
    elif assembler == "idba":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "idba")))
    # setup starting values
    previous_match = None
    contig_depth = []
    contig_data = OrderedDict()
    overall_coverage = []
    overall_length = []
    overall_count = 1
    overall_contigs = {}
    pbc = os.path.join(
        assembly_pth,
        '{}-TRIMMED-per-base-coverage.txt.gz'.format(sample)
    )
    pcc = os.path.join(
        assembly_pth,
        '{}-TRIMMED-per-contig-coverage.txt'.format(sample)
    )
    upcc = os.path.join(
        assembly_pth,
        '{}-UNTRIMMED-per-contig-coverage.txt'.format(sample)
    )
    with open(coverage, 'rU') as infile:
        with gzip.open(pbc, 'w') as per_base_cov:
            with open(pcc, 'w') as per_contig_cov:
                with open(upcc, 'w') as unt_per_contig_cov:
                    # read header line
                    gatk_header = infile.readline()
                    # write headers to outfiles
                    per_contig_cov.write("name\tbeginning-length\tbeginning-mean-cov\ttrim-start\ttrim-end\tend-length\tend-mean-cov\n")
                    unt_per_contig_cov.write("name\tbeginning-length\tbeginning-mean-cov\n")
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
                            metadata = compute_coverage_metrics(contig_depth, trim=True)
                            unt_per_contig_cov.write("{}\t{}\t{}\n".format(
                                    previous_match,
                                    metadata["beginning-length"],
                                    metadata["beginning-mean-cov"]
                                ))
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


def get_untrimmed_coverage_from_output(log, sample, assembly_pth, coverage, assembler):
    log.info("Screening contigs for coverage")
    if assembler == "trinity":
        regex = re.compile("({}).*:(\d+)".format(get_user_param("headers", "trinity")))
    elif assembler == "velvet":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "velvet")))
    elif assembler == "abyss":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "abyss")))
    elif assembler == "idba":
        regex = re.compile("({}.*):(\d+)".format(get_user_param("headers", "idba")))
    # setup starting values
    previous_match = None
    contig_depth = []
    contig_data = OrderedDict()
    overall_coverage = []
    overall_length = []
    overall_count = 1
    overall_contigs = {}
    upcc = os.path.join(
        assembly_pth,
        '{}-UNTRIMMED-per-contig-coverage.txt'.format(sample)
    )
    with open(coverage, 'rU') as infile:
        with open(upcc, 'w') as unt_per_contig_cov:
            # read header line
            gatk_header = infile.readline()
            # write headers to outfiles
            unt_per_contig_cov.write("name\tbeginning-length\tbeginning-mean-cov\n")
            for line in infile:
                ls = line.split()
                search = regex.search(ls[0])
                match_name, pos = search.groups()
                if previous_match is None or match_name == previous_match:
                    # hold onto current match_name
                    previous_match = match_name
                    # compute metrics on current position
                    #contig_data[int(pos)] = line
                    contig_depth.append(int(ls[1]))
                elif match_name != previous_match:
                    metadata = compute_coverage_metrics(contig_depth, trim=False)
                    unt_per_contig_cov.write("{}\t{}\t{}\n".format(
                            previous_match,
                            metadata["beginning-length"],
                            metadata["beginning-mean-cov"]
                        ))
                    overall_contigs[previous_match] = metadata
                    overall_count += 1
                    overall_coverage.append(metadata["beginning-mean-cov"])
                    overall_length.append(metadata["beginning-length"])
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


def remove_coverage_files(log, assembly_pth, coverage):
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
