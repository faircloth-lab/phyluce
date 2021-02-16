#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 06 July 2018 15:37 CDT (-0500)
"""

import os
import glob
import shutil
import subprocess
import configparser

# from phyluce.tests.common import get_contig_lengths_and_counts

import pytest
from Bio import SeqIO

import pdb


@pytest.fixture(autouse=True)
def cleanup_files(request):
    """cleanup extraneous log files"""

    def clean():
        log_files = os.path.join(
            request.config.rootdir, "phyluce", "tests", "*.log"
        )
        for file in glob.glob(log_files):
            os.remove(file)

    request.addfinalizer(clean)


@pytest.fixture(scope="module")
def o_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-observed"
    )
    os.mkdir(directory)

    def clean():
        shutil.rmtree(directory)

    request.addfinalizer(clean)
    return directory


@pytest.fixture(scope="module")
def e_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-expected"
    )
    return directory


@pytest.fixture(scope="module")
def conf_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-conf"
    )
    return directory


@pytest.fixture(scope="module")
def raw_dir(request):
    directory = os.path.join(
        request.config.rootdir,
        "phyluce",
        "tests",
        "test-expected",
        "raw-reads-short",
        "alligator-mississippiensis",
    )
    return directory


def get_match_count_cmd(
    o_dir, e_dir, conf_dir, output_config, request, incomplete=False
):
    program = "bin/assembly/phyluce_assembly_get_match_counts"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--locus-db",
        os.path.join(e_dir, "probe-match", "probe.matches.sqlite"),
        "--taxon-list-config",
        os.path.join(conf_dir, "taxon-set.conf"),
        "--taxon-group",
        "all",
        "--output",
        output_config,
        "--log-path",
        o_dir,
    ]
    if not incomplete:
        return cmd
    else:
        cmd.append("--incomplete-matrix")
        return cmd


def test_get_fastq_lengths(o_dir, e_dir, raw_dir, request):
    program = "bin/assembly/phyluce_assembly_get_fastq_lengths"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--input",
        raw_dir,
        "--csv",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    stdout_str = stdout.decode("utf-8")
    stdout_str_split = stdout_str.strip().split(",")[1:]
    expected = "7404,677024,91.44030253916802,0.1993821226016458,40,100,100.0"
    expected_str_split = expected.split(",")
    assert stdout_str_split == expected_str_split


def test_get_match_counts_complete(o_dir, e_dir, conf_dir, request):
    output_config = os.path.join(o_dir, "taxon-set.conf")
    cmd = get_match_count_cmd(o_dir, e_dir, conf_dir, output_config, request)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    obs_config = configparser.RawConfigParser(allow_no_value=True)
    obs_config.optionxform = str
    obs_config.read(output_config)
    expected_config = os.path.join(e_dir, "taxon-set.complete.conf")
    exp_config = configparser.RawConfigParser(allow_no_value=True)
    exp_config.optionxform = str
    exp_config.read(expected_config)
    assert obs_config == exp_config


def test_get_match_counts_incomplete(o_dir, e_dir, conf_dir, request):
    output_config = os.path.join(o_dir, "taxon-set.conf")
    cmd = get_match_count_cmd(
        o_dir, e_dir, conf_dir, output_config, request, incomplete=True
    )
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    obs_config = configparser.RawConfigParser(allow_no_value=True)
    obs_config.optionxform = str
    obs_config.read(output_config)
    expected_config = os.path.join(e_dir, "taxon-set.incomplete.conf")
    exp_config = configparser.RawConfigParser(allow_no_value=True)
    exp_config.optionxform = str
    exp_config.read(expected_config)
    assert obs_config == exp_config


def get_fastas_cmd(o_dir, e_dir, o_file, request, incomplete=False):
    program = "bin/assembly/phyluce_assembly_get_fastas_from_match_counts"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--locus-db",
        os.path.join(e_dir, "probe-match", "probe.matches.sqlite"),
        "--contigs",
        os.path.join(e_dir, "spades", "contigs"),
        "--locus-db",
        os.path.join(e_dir, "probe-match", "probe.matches.sqlite"),
        "--match-count-output",
        os.path.join(e_dir, "taxon-set.complete.conf"),
    ]
    if not incomplete:
        cmd.extend(["--output", o_file, "--log-path", o_dir])
    else:
        cmd.extend(
            [
                "--output",
                o_file,
                "--log-path",
                o_dir,
                "--incomplete-matrix",
                os.path.join(o_dir, "taxon-set.incomplete"),
            ]
        )
    return cmd


def test_get_fastas_complete(o_dir, e_dir, request):
    o_file = os.path.join(o_dir, "taxon-set.complete.fasta")
    cmd = get_fastas_cmd(o_dir, e_dir, o_file, request, incomplete=False)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    # read in the new outfile
    observed_sequences = SeqIO.to_dict(SeqIO.parse(o_file, "fasta"))
    # read in the expected outfile
    expected_sequences = SeqIO.to_dict(
        SeqIO.parse(os.path.join(e_dir, "taxon-set.complete.fasta"), "fasta")
    )
    # compare
    for k, v in observed_sequences.items():
        assert v.seq == expected_sequences[k].seq


def test_get_fastas_incomplete(o_dir, e_dir, request):
    o_file = os.path.join(o_dir, "taxon-set.incomplete.fasta")
    cmd = get_fastas_cmd(o_dir, e_dir, o_file, request, incomplete=True)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    # read in the new outfile
    observed_sequences = SeqIO.to_dict(SeqIO.parse(o_file, "fasta"))
    # read in the expected outfile
    expected_sequences = SeqIO.to_dict(
        SeqIO.parse(os.path.join(e_dir, "taxon-set.incomplete.fasta"), "fasta")
    )
    # compare
    for k, v in observed_sequences.items():
        assert v.seq == expected_sequences[k].seq


def test_explode_get_fastas_file_by_taxon(o_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_explode_get_fastas_file"
    output = os.path.join(o_dir, "exploded-by-taxa")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--input",
        os.path.join(e_dir, "taxon-set.complete.fasta"),
        "--output",
        output,
        "--by-taxon",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    # iterate over observed and expected data and get file stats
    expected = os.path.join(e_dir, "exploded-by-taxa")
    for taxon in [
        "alligator-mississippiensis",
        "gallus-gallus",
        "peromyscus-maniculatus",
        "rana-sphenocephafa",
    ]:
        fname = "{}.unaligned.fasta".format(taxon)
        observed_sequences = SeqIO.to_dict(
            SeqIO.parse(os.path.join(output, fname), "fasta")
        )
        expected_sequences = SeqIO.to_dict(
            SeqIO.parse(os.path.join(expected, fname), "fasta")
        )
        for k, v in observed_sequences.items():
            # assert all are one taxon
            assert taxon.replace("-", "_") in k
            # assert that obs == expected
            assert v.seq == expected_sequences[k].seq


def test_explode_get_fastas_file_by_locus(o_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_explode_get_fastas_file"
    output = os.path.join(o_dir, "exploded-by-locus")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--input",
        os.path.join(e_dir, "taxon-set.complete.fasta"),
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    # iterate over observed and expected data and get file stats
    expected = os.path.join(e_dir, "exploded-by-locus")
    for locus in [
        "uce-1732",
        "uce-2120",
        "uce-3046",
        "uce-4179",
        "uce-553",
        "uce-7014",
    ]:
        fname = "{}.unaligned.fasta".format(locus)
        observed_sequences = SeqIO.to_dict(
            SeqIO.parse(os.path.join(output, fname), "fasta")
        )
        expected_sequences = SeqIO.to_dict(
            SeqIO.parse(os.path.join(expected, fname), "fasta")
        )
        for k, v in observed_sequences.items():
            # assert all are one taxon
            assert locus in k
            # assert that obs == expected
            assert v.seq == expected_sequences[k].seq
