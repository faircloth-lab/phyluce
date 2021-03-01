#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
(c) 2021 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 2021-02-14 T10:44:15-06:00
"""


import os
import re
import glob
import shutil
import platform
import subprocess
import logging

import pytest
from Bio import AlignIO
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


@pytest.mark.skipif(
    platform.processor() == "arm64", reason="Won't run on arm64"
)
def test_align_gblocks_trim(o_dir, e_dir, request):
    program = (
        "bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed"
    )
    output = os.path.join(o_dir, "mafft-gblocks")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-gblocks", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_trimal_trim(o_dir, e_dir, request):
    program = (
        "bin/align/phyluce_align_get_trimal_trimmed_alignments_from_untrimmed"
    )
    output = os.path.join(o_dir, "mafft-trimal")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-trimal", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_edge_trim(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_trimmed_alignments_from_untrimmed"
    output = os.path.join(o_dir, "mafft-edge-trim")
    # note that thus only uses alignemnts with an odd
    # number of taxa so ties in base composition at a
    # column do not cause random differences in expected output
    # this also completes testing of generic_align and seqalign
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-for-edge-trim"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        print(name)
        expected_file = os.path.join(e_dir, "mafft-edge-trim", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_missing_data_designators(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_add_missing_data_designators"
    output = os.path.join(o_dir, "mafft-missing-data-designators")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--match-count-output",
        os.path.join(e_dir, "taxon-set.incomplete.conf"),
        "--incomplete-matrix",
        os.path.join(e_dir, "taxon-set.incomplete"),
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        print(name)
        expected_file = os.path.join(
            e_dir, "mafft-missing-data-designators", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_convert_degen_bases(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_convert_degen_bases"
    output = os.path.join(o_dir, "mafft-degen-bases-converted")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-degen-bases"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        print(name)
        expected_file = os.path.join(
            e_dir, "mafft-degen-bases-converted", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_convert_align_mafft_fasta_to_nexus(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_convert_one_align_to_another"
    output = os.path.join(o_dir, "mafft-fasta-to-nexus")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-fasta-to-nexus", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_convert_align_mafft_fasta_to_phylip_relaxed(
    o_dir, e_dir, request
):
    program = "bin/align/phyluce_align_convert_one_align_to_another"
    output = os.path.join(o_dir, "mafft-fasta-to-phylip-relaxed")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "phylip-relaxed",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-fasta-to-phylip-relaxed", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_convert_align_mafft_nexus_to_fasta(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_convert_one_align_to_another"
    output = os.path.join(o_dir, "mafft-nexus-to-fasta")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-fasta-to-nexus"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--output-format",
        "fasta",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-nexus-to-fasta", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_explode_alignments_by_taxon(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_explode_alignments"
    output = os.path.join(o_dir, "mafft-exploded-by-taxon")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--by-taxon",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-exploded-by-taxon", name)
        observed = SeqIO.to_dict(SeqIO.parse(output_file, "fasta"))
        expected = SeqIO.to_dict(SeqIO.parse(expected_file, "fasta"))
        for name, observed in observed.items():
            assert expected[name].seq == observed.seq


def test_align_remove_locus_name(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_remove_locus_name_from_files"
    output = os.path.join(o_dir, "mafft-gblocks-clean")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-gblocks-clean", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_extract_taxa_from_alignments_exclude(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_extract_taxa_from_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-drop-gallus-gallus")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--exclude",
        "gallus_gallus",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-drop-gallus-gallus", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_extract_taxa_from_alignments_include(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_extract_taxa_from_alignments"
    output = os.path.join(
        o_dir, "mafft-gblocks-clean-keep-gallus-and-peromyscus"
    )
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--include",
        "gallus_gallus",
        "peromyscus_maniculatus",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    # pdb.set_trace()
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-keep-gallus-and-peromyscus", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_extract_taxon_fasta_from_alignments(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_extract_taxon_fasta_from_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-gallus.fasta")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--taxon",
        "gallus_gallus",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    assert output, "There are is no output"
    expected_file = os.path.join(e_dir, "mafft-gblocks-clean-gallus.fasta")
    observed = SeqIO.to_dict(SeqIO.parse(output, "fasta"))
    expected = SeqIO.to_dict(SeqIO.parse(expected_file, "fasta"))
    for name, observed in observed.items():
        assert expected[name].seq == observed.seq


def test_align_filter_alignments(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_filter_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-filtered-alignments")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--input-format",
        "nexus",
        "--containing-data-for",
        "gallus_gallus",
        "--min-length",
        "600",
        "--min-taxa",
        "3",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = [
        os.path.basename(i) for i in glob.glob(os.path.join(output, "*"))
    ]
    assert output_files, "There are no output files"
    expected_files = [
        os.path.basename(i)
        for i in glob.glob(
            os.path.join(e_dir, "mafft-gblocks-filtered-alignments", "*")
        )
    ]
    assert set(output_files) == set(expected_files)


def test_align_concatenate_to_nexus(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_concatenate_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-concat-nexus")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--nexus",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-concat-nexus", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_concatenate_to_phylip(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_concatenate_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-concat-phylip")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--output",
        output,
        "--phylip",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-concat-phylip", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_concatenate_fasta_to_phylip(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_concatenate_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-fasta-concat")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean-fasta"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--phylip",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-fasta-concat", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_get_align_summary_data(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_align_summary_data"
    output = os.path.join(o_dir, "gblocks-clean-align-summary.csv")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output-stats",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    assert output, "There are no output files"
    output_dict = {}
    with (open(output)) as output_file:
        for line in output_file:
            ls = line.strip().split(",")
            output_dict[ls[0]] = ",".join(ls[1:])
    expected = os.path.join(e_dir, "gblocks-clean-align-summary.csv")
    expected_dict = {}
    with (open(expected)) as expected_file:
        for line in expected_file:
            ls = line.strip().split(",")
            expected_dict[ls[0]] = ",".join(ls[1:])
    for k, v in output_dict.items():
        assert expected_dict[k] == v


def test_align_get_informative_sites(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_informative_sites"
    output = os.path.join(o_dir, "mafft-gblocks-clean-informative-sites.csv")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    assert output, "There are no output files"
    output_dict = {}
    with (open(output)) as output_file:
        for line in output_file:
            ls = line.strip().split(",")
            output_dict[ls[0]] = ",".join(ls[1:])
    expected = os.path.join(e_dir, "mafft-gblocks-clean-informative-sites.csv")
    expected_dict = {}
    with (open(expected)) as expected_file:
        for line in expected_file:
            ls = line.strip().split(",")
            expected_dict[ls[0]] = ",".join(ls[1:])
    for k, v in output_dict.items():
        assert expected_dict[k] == v


def test_align_get_only_loci_with_min_taxa(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_only_loci_with_min_taxa"
    output = os.path.join(o_dir, "mafft-gblocks-clean-75p")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output",
        output,
        "--taxa",
        "4",
        "--percent",
        "0.75",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = [
        os.path.basename(i) for i in glob.glob(os.path.join(output, "*"))
    ]
    expected_files = [
        os.path.basename(i)
        for i in glob.glob(os.path.join(e_dir, "mafft-gblocks-clean-75p", "*"))
    ]
    assert set(output_files) == set(expected_files)


def test_align_get_ry_recoded_alignments(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_ry_recoded_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-ry")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-gblocks-clean-ry", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_get_ry_recoded_alignments_binary(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_ry_recoded_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-ry-binary")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--output",
        output,
        "--binary",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-gblocks-clean-ry-binary", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_get_taxon_locus_counts_in_alignments(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_get_taxon_locus_counts_in_alignments"
    output = os.path.join(o_dir, "mafft-gblocks-clean-taxon-counts.csv")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean"),
        "--input-format",
        "nexus",
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    assert output, "There are no output files"
    expected_file = os.path.join(e_dir, "mafft-gblocks-clean-taxon-counts.csv")
    observed = open(output).read()
    expected = open(expected_file).read()
    assert observed == expected


def test_align_remove_empty_taxa(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_remove_empty_taxa"
    output = os.path.join(o_dir, "mafft-missing-data-designators-removed")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-missing-data-designators"),
        "--input-format",
        "nexus",
        "--output-format",
        "nexus",
        "--output",
        output,
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(
            e_dir, "mafft-missing-data-designators-removed", name
        )
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_align_screen_alignments_for_problems(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_screen_alignments_for_problems"
    output = os.path.join(o_dir, "mafft-gblocks-clean-problems-screened")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft-gblocks-clean-problems"),
        "--input-format",
        "nexus",
        "--output",
        output,
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = [
        os.path.basename(i) for i in glob.glob(os.path.join(output, "*"))
    ]
    expected_files = [
        os.path.basename(i)
        for i in glob.glob(
            os.path.join(e_dir, "mafft-gblocks-clean-problems-screened", "*")
        )
    ]
    assert set(output_files) == set(expected_files)


def test_align_split_concat_nexus_to_loci(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_split_concat_nexus_to_loci"
    output = os.path.join(o_dir, "mafft-gblocks-clean-concat-split-nexus")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--nexus",
        os.path.join(
            e_dir,
            "mafft-gblocks-clean-concat-nexus",
            "mafft-gblocks-clean-concat-nexus.nexus",
        ),
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    # make sure we get back the alignments that we started with when
    # concatenating
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-gblocks-clean", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected
