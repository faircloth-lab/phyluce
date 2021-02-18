#!/usr/bin/env python
# encoding: utf-8
"""
File: multi_lastz.py
Author: Brant Faircloth

Created by Brant Faircloth on 31 August 2012 14:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:

"""

import sys
import os
import tempfile
import subprocess
import multiprocessing

from phyluce.pth import get_user_path

from bx.seq import twobit

import pdb


def run_lastz(work):
    """ """
    unit, probes, coverage, identity = work
    temp_fd, temp_out = tempfile.mkstemp(suffix=".lastz")
    os.close(temp_fd)
    cmd = lastz_params(unit, probes, coverage, identity, temp_out)
    lzstdout, lztstderr = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate(None)
    if lztstderr:
        raise IOError("Lastz returned:\n{0}".format(lztstderr))
    # don't keep empty files
    if os.stat(temp_out)[6] == 0:
        os.remove(temp_out)
        return None
    else:
        sys.stdout.write(".")
        sys.stdout.flush()
        return temp_out


def lastz_params(target, query, coverage, identity, outfile):
    output_format = "general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity,coverage"
    cmd = [
        get_user_path("binaries", "lastz"),
        "{0}[multiple]".format(target),
        "{0}[nameparse=full]".format(query),
        "--strand=both",
        "--seed=12of19",
        "--transition",
        "--nogfextend",
        "--nochain",
        "--gap=400,30",
        "--xdrop=910",
        "--ydrop=8370",
        "--hspthresh=3000",
        "--gappedthresh=3000",
        "--noentropy",
        "--coverage={0}".format(coverage),
        "--identity={0}".format(identity),
        "--output={0}".format(outfile),
        "--format={0}".format(output_format),
    ]
    return cmd


def chunk_scaffolds(log, target, size):
    chromos = []
    # split target file into `options.size` (~10 Mbp) chunks
    temp_fd, temp_out = tempfile.mkstemp(suffix=".fasta")
    os.close(temp_fd)
    temp_out_handle = open(temp_out, "w")
    with open(target, "rb") as f:
        tb = twobit.TwoBitFile(f)
        sequence_length = 0
        tb_key_len = len(tb.keys()) - 1
        log.info("Running against {}".format(os.path.basename(target)))
        log.info(
            "Running with the --huge option.  Chunking files into {0} bp...".format(
                size
            )
        )
        for sequence_count, seq in enumerate(tb.keys()):
            sequence = tb[seq][0:]
            sequence_length += len(sequence)
            # write it to the outfile
            temp_out_handle.write(
                ">{0}\n{1}\n".format(seq.decode("utf-8"), sequence)
            )
            if sequence_length > size:
                temp_out_handle.close()
                # put tempfile name on stack
                chromos.append(temp_out)
                # open a new temp file
                temp_fd, temp_out = tempfile.mkstemp(suffix=".fasta")
                os.close(temp_fd)
                temp_out_handle = open(temp_out, "w")
                # reset sequence length
                sequence_length = 0
            # if we hit the end of the twobit file
            elif sequence_count >= tb_key_len:
                temp_out_handle.close()
                # put tempfile name on stack
                chromos.append(temp_out)
            else:
                pass
        temp_out_handle.close()
    return chromos


def multi_lastz_runner(
    log,
    output,
    cores,
    target,
    query,
    huge,
    coverage=83,
    identity=92.5,
    size=10000000,
):
    if not huge:
        # get individual records from the 2bit file
        with open(target, "rb") as f:
            tb = twobit.TwoBitFile(f)
            tb_keys = [i.decode("utf-8") for i in tb.keys()]
        chromos = [os.path.join(target, c) for c in tb_keys]
    else:
        chromos = chunk_scaffolds(log, target, size)
    work = [[chromo, query, coverage, identity] for chromo in chromos]
    # pdb.set_trace()
    log.info("Running the targets against %s queries..." % len(chromos))
    if cores == 1:
        results = list(map(run_lastz, work))
    else:
        pool = multiprocessing.Pool(cores)
        results = pool.map(run_lastz, work)
    print("")
    log.info("Writing the results file...")
    outp = open(output, "wb")
    for t_file in results:
        if t_file is not None:
            sys.stdout.write(".")
            sys.stdout.flush()
            # read the file
            outp.write(open(t_file, "rb").read())
            # cleanup the lastz output files
            os.remove(t_file)
    outp.close()
    print("")
    if huge:
        log.info("Cleaning up the chunked files...")
        for t_file in chromos:
            # cleanup the chunked files
            os.remove(t_file)
