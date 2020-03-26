#!/usr/bin/env python
# encoding: utf-8

"""
Lastz.py

Created by Brant Faircloth on 01 May 2010 19:01 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  This is a helper class to simplify calls to lastz.

USAGE:

    import Lastz

    lastz = Lastz.Align(target, query, matchcount, identity, output)
    lastz_stderr, lastz_stdout = lastz.run()
    results_file = lastz.output

"""

import os
import tempfile
import subprocess
from collections import namedtuple
from collections import defaultdict

from phyluce.pth import get_user_path

import pdb


class SimpleAlign:
    """docstring for lastz"""

    def __init__(self, target, query, out=False):
        # if not an output file, create a temp file to hold output
        if not out:
            fd, self.output = tempfile.mkstemp(suffix=".lastz")
            os.close(fd)
        else:
            self.output = out
        self.cli = "{3} {0}[multiple,nameparse=full] {1}[nameparse=full]\
                --output={2} \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity".format(
            target, query, self.output, get_user_path("binaries", "lastz")
        )

    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).communicate(None)
        return lastz_stdout, lastz_stderr


class Align:
    """docstring for lastz"""

    def __init__(
        self, target, query, coverage, identity, out=False, min_match=None
    ):
        # if not an output file, create a temp file to hold output
        if not out:
            fd, self.output = tempfile.mkstemp(suffix=".lastz")
            os.close(fd)
        else:
            self.output = out
        if identity and not min_match:
            self.cli = "{5} {0}[multiple,nameparse=full] {1}[nameparse=full]\
                --strand=both \
                --seed=12of19 \
                --transition \
                --nogfextend \
                --nochain \
                --gap=400,30 \
                --xdrop=910 \
                --ydrop=8370 \
                --hspthresh=3000 \
                --gappedthresh=3000 \
                --noentropy \
                --coverage={2} \
                --identity={3} \
                --output={4} \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity".format(
                target,
                query,
                coverage,
                identity,
                self.output,
                get_user_path("binaries", "lastz"),
            )
        elif min_match:
            self.cli = "{5} {0}[multiple,nameparse=full] {1}[nameparse=full]\
                --strand=both \
                --seed=12of19 \
                --transition \
                --nogfextend \
                --nochain \
                --gap=400,30 \
                --xdrop=910 \
                --ydrop=8370 \
                --hspthresh=3000 \
                --gappedthresh=3000 \
                --noentropy \
                --matchcount={2} \
                --identity={3} \
                --output={4} \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity".format(
                target,
                query,
                min_match,
                identity,
                self.output,
                get_user_path("binaries", "lastz"),
            )

    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).communicate(None)
        return lastz_stdout, lastz_stderr


class Reader:
    """read a lastz file and return an iterator over that file"""

    def __init__(self, lastz_file, long_format=False):
        self.file = open(lastz_file, "rU")
        self.long_format = long_format

    def __del__(self):
        """close files"""
        self.file.close()

    def __iter__(self):
        """iterator"""
        while True:
            yield next(self)

    def __next__(self):
        """read next fastq sequence and return as named tuple"""
        lastz_result = self.file.readline()
        if not lastz_result:
            raise StopIteration
        if not self.long_format:
            Lastz = namedtuple(
                "Lastz",
                "score,name1,strand1,zstart1,end1,length1,name2,"
                + "strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,"
                + "continuity,percent_continuity",
            )
        else:
            Lastz = namedtuple(
                "Lastz",
                "score,name1,strand1,zstart1,end1,length1,name2,"
                + "strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,"
                + "continuity,percent_continuity,coverage,percent_coverage",
            )
        lastz_result_split = lastz_result.strip("\n").split("\t")
        for k, v in enumerate(lastz_result_split):
            if k in [3, 4, 5, 8, 9, 10]:
                lastz_result_split[k] = int(v)
            elif "%" in v:
                lastz_result_split[k] = float(v.strip("%"))
        lastz_result_split[1] = lastz_result_split[1].lstrip(">")
        lastz_result_split[6] = lastz_result_split[6].lstrip(">")
        return Lastz._make(lastz_result_split)


if __name__ == "__main__":
    pass
