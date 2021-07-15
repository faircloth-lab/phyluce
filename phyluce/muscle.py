#!/usr/bin/env python
# encoding: utf-8
"""
File: muscle.py
Author: Brant Faircloth

Created by Brant Faircloth on 10 March 2012 10:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:  Alignement class wrapping MUSCLE aligner
(http://www.drive5.com/muscle/)

"""

import os
import tempfile
import subprocess

from Bio import AlignIO

from phyluce.pth import get_user_path
from phyluce.generic_align import GenericAlign


class Align(GenericAlign):
    """MUSCLE alignment class.  Subclass of GenericAlign which
    contains a majority of the alignment-related helper functions
    (trimming, etc.)"""

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)

    def run_alignment(self, clean=True):
        """muscle"""
        # create results file
        fd, aln = tempfile.mkstemp(suffix=".muscle")
        os.close(fd)
        # run MUSCLE on the temp file
        cmd = [
            get_user_path("binaries", "muscle"),
            "-in",
            self.input,
            "-out",
            aln,
        ]
        proc = subprocess.Popen(
            cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        self.alignment = AlignIO.read(open(aln, "rU"), "fasta")
        # we now need to set the molecule type for biopython
        # due to removal of seq.alphabet
        for seq in self.alignment:
            seq.annotations = {"molecule_type": "DNA"}
        # cleanup temp files
        if clean:
            self._clean(aln)


if __name__ == "__main__":
    pass
