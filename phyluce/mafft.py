#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 0 March 2012 09:03 PST (-0800)
"""


import os
import tempfile
import subprocess

from Bio import AlignIO

from phyluce.pth import get_user_path
from phyluce.generic_align import GenericAlign

# import pdb


class Align(GenericAlign):
    """MAFFT alignment class.  Subclass of GenericAlign which
    contains a majority of the alignment-related helper functions
    (trimming, etc.)"""

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)

    def run_alignment(self, clean=True):
        # create results file
        fd, aln = tempfile.mkstemp(suffix=".mafft")
        os.close(fd)
        aln_stdout = open(aln, "w")
        # run MAFFT on the temp file
        cmd = [
            get_user_path("binaries", "mafft"),
            "--adjustdirection",
            "--maxiterate",
            "1000",
            self.input,
        ]
        # just pass all ENV params
        proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=aln_stdout)
        proc.communicate()
        aln_stdout.close()
        self.alignment = AlignIO.read(open(aln, "rU"), "fasta")
        # we now need to set the molecule type for biopython
        # due to removal of seq.alphabet
        for seq in self.alignment:
            seq.annotations = {"molecule_type": "DNA"}
        if clean:
            self._clean(aln)


if __name__ == "__main__":
    pass
