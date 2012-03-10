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
from Bio.Alphabet import IUPAC, Gapped

from phyluce.helpers import which
from phyluce.generic_align import GenericAlign


class Align(GenericAlign):
    """ MUSCLE alignment class.  Subclass of GenericAlign which
    contains a majority of the alignment-related helper functions
    (trimming, etc.) """

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)

    def run_alignment(self, clean=True, consensus=True):
        """ muscle """
        # dialign requires ENV variable be set for dialign_dir
        muscle = which("muscle")
        # create results file
        fd, aln = tempfile.mkstemp(suffix='.muscle')
        os.close(fd)
        # run MUSCLE on the temp file
        cmd = [muscle, "-in", self.input, "-out", aln]
        proc = subprocess.Popen(cmd,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE
            )
        stdout, stderr = proc.communicate()
        self.alignment = AlignIO.read(open(aln, 'rU'), \
                "fasta", alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
        # build a dumb consensus
        if consensus:
            self.alignment_summary, self.alignment_consensus = \
                self._alignment_summary(self.alignment)
        # cleanup temp files
        if clean:
            self._clean(aln)


if __name__ == '__main__':
    pass
