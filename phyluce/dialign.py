#!/usr/bin/env python
# encoding: utf-8
"""
File: dialign.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 March 2012 11:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import tempfile
import subprocess

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

from phyluce.helpers import which
from phyluce.generic_align import GenericAlign

import pdb



class Align(GenericAlign):
    """ text """

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)

    def run_alignment(self, clean=True, consensus=True):
        # dialign requires ENV variable be set for dialign_dir
        daln = which("dialign2-2")
        daln = os.path.join(os.path.split(daln)[0], 'dialign2_dir')
        os.environ["DIALIGN2_DIR"] = "/Users/bcf/Bin/dialign2_dir/"
        # create results file
        fd, aln = tempfile.mkstemp(suffix='.dialign')
        os.close(fd)
        # dialign makes an extra file for fasta output
        fasta = "{}.{}".format(aln, 'fa')
        # run MUSCLE on the temp file
        cmd = ["dialign2-2", "-fa", "-fn", aln, "-n", self.input]
        # just pass all ENV params
        proc = subprocess.Popen(cmd,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                env=os.environ
            )
        stdout, stderr = proc.communicate()
        self.alignment = AlignIO.read(open(fasta, 'rU'), "fasta", \
                alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
        # build a dumb consensus
        if consensus:
            self.alignment_summary, self.alignment_consensus = \
                self._alignment_summary(self.alignment)
        # cleanup temp files
        if clean:
            self._clean([aln, fasta])


if __name__ == '__main__':
    pass
