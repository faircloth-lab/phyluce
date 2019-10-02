#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 0 March 2012 09:03 PST (-0800)
"""


import logging
import os
import sys
import tempfile
import subprocess

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

from phyluce.pth import get_user_path
from phyluce.generic_align import GenericAlign

#import pdb


class Align(GenericAlign):
    """ MAFFT alignment class.  Subclass of GenericAlign which
    contains a majority of the alignment-related helper functions
    (trimming, etc.) """

    def __init__(self, input):
        """initialize, calling superclass __init__ also"""
        super(Align, self).__init__(input)
        self.log = logging.getLogger(__name__)
        console = logging.StreamHandler(sys.stdout)
        self.log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console.setFormatter(formatter)
        self.log.addHandler(console)


    def run_alignment(self, clean=True):
        # create results file
        fd, aln = tempfile.mkstemp(suffix='.mafft')
        os.close(fd)
        aln_stdout = open(aln, 'w')
        # run MAFFT on the temp file
        # Handle alternate mafft arguments
        mafft_args_env = os.getenv("PHYLUCE_MAFFT_ARGS")
        if mafft_args_env:
            mafft_args = mafft_args_env.strip().split()
            cmd = [get_user_path("binaries", "mafft")]
            cmd.extend(mafft_args)
            cmd.append(self.input)
        else:
            cmd = [get_user_path("binaries", "mafft"), "--adjustdirection", "--maxiterate", "1000", self.input]
        self.log.info("MAFFT Command: {}".format(cmd))
        # just pass all ENV params
        proc = subprocess.Popen(cmd,
                stderr=subprocess.PIPE,
                stdout=aln_stdout
            )
        stderr = proc.communicate()
        aln_stdout.close()
        self.alignment = AlignIO.read(open(aln, 'rU'), "fasta", \
                alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
        if clean:
            self._clean(aln)


if __name__ == '__main__':
    pass
