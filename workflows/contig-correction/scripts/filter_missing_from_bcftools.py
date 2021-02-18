#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq

for seq in SeqIO.parse(sys.stdin, "fasta"):
    seq.seq = Seq(str(seq.seq).replace(".", ""))
    if len(seq) > 50:
        sys.stdout.write(seq.format("fasta"))
