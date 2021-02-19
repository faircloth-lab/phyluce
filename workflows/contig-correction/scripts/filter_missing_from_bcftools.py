#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq

# import pdb

for input, output in zip(snakemake.input, snakemake.output):
    with open(output, "w") as outf:
        for seq in SeqIO.parse(open(input), "fasta"):
            # replace the missing data dot with IUPAC N
            temp_seq_str = str(seq.seq).replace(".", "N")
            # strip the leading Ns
            temp_seq_str = temp_seq_str.lstrip("N")
            # strip the trailing Ns
            temp_seq_str = temp_seq_str.rstrip("N")
            # finally, convert to sequence and remove those < 50 bp
            # after filtering
            seq.seq = Seq(temp_seq_str)
            if len(seq) > 50:
                outf.write(format(seq, "fasta"))
