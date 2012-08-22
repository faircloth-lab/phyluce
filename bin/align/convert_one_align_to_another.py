
import os
import glob
import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from phyluce.helpers import get_file_extensions

import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Align records in a file of UCE fastas""")
    parser.add_argument('indir',
            help='The file containing fasta reads associated with UCE loci')
    parser.add_argument('outdir',
            help='A directory for the output.')
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format"""
        )
    parser.add_argument(
            "--output-format",
            dest="output_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='nexus',
            help="""The input alignment format"""
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    extensions = get_file_extensions(input_format)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    return files


def main():
    args = get_args()
    files = get_files(args.indir, args.input_format)
    #pdb.set_trace()
    for count, f in enumerate(files):
        align = AlignIO.read(f, args.input_format, alphabet=Gapped(IUPAC.ambiguous_dna))
        new_name = os.path.splitext(os.path.split(f)[1])[0] + '.{0}'.format(args.output_format)
        outf = open(os.path.join(args.outdir, new_name), 'w')
        AlignIO.write(align, outf, args.output_format)
        outf.close()
        print count

if __name__ == '__main__':
    main()
