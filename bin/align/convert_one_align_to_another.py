
import os
import sys
import glob
import argparse
from Bio import AlignIO
from multiprocessing import Pool
from Bio.Alphabet import IUPAC, Gapped
from phyluce.helpers import get_file_extensions, is_dir, FullPaths

import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Align records in a file of UCE fastas""")
    parser.add_argument('indir',
            type=is_dir,
            action=FullPaths,
            help='The file containing fasta reads associated with UCE loci')
    parser.add_argument('outdir',
            type=is_dir,
            action=FullPaths,
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
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of compute cores to use"""
        )
    parser.add_argument(
            "--wrap",
            action="store_true",
            default=False,
            help="""Wrap fasta alignments""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    extensions = get_file_extensions(input_format)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    # ensure we collapse duplicate filenames
    return list(set(files))


def convert_files_worker(params):
    f, args = params
    align = AlignIO.read(f, args.input_format, alphabet=Gapped(IUPAC.ambiguous_dna))
    new_name = os.path.splitext(os.path.split(f)[1])[0] + '.{0}'.format(args.output_format)
    outf = open(os.path.join(args.outdir, new_name), 'w')
    AlignIO.write(align, outf, args.output_format)
    outf.close()
    sys.stdout.write('.')
    sys.stdout.flush()


def main():
    args = get_args()
    files = get_files(args.indir, args.input_format)
    params = [[f, args] for f in files]
    sys.stdout.write('Converting')
    sys.stdout.flush()
    if args.cores > 1:
        pool = Pool(args.cores)
        pool.map(convert_files_worker, params)
    else:
        map(convert_files_worker, params)

if __name__ == '__main__':
    main()
