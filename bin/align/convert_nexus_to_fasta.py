
import os
import glob
import argparse
from Bio import AlignIO

#import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Align records in a file of UCE fastas""")
    parser.add_argument('indir',
            help='The file containing fasta reads associated with UCE loci')
    parser.add_argument('outdir',
            help='A directory for the output.')
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))


def main():
    args = get_args()
    files = get_files(args.indir)
    # iterate through all the files to determine the taxa present
    #taxa = get_and_sort_all_taxa(files)
    #pdb.set_trace()
    for count, f in enumerate(files):
        align = AlignIO.parse(f, "nexus")
        #pdb.set_trace()
        new_name = os.path.splitext(os.path.split(f)[1])[0] + '.fasta'
        #pdb.set_trace()
        outf = open(os.path.join(args.outdir, new_name), 'w')
        #AlignIO.write(align, open(outf, 'w'), 'fasta', wrap = None)
        td = {}
        for a in align:
            for seq in a:
                td[seq.name] = seq.seq
        for o in sorted(td.keys()):
                outf.write(">{}\n{}\n".format(o, td[o]))
        print count
        outf.close()

if __name__ == '__main__':
    main()
