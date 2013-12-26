
import os
import sys
import glob
import argparse
import ConfigParser
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment as Alignment
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
            choices=['fasta', 'nexus', 'phylip', 'phylip-relaxed', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format"""
        )
    parser.add_argument(
            "--output-format",
            dest="output_format",
            choices=['fasta', 'nexus', 'phylip', 'phylip-relaxed', 'clustal', 'emboss', 'stockholm'],
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
            "--shorten-names",
            dest='shorten_name',
            action="store_true",
            default=False,
            help="""Convert names to a 6 or 7 character representation""",
        )
    parser.add_argument(
            "--name-conf",
            action=FullPaths,
            type=str,
            help="""A config-formatted file containing full-name:shortname mappings""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    extensions = get_file_extensions(input_format)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    # ensure we collapse duplicate filenames
    return list(set(files))


def test_if_name_in_keys(name, keys):
    if name in keys:
        for i in xrange(100):
            name = "{0}{1}".format(name, i)
            if name not in keys:
                break
            else:
                continue
    return name


def shorten_name(args, aln):
    aln = AlignIO.read(aln, args.input_format)
    names = {}
    for seq in aln:
        if "-" in seq.id:
            split_name = seq.id.split("-")
        elif "_" in seq.id:
            split_name = seq.id.split("_")
        elif " " in seq.id:
            split_name = seq.id.split(" ")
        else:
            split_name = None
        if split_name is not None:
            f3, l3 = split_name[0][0:3].title(), split_name[1][0:3].title()
            new_name = "{0}{1}".format(f3, l3)
        else:
            new_name = seq.id[:6]
        new_name = test_if_name_in_keys(new_name, names.keys())
        names[seq.id] = new_name
        seq.id, seq.name = new_name, new_name
    return names


def rename_alignment_taxa(aln, name_map):
    new_align = Alignment([], alphabet=Gapped(IUPAC.unambiguous_dna, "-"))
    for seq in aln:
        seq.id, seq.name = name_map[seq.id], name_map[seq.id]
        new_align.append(seq)
    return new_align


def convert_files_worker(params):
    f, args, name_map = params
    aln = AlignIO.read(f, args.input_format, alphabet=Gapped(IUPAC.ambiguous_dna))
    if args.shorten_name:
        aln = rename_alignment_taxa(aln, name_map)
    new_name = os.path.splitext(os.path.split(f)[1])[0] + '.{0}'.format(args.output_format)
    outf = open(os.path.join(args.outdir, new_name), 'w')
    AlignIO.write(aln, outf, args.output_format)
    outf.close()
    sys.stdout.write('.')
    sys.stdout.flush()


def main():
    args = get_args()
    files = get_files(args.indir, args.input_format)
    if args.shorten_name and not args.name_conf:
        name_map = shorten_name(args, files[0])
    elif args.shorten_name and args.name_conf:
        conf = ConfigParser.ConfigParser()
        conf.readfp(open(args.name_conf))
        name_map = dict(conf.items('taxa'))
    else:
        name_map = None
    params = [[f, args, name_map] for f in files]
    sys.stdout.write('Converting')
    sys.stdout.flush()
    if args.cores > 1:
        pool = Pool(args.cores)
        pool.map(convert_files_worker, params)
    else:
        map(convert_files_worker, params)
    if args.shorten_name:
        print "\n\nTaxa renamed (from) => (to):"
        for k, v in name_map.iteritems():
            print "\t{0} => {1}".format(k, v)

if __name__ == '__main__':
    main()
