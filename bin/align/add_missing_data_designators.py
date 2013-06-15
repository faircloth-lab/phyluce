#!/usr/bin/env python
# encoding: utf-8

"""
add_gaps_for_missing_taxa.py

Created by Brant Faircloth on 05 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import copy
import glob
import argparse
import ConfigParser
from collections import defaultdict
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument(
            'input',
            help="Alignment files to process"
        )
    parser.add_argument(
            'output',
            nargs='?',
            default=sys.stdout,
            help="Output dir for alignment files"
        )
    parser.add_argument(
            'config',
            default=None,
            help='The config file containing match information'
        )
    parser.add_argument(
            '--notstrict',
            default=None,
            help="Name of conf file containing notstrict species and data"
            )
    parser.add_argument(
            '--min-taxa',
            dest='min_taxa',
            help="The minimum number of taxa to keep (default = 3)",
            default=3,
            type=int
        )
    parser.add_argument(
            '--verbatim',
            action="store_true",
            default=False,
            help="""Do not parse species names at all - use them verbatim""",
        )
    return parser.parse_args()


def get_names_from_config(config, group):
    try:
        return [i[0].rstrip('*') for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None

def record_formatter(seq, name):
        """return a string formatted as a biopython sequence record"""
        return SeqRecord(Seq(seq, Gapped(IUPAC.ambiguous_dna, "-?")),
            id=name,
            name=name,
            description=name)

def add_gaps_to_align(organisms, missing, align, verbatim=False, min_taxa=3):
    local_organisms = copy.deepcopy(organisms)
    for a in align:
        if len(a) < min_taxa:
            new_align = None
            break
        elif len(a) >= min_taxa:
            new_align = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna, "-?"))
            overall_length = len(a[0])
            for seq in a:
                # strip any reversal characters from mafft
                seq.name = seq.name.lstrip('_R_')
                if not verbatim:
                    new_seq_name = '_'.join(seq.name.split('_')[1:])
                else:
                    new_seq_name = seq.name.lower()
                new_align.append(record_formatter(str(seq.seq), new_seq_name))
                local_organisms.remove(new_seq_name)
            for org in local_organisms:
                if not verbatim:
                    loc = '_'.join(seq.name.split('_')[:1])
                else:
                    loc = seq.name
                if missing:
                    try:
                        assert loc in missing[org], "Locus missing"
                    except:
                        assert loc in missing['{}*'.format(org)], "Locus missing"
                missing_string = '?' * overall_length
                new_align.append(record_formatter(missing_string, org))
    return new_align


def get_missing_loci_from_conf_file(config):
    missing = defaultdict(list)
    for sec in config.sections():
    #    pdb.set_trace()
        for item in config.items(sec):
            missing[sec].append(item[0])
    return missing


def main():
    args = get_args()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.config)
    if args.notstrict:
        notstrict = ConfigParser.RawConfigParser(allow_no_value=True)
        notstrict.read(args.notstrict)
        missing = get_missing_loci_from_conf_file(notstrict)
    else:
        missing = None
    organisms = get_names_from_config(config, 'Organisms')
    for count, nex in enumerate(glob.glob(os.path.join(args.input, '*.nex*'))):
        align = AlignIO.parse(nex, "nexus")
        new_align = add_gaps_to_align(organisms, missing, align, args.verbatim, args.min_taxa)
        if new_align is not None:
            outf = os.path.join(args.output, os.path.basename(nex))
            AlignIO.write(new_align, open(outf, 'w'), 'nexus')
            print count
        else:
            print "{0} Skipping {1} - < {2} taxon".format(count, nex,
                    args.min_taxa)

if __name__ == '__main__':
    main()
