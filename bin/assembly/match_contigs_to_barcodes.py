#!/usr/bin/env python
# encoding: utf-8

"""
match_contigs_to_probes.py

Created by Brant Faircloth on 02 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import glob
import argparse
from xml.etree import ElementTree

import requests


from phyluce import lastz
from phyluce.helpers import is_dir, is_file, FullPaths, CreateDir
from phyluce.log import setup_logging

from Bio import SeqIO

import pdb


def get_args():
    parser = argparse.ArgumentParser(description=""
    )
    parser.add_argument(
        '--contigs',
        required=True,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the assembled contigs in which you are searching for UCE loci."
    )
    parser.add_argument(
        '--barcodes',
        required=True,
        type=is_file,
        action=FullPaths,
        help="A FASTA-formatted file of barcode(s)."
    )
    parser.add_argument(
        '--output',
        required=True,
        action=CreateDir,
        help="The directory in which to store the resulting SQL database and LASTZ files."
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    log, my_name = setup_logging(args)
    fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
    for contig in fasta_files:
        # get a name for the output file
        critter = os.path.basename(contig).split('.')[0].replace('-', "_")
        text = " Processing {} ".format(critter)
        log.info(text.center(65, "-"))
        log.info("Parsing FASTA {}".format(os.path.basename(contig)))
        # go ahead and parse the fasta into a dictionary
        record_dict = SeqIO.to_dict(SeqIO.parse(
            open(contig, "rU"),
            "fasta"
        ))
        log.info("Running LASTZ".format(critter))
        output = os.path.join(
            args.output,
            os.path.splitext(os.path.basename(contig))[0] + '.lastz'
        )
        # align
        alignment = lastz.SimpleAlign(
            args.barcodes,
            contig,
            output
        )
        lzstdout, lztstderr = alignment.run()
        if lztstderr:
            raise EnvironmentError("lastz: {}".format(lztstderr))
        else:
            fasta_slices = []
            log.info("Running match against BOLD")
            for lz in lastz.Reader(output):
                try:
                    matching_contig = lz.name2.split(' ')[0]
                    matching_contig_strand = lz.strand2
                    matching_contig_sequence = record_dict[matching_contig]
                    if matching_contig_strand == "+":
                        slice = matching_contig_sequence[lz.zstart2:lz.end2]
                        matches_barcode_sequence = str(slice.seq)
                    elif matching_contig_strand == "-":
                        rev_matching_contig_sequence = matching_contig_sequence.reverse_complement()
                        slice = rev_matching_contig_sequence[lz.zstart2:lz.end2]
                        matches_barcode_sequence = str(slice.seq)
                    # keep fasta slice
                    slice.id = matching_contig
                    fasta_slices.append(slice)
                    payload = {
                        "sequence":matches_barcode_sequence,
                        "db":"COX1_SPECIES"
                    }
                    r = requests.get("http://boldsystems.org/index.php/Ids_xml", params=payload)
                    root = ElementTree.fromstring(r.content)
                    id = [elem.text for elem in root.findall('match/ID')]
                    matches = [elem.text for elem in root.findall('match/taxonomicidentification')]
                    similarity = [elem.text for elem in root.findall('match/similarity')]
                    log.info("\tBest BOLD systems match for locus {0}: {1} [{2}] ({3})".format(
                        matching_contig,
                        matches[0],
                        id[0],
                        similarity[0]
                    ))
                except IndexError:
                    log.warn("Did not find a match for locus {}".format(matching_contig))
            output = os.path.join(
                args.output,
                os.path.splitext(os.path.basename(contig))[0] + '.slices.fasta'
            )
            with open(output, "w") as outf:
                SeqIO.write(fasta_slices, outf, "fasta")
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))




if __name__ == '__main__':
    main()
