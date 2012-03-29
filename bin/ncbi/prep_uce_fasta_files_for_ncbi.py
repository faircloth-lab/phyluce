#!/usr/bin/env python
# encoding: utf-8

"""
File: prep_uce_fasta_files_for_ncbi.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 February 2012 18:00 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Takes a fasta of contigs and the UCE match database
and outputs a merged fasta file with NCBI header info containing 
those contigs matching UCE loci in the species checked - ready for 
import to Sequin for => XML.

"""

import os
import glob
import sqlite3
import argparse
from seqtools.sequence.fasta import FastaReader, FastaWriter

import pdb

def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description = """use fasta + sqlite to format files for NCBI""")
    parser.add_argument('fastas', help="The directory containing fasta files")
    parser.add_argument('db', help="The db containing UCE matches")
    parser.add_argument('outfile', help="The outfile fasta to hold results")
    parser.add_argument('--fish', help="If working with fish data")
    return parser.parse_args()

def main():
    args = get_args()
    pth = os.path.join(args.fastas, "*.fasta")
    outf = FastaWriter(args.outfile)
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    counter = 0
    for infile in glob.glob(pth):
        sp = os.path.basename(infile).split('.')[0].replace('-','_')
        species = sp.replace('_',' ').capitalize()
        print "Working on {}".format(species)
        partial = species.split(' ')[0].lower()[:3]
        for read in FastaReader(infile):
            # check for header match, if match get locus name for header
            nn = read.identifier.split("_")[:2]
            nn = "{}_{}".format(nn[0].strip('>').lower(), nn[1].lower())
            query = "SELECT uce FROM match_map WHERE {0} = '{1}(+)' OR {0} = '{1}(-)'".format(sp, nn)
            cur.execute(query)
            result = cur.fetchall()
            #pdb.set_trace()
            if result:
                assert len(result) == 1, "More than 1 result"
                #pdb.set_trace()
                if args.fish:
                    uce = result[0][0].split('_')[0]
                else:
                    uce = result[0][0]
                read.identifier = """{3}{2} [organism={0}] [molecule=DNA] [moltype=genomic] [location=genomic] [note=ultra conserved element locus {1}] {0} ultra-conserved element locus {1}.""".format(species, uce, partial, counter)
                # write all to a common fasta
                outf.write(read)
                # if not match, pass
                counter += 1
            else:
                pass
    outf.close()

if __name__ == '__main__':
    main()