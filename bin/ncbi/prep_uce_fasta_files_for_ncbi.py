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
import ConfigParser
from seqtools.sequence.fasta import FastaReader, FastaWriter
from phyluce.helpers import is_dir, is_file, FullPaths

import pdb


def get_args():
    """get arguments (config file location)"""
    parser = argparse.ArgumentParser(description="""use fasta + sqlite to format files for NCBI""")
    parser.add_argument(
        'fastas',
        action=FullPaths,
        type=is_dir,
        help="The directory containing fasta files"
        )
    parser.add_argument(
        'db',
        action=FullPaths,
        type=is_file,
        help="The db containing UCE matches"
        )
    parser.add_argument(
        'conf',
        action=FullPaths,
        type=is_file,
        help="A configuration file containing general and sample-specific metadata"
        )
    parser.add_argument(
        'outfile',
        action=FullPaths,
        help="The outfile fasta to hold results"
        )
    parser.add_argument(
        '--fish',
        help="If working with fish data"
        )
    parser.add_argument(
        '--start-value',
        dest="start_value",
        type=int,
        default=0,
        help="The starting index value to append to sequences"
        )
    return parser.parse_args()


def get_excludes(conf, sec):
    if conf.has_section(sec):
        excludes = [i[0] for i in conf.items(sec)]
    else:
        return [None]
    return excludes


def get_metadata(conf):
    return dict(conf.items('metadata'))


def get_vouchers(conf):
    if conf.has_section('vouchers'):
        return dict(conf.items('vouchers'))
    else:
        return None


def get_remaps(conf):
    if conf.has_section('remap'):
        return {i[0].replace(' ', '_'): i[1] for i in conf.items('remap')}
    else:
        return None


def get_species_name(infile, remap):
    sp = os.path.basename(infile).split('.')[0].replace('-', '_')
    if remap is not None and sp in remap.keys():
        oldname = sp
        sp = remap[sp]
    else:
        oldname = sp
    species = sp.replace('_', ' ').capitalize()
    partial = species.split(' ')[0].lower()[:3]
    return sp, species, partial, oldname


def get_node_name(read):
    # check for header match, if match get locus name for header
    nn_split = read.identifier.split("_")[:2]
    nn = "{}_{}".format(nn_split[0].strip('>').lower(), nn_split[1].lower())
    return nn


def get_new_identifier(species, uce, partial, counter, metadata, vouchers):
    title = "{0} ultra-conserved element locus {1}".format(species, uce)
    note = metadata['note'].format(uce)
    metadata = {
                "counter":counter,
                "partial":partial,
                "title":title,
                "organism":"{0}".format(species),
                "moltype":"{0}".format(metadata['moltype']),
                "location":"{0}".format(metadata['location']),
                "note":note,
                "specimen_voucher":"{0}".format(vouchers[species.lower()].replace(" ", ":"))
            }
    #pdb.set_trace()
    new_identifier = "{counter}{partial} [organism={organism}] [moltype={moltype}] [location={location}] [note={note}] [specimen-voucher={specimen_voucher}] {title}".format(**metadata)
    return new_identifier


def main():
    args = get_args()
    conf = ConfigParser.ConfigParser(allow_no_value=True)
    conf.read(args.conf)
    # get metadata from conf file
    taxon_excludes = get_excludes(conf, "exclude taxa")
    locus_excludes = get_excludes(conf, "exclude loci")
    metadata = get_metadata(conf)
    vouchers = get_vouchers(conf)
    #pdb.set_trace()
    remap = get_remaps(conf)
    # get fasta and db locations
    pth = os.path.join(args.fastas, "*.fasta")
    outf = FastaWriter(args.outfile)
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    counter = args.start_value
    # iterate over fasta files
    for infile in glob.glob(pth):
        sp, species, partial, oldname = get_species_name(infile, remap)
        if species.lower() not in taxon_excludes:
            print "Working on {}".format(species)
            for read in FastaReader(infile):
                nodename = get_node_name(read)
                query = "SELECT uce FROM match_map WHERE {0} = '{1}(+)' OR {0} = '{1}(-)'".format(oldname, nodename)
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    # ensure we get only 1 result
                    assert len(result) == 1, "More than 1 result"
                    # if getting fish data TODO: deprecate
                    if args.fish:
                        uce = result[0][0].split('_')[0]
                    else:
                        uce = result[0][0]
                    if uce not in locus_excludes:
                        read.identifier = get_new_identifier(species, uce, partial, counter, metadata, vouchers)
                        #read.identifier = """{3}{2} [organism={0}] [molecule=DNA] [moltype=genomic] [location=genomic] [note=ultra conserved element locus {1}] {0} ultra-conserved element locus {1}.""".format(species, uce, partial, counter)
                        # write all to a common fasta
                        outf.write(read)
                        # if not match, pass
                        counter += 1
                else:
                    pass
        else:
            print "Skipping {0}".format(species)
    outf.close()

if __name__ == '__main__':
    main()