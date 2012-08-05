"""
File: get_AMOS_coverage_for_contigs.py
Author: Brant Faircloth

Created by Brant Faircloth on 27 February 2012 13:02 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Get overall contig coverage and UCE-specific coverage for 
velvet contigs output during the assembly process.  Assumes that directories
are generally of the following structure:

genus-species-1/
    auto_data_69/
        velvet_asm.afg
genus-species-2/
    auto_data_71/
        velvet_asm.afg
genus-species-3/
    auto_data_69/
        velvet_asm.afg

The directory structure can also be (output by assemblo.py, in this package):

genus-species-1/
    assembly/
        auto_data_69/
            velvet_asm.afg
genus-species-2/
    assembly/
        auto_data_71/
            velvet_asm.afg
genus-species-3/
    assembly/
        auto_data_69/
            velvet_asm.afg

If run in TLD, will walk file-tree, working on any file in a subdir that
ends with *.afg.  Must pass UCE match database to compute coverage across
UCE loci.  Probably you need to be strict about file structure above, as
dealing with filenames, etc. is rather fragile.

Execution:

    python get_AMOS_coverage_for_contigs.py ./ results.csv 
        --db probe.matches.sqlite
        
"""

import os
import sys
import shutil
import sqlite3
import argparse
import subprocess

import pdb

class Unbuffered:
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def get_args():
    parser = argparse.ArgumentParser(description="""Determine the coverage
         of contigs using AMOS.""")
    parser.add_argument('contigs',
            help='The folder containing contigs')
    parser.add_argument('outfile', type=argparse.FileType('w'),
            help='A file for the output.')
    parser.add_argument('--output',
            help='The folder containing contigs')
    parser.add_argument('--db',
            help='Path to the database containing UCE matches')
    return parser.parse_args()


def create_bnk(amos, infile, bank):
    print "\tCreating AMOS bnk file..."
    cmd = [os.path.join(amos, 'Bank/bank-transact'), '-m', infile, '-b', bank, '-c']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    output = process.communicate()
    return output


def get_overall_cvg(amos, infile, bank, outfile):
    print "\tGetting read depth..."
    cmd = [os.path.join(amos, 'Validation/analyze-read-depth'), bank, '-d']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    results = stdout.split()
    outfile.write("{},{},{},{},".format(infile, results[1], results[2], results[3]))


def get_cvg_stat(amos, infile, bank):
    print "\tGetting all loci cvgStat..."
    cmd = [os.path.join(amos, 'Utils/cvgStat'), '-b', bank]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout


def get_iids_from_cvg_stat(cvg, nodenames, fullpath):
    """parse cvg st output into iid file"""
    d = {}
    outpath = os.path.dirname(fullpath)
    taxon = os.path.dirname(os.path.dirname(fullpath)).replace('-', '_')
    if '/assembly' in taxon:
        taxon = os.path.dirname(taxon)
    outfile_name = os.path.join(outpath, "{}.iid".format(taxon))
    outfile = open(outfile_name, 'w')
    for line in cvg.split('\n'):
        if line.startswith('>'):
            ls = line.strip()
            lsp = ls.split(' ')
            contig = lsp[0].lstrip('>').split('-')[0]
            iid = lsp[2].split(':')[1]
            if contig in nodenames:
                d[contig] = iid
    for k in sorted(d.keys()):
        outfile.write("{0}\n".format(d[k]))
    outfile.close()
    return outfile_name


def get_loci_list(db, fullpath):
    """Get list of UCE loci for taxon from sqlite"""
    #pdb.set_trace()
    outpath = os.path.dirname(fullpath)
    taxon = os.path.dirname(os.path.dirname(fullpath)).replace('-', '_').strip('./')
    if '/assembly' in taxon:
        taxon = os.path.dirname(taxon)
    outfile = open(os.path.join(outpath, "{}.loci".format(taxon)), 'w')
    conn = sqlite3.connect(db)
    cur = conn.cursor()
    nodenames = []
    query = "SELECT {0} FROM match_map WHERE {0} IS NOT NULL".format(taxon)
    for result in cur.execute(query):
        outfile.write("{}\n".format(result[0]))
        if result[0].startswith('node'):
            nodenames.append(result[0].split('_')[1].split('(')[0])
    outfile.close()
    nodenames = set(nodenames)
    print "\t{0} distinct UCE loci".format(len(nodenames))
    return nodenames


def get_uce_cvg(amos, infile, bank, iid, outfile):
    print "\tGetting UCE cvgStat..."
    cmd = [os.path.join(amos, 'Validation/analyze-read-depth'), bank, '-d', '-I', iid]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdsplit = stdout.strip().split()
    #pdb.set_trace()
    outfile.write("{0},{1},{2}\n".format(stdsplit[1], stdsplit[2], stdsplit[3]))


def main():
    args = get_args()
    amos = "/Users/bcf/Source/amos-3.0.0/src/"
    args.outfile.write("filename,contigs,bp-in-contigs,contig-coverage,uce-contigs,bp-in-uce-contigs,uce-coverage\n")
    for root, dirs, infiles in os.walk(args.contigs):
        for infile in infiles:
            fullpath = os.path.join(root, infile)
            if infile == "velvet_asm.afg":
                print "Working on {}".format(fullpath)
                if not args.output:
                    outdir = root
                else:
                    outdir = args.output
                bank = os.path.join(outdir, "velvet_asm.bnk")
                # make sure we don't overwrite bnk files when we don't need to
                if os.path.exists(bank):
                    answer = raw_input("\tBNK file exists, overwrite [Y/n]? ")
                    if answer == "Y":
                        shutil.rmtree(bank)
                        create_bnk(amos, fullpath, bank)
                    else:
                        answer = raw_input("\tWould you like to proceed" \
                            + " with the existing BNK file [Y/n]? ")
                        if answer != "Y":
                            sys.exit()
                else:
                    create_bnk(amos, fullpath, bank)
                # get the coverage
                get_overall_cvg(amos, fullpath, bank, args.outfile)
                # for UCE loci, get the coverage of those loci, only:
                if args.db:
                    cvg = get_cvg_stat(amos, fullpath, bank)
                    nodenames = get_loci_list(args.db, fullpath)
                    iid = get_iids_from_cvg_stat(cvg, nodenames, fullpath)
                    uce_cvg = get_uce_cvg(amos, infile, bank, iid, args.outfile)                  
                else:
                    pass
    args.outfile.close()

if __name__ == '__main__':
    main()
