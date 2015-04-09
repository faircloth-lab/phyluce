#!/usr/bin/env python
# encoding: utf-8

"""
get_iids_from_cvg_stat.py

Created by Brant Faircloth on 28 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import argparse

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get the contig iid from the cvgStat file output by AMOS')
    parser.add_argument('cvg', help='The cvgStat output file', type=argparse.FileType('rU'))
    parser.add_argument('nodes', help='The taxa-specific nodes matching to probes', type=argparse.FileType('rU'))
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()

def main():
    args = get_args()
    # get the node names for given taxa
    nodenames = set([node.split('_')[1].split('(')[0] for node in args.nodes if node.startswith('node')])
    print "{0} distinct loci".format(len(nodenames))
    d = {}
    for line in args.cvg:
        if line.startswith('>'):
            ls = line.strip()
            lsp = ls.split(' ')
            contig = lsp[0].lstrip('>').split('-')[0]
            iid = lsp[2].split(':')[1]
            if contig in nodenames:
                d[contig] = iid
    for k in sorted(d.keys()):
        args.output.write("{0}\n".format(d[k]))
    args.output.close()
    
if __name__ == '__main__':
    main()
    