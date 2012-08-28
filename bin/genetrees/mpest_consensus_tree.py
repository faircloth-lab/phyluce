#!/usr/bin/env python
# encoding: utf-8
"""
consensus_tree.py

Created by Nick Crawford on 2010-05-19.
Copyright (c) 2010

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com

PURPOSE:  compute the majority tree (> 0.5), given a list of trees on input.  
If inputting mpest trees, will also take a control file and label the leaves
given the control file contents.

USAGE:  python ../../consensus_tree.py --input=917_loci_19_species_mpest.tree \
            --control-file=917_loci_19_species.control 
            --output=917_loci_19_species_mpest_consensus.tree
"""

import os
import sys
import glob
import optparse
import dendropy

import pdb


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the configuration file.', 
metavar='FILE')
    p.add_option('--control-file', dest = 'control', action='store', 
type='string', default = None, help='The path to the configuration file.', 
metavar='FILE')
    p.add_option('--output', dest = 'output', action='store', 
type='string', default = None, help='The path to the output file.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    options.input  = os.path.abspath(os.path.expanduser(options.input))
    options.output  = os.path.abspath(os.path.expanduser(options.output))
    if not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.input):
        print "You must provide a valid path to the configuration file."
        p.print_help()
        sys.exit(2)
    return options, arg


def make_consensus(treefile, min_freq=0.5):
    """generatae a consensus tree from the input"""
    trees = dendropy.TreeList()
    for tree_file in [treefile]:
        trees.read_from_path(tree_file, "newick")
    con_tree = trees.consensus(min_freq)
    return con_tree


def newick2nexus(in_dir):
    """some doc"""
    alignments = glob.glob(os.path.join(in_dir, '*.tre'))
    phyml_trees = dendropy.TreeList()
    for tree_file in alignments:
        phyml_trees.read_from_path(tree_file, 'newick',)
    return phyml_trees


def get_species_dict(control):
    """docstring for get_species_dict"""
    control_file = open(control, 'rU')
    [control_file.readline() for i in xrange(3)]
    species = []
    for line in control_file.readlines():
        ls = line.strip().split('\t')
        if ls[0] == '0':
            break
        else:
            species.append(ls[0])
    return dict(zip(xrange(1, len(species) + 1), species))


def main():
    options, arg = interface()
    if options.control:
        # get a dict of the species in the file
        sp_dict = get_species_dict(options.control)
        cons = make_consensus(options.input)
        # rename the leaves with something other than 1,2,3
        for leaf in cons.leaf_nodes():
            current = leaf.taxon.label
            new = sp_dict[int(current)]
            leaf.taxon.label = new
    else:
        cons = make_consensus(options.input)
    if options.output:
        outf = open(options.output, 'w')
        outf.write(cons.as_string('newick'))
        outf.close()
    else:
        print cons.as_string('newick')


if __name__ == '__main__':
    main()
