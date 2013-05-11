#!/usr/bin/env python
# encoding: utf-8
"""
get_mpest_format.py

Created by Nick Crawford on 2010-05-16.
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

PURPOSE:  Given directory or file input, read the contents of the either and 
root the trees within (using --root and --outgroup= options).  Additionally, 
build a control file for our slightly customized version of mpest.

USAGE:  python ../../get_mpest_format.py --input=408Loci_25Species.PhyML.trees \
            --output=408_loci_25_species.tree \
            --root --outgroup=anoCar2 
            --build-control
"""

import os
import sys
import glob
import optparse
from ete2 import Tree

import pdb

def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"
    
    p = optparse.OptionParser(usage)
    
    p.add_option('--input', dest = 'input', action='store', \
type='string', default = None, help='The path to the input file/directory')
    p.add_option('--output', dest = 'output', action='store', \
type='string', default = None, help='The path to the output file', \
metavar='FILE')
    p.add_option('--outgroup', dest = 'outgroup', action='store', \
type='string', default = None, help='The outgroup', \
metavar='FILE')
    p.add_option('--extension', dest = 'extension', action='store', \
type='string', default = '*.tre', help='The extension of tree files')
    p.add_option('--root', dest = 'root', action='store_true', default=False, 
help='Root trees in input directory')
    p.add_option('--build-control', dest = 'control', action='store_true', default=False, 
help='Build mpest control files for trees in input directory')

    (options,arg) = p.parse_args()
    options.input  = os.path.abspath(os.path.expanduser(options.input))
    options.output = os.path.abspath(os.path.expanduser(options.output))
    return options, arg


def create_rooted_trees_from_dir(paths, fout, outgroup, control):
    """provide paths of phyml bootstrap replicates """
    #pdb.set_trace()
    for count, path in enumerate(paths):
        base_path, tree_file_name = os.path.split(path)
        rooted_tree_name = path.split('.')[0] + '.PhyML.rooted.trees'
        fout = open(rooted_tree_name, 'w')
        #pdb.set_trace()
        fin = open(path)
        for tree in fin:
            tree = tree.strip()
            tree = Tree(tree)
            tree.set_outgroup(outgroup)
            newick = tree.write(format=5) + '\n'
        fout.write(newick)
        print count + 1
        fout.close()
        if control:
            #pdb.set_trace()
            print rooted_tree_name
        create_control_files(rooted_tree_name)


def create_rooted_trees_from_file(input, fout, outgroup):
    """docstring for create_rooted_trees_from_file"""
    fin = open(input, 'rU')
    fout = open(fout, 'w')
    for count, line in enumerate(fin):
        line = line.split('[&U]')[-1].strip()
        tree = Tree(line)
        tree.set_outgroup(outgroup)
        newick = tree.write(format=5) + '\n'
        fout.write(newick)
        print count + 1
    fout.close()
    fin.close()


def create_control_files(path):
    """explaination"""
    #pdb.set_trace()
    # set output filename for control file
    base_path, tree_file = os.path.split(path)
    control_file = '.'.join([tree_file.split('.')[0], 'control'])
    fout = open(os.path.join(base_path, control_file), 'w')
    # open up the input tree file and get the first tree
    tree_file_contents = open(path, 'r').readlines()
    tree = tree_file_contents[0].strip()
    tree = Tree(tree)
    taxa_names = tree.get_leaf_names()
    #pdb.set_trace()
    taxa_string = ['{0}\t1\t{0}'.format(taxa) for taxa in taxa_names]
    template_info = {'filename': path,
                    'tree_count': len(tree_file_contents),
                    'numb_taxa': len(taxa_names),
                    'taxa_details': '\n'.join(taxa_string)}
    # needs to go to format()
    template = "%(filename)s\n0\n%(tree_count)s %(numb_taxa)s\n%(taxa_details)s\n0\n" % template_info
    fout.write(template)
    fout.close()


def main():
    options, arg = interface()
    if options.root:
        # if we're rooting we need an outgroup.  this is now CLI required,
        # so i don't root shit without thinking what the root is
        if not options.outgroup:
            print "You must provide an outgroup (--outgroup=)"
            sys.exit()
        if os.path.isdir(options.input):
            raw_phyml_trees = glob.glob(os.path.join(options.input, '*.phylip_phyml_tree.txt'))
            create_rooted_trees_from_dir(raw_phyml_trees, options.output, options.outgroup, options.control)
            sys.exit()
        elif os.path.isfile(options.input):
            create_rooted_trees_from_file(options.input, options.output, options.outgroup)
    if options.control:
        #pdb.set_trace()
        create_control_files(options.output)
if __name__ == '__main__':
    main()
