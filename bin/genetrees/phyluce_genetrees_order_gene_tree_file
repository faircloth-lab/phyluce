#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 01 June 2015 11:40 CDT (-0500)
"""

import dendropy
import argparse
import ConfigParser
from phyluce.helpers import is_file, is_dir, FullPaths
from phyluce.log import setup_logging

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Root some genetrees"""
    )
    parser.add_argument(
        "--treefile",
        required=True,
        action=FullPaths,
        type=is_file,
        help="""The treefile to input"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        help="""The output treefile to create"""
    )
    parser.add_argument(
        "--input-format",
        choices=["nexus", "newick", "nexml", "fasta", "phylip"],
        default="newick",
        help="""The format of the input data"""
    )
    parser.add_argument(
        "--output-format",
        choices=["nexus", "newick", "nexml", "fasta", "phylip"],
        default="newick",
        help="""The format of the input data"""
    )
    parser.add_argument(
        "--root",
        required=True,
        type=str,
        help="""The taxon on which to root trees"""
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
    return parser.parse_args()


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # read in trees
    log.info("Reading the treefile...")
    tree_list = dendropy.TreeList.get_from_path(args.treefile, schema=args.input_format)
    log.info("Found {} trees.".format(len(tree_list)))
    # replace underscores in outgroup names
    root_taxon = args.root.replace("_", " ")
    trees_to_keep = []
    drop_count = 0
    log.info("Rerooting trees and sorting by count (high to low)...")
    for tree in tree_list:
        # get taxon list for trees - dont use treeset because taxon names include
        # all names in set
        tree_taxa = set([node.label for node in [node.taxon for node in tree] if node is not None])
        # ensure outgroup is in taxon list
        if root_taxon in tree_taxa:
            mrca = tree.mrca(taxon_labels=[root_taxon])
            tree.reroot_at_edge(mrca.edge, update_splits=False)
            #outgroup = tree.find_node_with_taxon_label(root_taxon)
            #tree.to_outgroup_position(outgroup, update_splits=False)
            trees_to_keep.append([tree, len(tree_taxa)])
        else:
            drop_count += 1
    trees_to_keep = sorted(trees_to_keep, key=lambda x: x[1])[::-1]
    sorted_trees_to_keep = dendropy.TreeList()
    for tree in trees_to_keep:
        sorted_trees_to_keep.append(tree[0])
    log.info("Dropped {} trees.  Retained {} trees".format(drop_count, len(sorted_trees_to_keep)))
    log.info("Writing trees...")
    sorted_trees_to_keep.write_to_path(args.output, args.output_format, preserve_spaces=False)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))


if __name__ == '__main__':
    main()
