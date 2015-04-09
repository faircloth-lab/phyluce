#!/usr/bin/env python
# encoding: utf-8
"""
Phybase.py

Created by Nicholas Crawford and Brant C. Faircloth Copyright (c) 2010 Nicholas Crawford and
Brant C. Faircloth. All rights reserved.

Parses gene trees from directories of tree files. Submits the trees to R and runs phybase to
generate species trees. Both Steac and Star trees are generated. Seperate functions are included
to parse NJ trees from Paup and ML trees from PhyMl as the tree formats are slightly different.

Dependencies:

   dendropy - http://packages.python.org/DendroPy/tutorial/index.html
   rpy2 - http://rpy.sourceforge.net/rpy2.html
   phybase - http://cars.desu.edu/faculty/lliu/research/phybase.html

Future directions:

   - add commandline (optparse)
   - add NJtrees fuction
   - Parsing of trees from Garli
   - Bootstrapping with Phybase
   - Nexus output of Steac and Star trees
      - also png/svg output?
"""

import os
import sys
import gzip
import argparse
import dendropy
import rpy2.robjects as robjects
from multiprocessing import Pool

from phyluce.helpers import FullPaths, is_file

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Phybase.py calculates species trees from gene trees."""
    )
    parser.add_argument(
        "--input",
        required=True,
        action=FullPaths,
        type=is_file,
        help="""The path to the input genetrees or bootreps.  Must be parsed""" +
        """ using reformat_raxml_output.py."""
    )
    parser.add_argument(
        "--outgroup",
        type=str,
        default=None,
        help="""The name of the outgroup on which the trees are rooted.""" +
        """Not needed for NJst"""
    )
    parser.add_argument(
        "--method",
        choices=["star", "njst"],
        default="njst",
        help="""The method to use."""
    )
    parser.add_argument(
        '--cores',
        type = int,
        default = 1,
        help='The number of compute cores to use'
    )
    t = parser.add_mutually_exclusive_group(required=True)
    t.add_argument(
        '--genetrees',
        action='store_true',
        help='Choose this if you want to generate the species tree from gene trees.'
    )
    t.add_argument(
        '--bootreps',
        action='store_true',
        help="""Choose this if you want to generate support values for bootrepped""" +
        """ gene trees"""
    )

    args = parser.parse_args()

    if args.method == 'star' and args.outgroup == None:
        raise ValueError("'star' method requires you to specify an outgroup")
        sys.exit()

    return args


def remove_support_values(str_newick_tree):
    tree = dendropy.Tree()
    tree.read_from_string(str_newick_tree, 'newick')
    for nd in tree:
        nd.label = None
    return tree


def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    #comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            #comma_back_paren_s = 1
            num = '%f' % float(num)
            new_tree += ":" + num
            colon_s = 0
            num = ''
        if colon_s != 0:
            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + ";"
    return new_tree


def clean_phyml_tree(tree, upper=False):
    tree = tree.strip()
    tree = remove_support_values(tree)
    if upper:
        for leaf in tree.leaf_nodes():
            leaf.taxon.label = leaf.taxon.label.upper()
    tree = tree.as_newick_string()
    # converts numbers in sci. notation
    # to decimals (e.g., 1e-22 = 0.000000)
    #tree = branch_lengths_2_decimals(tree)
    return tree


class Blackhole(object):
    """catch output and blackhole it"""

    def write(self, string):
        pass


class Phybase:
    def __init__(self, verbose=False):
        if not verbose:
            # suppress output
            stdout = sys.stdout
            sys.stdout = Blackhole()
        else:
            stdout = sys.stdout
        self.obj = robjects
        # load phybase
        self.obj.r['library']('phybase')
        sys.stdout = stdout
        self._create_nancdist()
        self._create_njst()

    def _create_nancdist(self):
        self.obj.r('''
            nancdist<-function(tree, taxaname)
            {
                ntaxa<-length(taxaname)
                nodematrix<-read.tree.nodes(tree,taxaname)$nodes
                if(is.rootedtree(nodematrix)) nodematrix<-unroottree(nodematrix)
                dist<-matrix(0, ntaxa,ntaxa)
                for(i in 1:(ntaxa-1))
                    for(j in (i+1):ntaxa)
                    {
                    anc1<-ancestor(i,nodematrix)
                    anc2<-ancestor(j,nodematrix)
                    n<-sum(which(t(matrix(rep(anc1,length(anc2)),ncol=length(anc2)))-anc2==0, arr.ind=TRUE)[1,])-3
                    if(n==-1) n<-0
                    dist[i,j]<-n
                    }
                    dist<-dist+t(dist)
                z<-list(dist=as.matrix, taxaname=as.vector)
                z$dist<-dist
                z$taxaname<-taxaname
                z
            }
        ''')

    def _create_njst(self):
        self.obj.r('''
        NJst<-function(genetrees, spname, taxaname, species.structure)
        {
            ntree<-length(genetrees)
            ntaxa<-length(taxaname)
            dist <- matrix(0, nrow = ntree, ncol = ntaxa * ntaxa)
            for(i in 1:ntree)
            {
                genetree1 <- read.tree.nodes(genetrees[i])
                thistreetaxa <- genetree1$names
                ntaxaofthistree <- length(thistreetaxa)
                thistreenode <- rep(-1, ntaxaofthistree)
                dist1<-matrix(0,ntaxa,ntaxa)
                for (j in 1:ntaxaofthistree)
                {
                    thistreenode[j] <- which(taxaname == thistreetaxa[j])
                    if (length(thistreenode[j]) == 0)
                    {
                        print(paste("wrong taxaname", thistreetaxa[j],"in gene", i))
                        return(0)
                    }
                }
                dist1[thistreenode, thistreenode]<-nancdist(genetrees[i],thistreetaxa)$dist
                dist[i,]<-as.numeric(dist1)
            }
            dist[dist == 0] <- NA
            dist2 <- matrix(apply(dist, 2, mean, na.rm = TRUE), ntaxa, ntaxa)
            diag(dist2) <- 0
            if (sum(is.nan(dist2)) > 0)
            {
                print("missing species!")
                dist2[is.nan(dist2)] <- 10000
            }
            speciesdistance <- pair.dist.mulseq(dist2, species.structure)
            tree<-write.tree(nj(speciesdistance))
            node2name(tree,name=spname)
        }
        ''')

    def clean_phybase_tree(self, tree):
        tree.strip("\"")
        tree = tree.split("\"")
        return tree[1]

    def run(self, trees, outgroup, all_taxa):
        """ generate Steac and Star trees from a list of trees. Requires Phybase and rpy2."""
        trees = self.obj.StrVector(trees)
        species_taxaname = self.obj.StrVector(all_taxa)
        # list of species in current tree
        species_spname = species_taxaname
        matrix_size = len(species_taxaname)
        species_structure = self.obj.r['diag'](1, matrix_size, matrix_size)
        star_sptree = self.obj.r['star.sptree'](trees, species_spname, species_taxaname,\
                                                species_structure, outgroup, 'nj')
        steac_sptree = self.obj.r['steac.sptree'](trees, species_spname, species_taxaname,\
                                                species_structure, outgroup, 'nj')
        star_sptree = self.clean_phybase_tree(str(star_sptree))
        steac_sptree = self.clean_phybase_tree(str(steac_sptree))
        return (star_sptree, steac_sptree)

    def njst(self, trees, all_taxa):
        """ generate NJst trees from a list of trees. Requires Phybase and rpy2."""
        # get trees
        trees = self.obj.StrVector(trees)
        species_taxaname = self.obj.StrVector(all_taxa)
        # list of species in current tree
        species_spname = species_taxaname
        matrix_size = len(species_taxaname)
        species_structure = self.obj.r['diag'](1, matrix_size, matrix_size)
        # run NJst
        njst_sptree = self.obj.r['NJst'](trees, species_spname, species_taxaname, species_structure)
        njst_sptree = self.clean_phybase_tree(str(njst_sptree))
        return (njst_sptree)


def consensus(tree_list, min_freq=0.5):
    trees = dendropy.TreeList()
    for tree in tree_list:
        t = dendropy.Tree()
        t.read_from_string(tree, 'newick')
        trees.append(t)
    con_tree = trees.consensus(min_freq)
    return con_tree.as_string('newick')


def get_taxa(trees):
    # need to iterate over all trees to ensure we have full list of taxon names
    taxa = set()
    for tree in trees:
        t = dendropy.Tree()
        if "=" in tree:
            tree = tree.split('=')[-1].strip()
        t.read_from_string(tree, 'newick', preserve_underscores=True)
        taxa.update(t.taxon_set.labels())
    return list(taxa)


def get_file_chunks(args):
    # LOOP THROUGH SORTED FILE
    if os.path.splitext(args.input)[-1] == '.gz':
        f = gzip.open(args.input)
    else:
        f = open(args.input)
    while 1:
        start = f.tell()
        line = f.readline()
        if not line:
            break
        bootrep, tree = line.split("\t")
        while line.split("\t")[0] == bootrep:
            line = f.readline()
        else:
            f.seek(-len(line), 1)
            yield start, f.tell() - start, args.input, args.outgroup, args.method
    f.close()


def get_chunk_data(start, stop, input):
    f = open(input)
    f.seek(start)
    if stop:
        data = f.read(stop)
    else:
        data = f.read()
    data = data.strip().split("\n")
    if data != ['']:
        return data


def worker(work):
    start, stop, input, outgroup, method = work
    data = get_chunk_data(start, stop, input)
    # instantiate Phybase
    phybase = Phybase()
    trees = [clean_phyml_tree(d.split("\t")[1]) for d in data]
    taxa = get_taxa(trees)
    if method == 'star':
        results = phybase.run(trees, outgroup, taxa)
    elif method == 'njst':
        results = phybase.njst(trees, taxa)
    sys.stdout.write(".")
    sys.stdout.flush()
    return results


def parse_bootreps(args):
    """docstring for parse_bootreps"""
    chunks = get_file_chunks(args)
    print "Running"
    if args.cores > 1:
        p = Pool(args.cores)
        results = p.map(worker, chunks)
	p.close()
    else:
        results = map(worker, chunks)
    if args.method == 'star':
        # SETUP OUTPUT FILES
        star_file = os.path.splitext(args.input)[0]
        star_file += '.star.bootrep.trees'
        star_fout = open(star_file, 'w')
        steac_file = os.path.splitext(args.input)[0]
        steac_file += '.steac.bootrep.trees'
        steac_fout = open(steac_file, 'w')
        all_star = []
        all_steac = []
        for trees in results:
            all_star.append(trees[0])
            star_fout.write("{}\n".format(trees[0]))
            all_steac.append(trees[1])
            steac_fout.write("{}\n".format(trees[1]))
        # Clean up files
        star_fout.close()
        steac_fout.close()
        for tree in [(all_star, ".star.bootrep.consensus.tre"), (all_steac, ".steac.bootrep.consensus.tre")]:
            template = "{0}".format(consensus(tree[0]))
            steac_star_cons_out = os.path.splitext(args.input)[0]
            steac_star_cons_out += tree[1]
            steac_star_cons_out = open(steac_star_cons_out, 'w')
            steac_star_cons_out.write(template)
            steac_star_cons_out.close()
    elif args.method == 'njst':
        # SETUP OUTPUT FILES
        njst_file = os.path.splitext(args.input)[0]
        njst_file += '.njst.bootrep.trees'
        njst_fout = open(njst_file, 'w')
        all_njst = []
        for tree in results:
            all_njst.append(tree)
            njst_fout.write("{0}\n".format(tree))
        # Clean up files
        njst_fout.close()
        template = "{0}".format(consensus(all_njst))
        njst_cons_out = os.path.splitext(args.input)[0]
        njst_cons_out += '.njst.bootrep.consensus.trees'
        njst_cons_out = open(njst_cons_out, 'w')
        njst_cons_out.write(template)
        njst_cons_out.close()


def get_genetree_chunks(args, is_nexus):
    """return generator from genetree input file where elements are trees"""
    if os.path.splitext(args.input)[-1] == '.gz':
        f = gzip.open(args.input)
    else:
        f = open(args.input)
    while 1:
        start = f.tell()
        line = f.readline()
        if not line:
            break
        yield start, f.tell() - start, args.input, is_nexus
    f.close()


def clean_genetree_worker(tree_offsets):
    """clean input genetree and return cleaned tree"""
    start, stop, input, is_nexus = tree_offsets
    sys.stdout.write('.')
    sys.stdout.flush()
    line = get_chunk_data(start, stop, input)[0]
    if is_nexus:
        if len(line.strip().split("=")) == 2:
            tree = line.strip().split("=")[-1]
            tree = tree.strip()
    else:
        tree = line.strip()
        if "=" in tree:
            tree = tree.split('=')[-1].strip()
    tree = tree.strip(";")
    tree = clean_phyml_tree(tree)
    return tree


def parse_genetrees(args):
    """parse a set of genetrees in serial or parallel fashion and run through PHYBASE"""
    is_nexus = False
    if args.input.endswith('.nex') or args.input.endswith('.nexus'):
        is_nexus = True
    chunks = get_genetree_chunks(args, is_nexus)
    print "Cleaning genetrees"
    if args.cores > 1:
        p = Pool(args.cores)
        trees = p.map(clean_genetree_worker, chunks)
	p.close()
    else:
        trees = map(clean_genetree_worker, chunks)
    # flush after sys.stdout
    print ""
    # get taxa from first tree
    taxa = get_taxa(trees)
    # instantiate Phybase instance and analyse trees
    phybase = Phybase()
    if args.method == 'star':
        star_tree, steac_tree = phybase.run(trees, args.outgroup, taxa)
        for cnt, tree in enumerate([(star_tree, ".star.species.tre"), (steac_tree, ".steac.species.tre")]):
            template = "{0}".format(tree[0])
            if cnt == 0:
                print "STAR = ", template
            else:
                print "STEAC = ", template
            star_steac_out = os.path.splitext(args.input)[0]
            star_steac_out += tree[1]
            star_steac_out = open(star_steac_out, 'w')
            star_steac_out.write(template)
            star_steac_out.close()
    elif args.method == 'njst':
        njst_tree = phybase.njst(trees, taxa)
        # just output newick
        template = "{0}".format(njst_tree)
        print "NJst = ", template
        njst_out = os.path.splitext(args.input)[0]
        njst_out += '.njst.species.tre'
        njst_out = open(njst_out, 'w')
        njst_out.write(template)
        njst_out.close()


def main():
    args = get_args()
    if args.genetrees:
        parse_genetrees(args)
    elif args.bootreps:
        parse_bootreps(args)

if __name__ == '__main__':
    main()
