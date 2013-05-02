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
import re
import sys
import gzip
import glob
import argparse
import dendropy
import textwrap
from operator import itemgetter
import rpy2.robjects as robjects
from multiprocessing import Pool

import pdb

def interface():
    """ create commandline interface for script"""
    
    description="""
    
Phybase.py calculates species trees from gene trees.
Phybase.py is actually a wrapper script that runs an R package
of the same name (Liu & Yu 2010).

Dependancies:
-------------

Phybase R package: basic functions for phylogenetic analysis
    Installation: 
        At the R prompt type 'install.packages("Phybase")'
    Website: http://cran.r-project.org/web/packages/phybase/index.html

DendroPy: phylogenetic computing library:
    Installation: sudo easy_install -U dendropy
    Website: http://packages.python.org/DendroPy/

Rpy2: simple and efficient access to R from Python
    Installation: sudo easy_install -U rpy2
    Website: http://rpy.sourceforge.net/rpy2.html
    
Argparse: present in python 2.7 and later (I think)

References:
-----------
Liu, L., & Yu, L. (2010). Phybase: an R package for species tree 
analysis. Bioinformatics (Oxford, England). 
doi:10.1093/bioinformatics/btq062 
"""
     
    p = argparse.ArgumentParser(description,)
    
    p.add_argument('--input-file','-i',
        help='Path to input file.')
    p.add_argument('--outgroup','-o',
        help='Name of outgroup.')
    p.add_argument('--method','-m',
        choices=['star', 'njst'],
        default='star',
        help='Phybase method to use.')
    p.add_argument('--genetrees','-g', action='store_true',
        help='Set this flag if the input is genetrees and you want the species tree.')
    p.add_argument('--bootstraps','-b', action='store_true',
        help='Set this flag if the input is bootstraps and you want multiple species trees.')
    p.add_argument('--print-taxa','-p', action='store_true',
        help='Print all the taxon names in the first tree.')
    p.add_argument('--sorted','-s', action='store_true',
        help='Set this flag if your bootstraps are sorted by key.')
    p.add_argument('--cores','-c', type = int, default = 1,
        help='The number of compute cores to use')

    args = p.parse_args()
    
    # check options for errors, etc.
    if args.input_file == None:
        print "Input directory required."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if args.genetrees == True and args.bootstraps == True:
        print "You must pick either genetrees or bootstraps,"
        print "but not both."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()

    if args.genetrees == False and args.bootstraps == False and args.sorted == False and args.print_taxa == False:
        print "You must select either --genetrees, --bootstraps, --sorted, or --print-taxa"
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if args.print_taxa != True and args.method == 'star' and args.outgroup == None:
        print 'You much define an outgroup'
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if args.outgroup:
        args.outgroup = args.outgroup.upper()

    return args


def cleanPhybaseTree(tree):
    tree.strip("\"")
    tree = tree.split("\"")
    return tree[1]


def remove_branch_lengths(str_newick_tree):
    tree = dendropy.Tree()
    tree.read_from_string(str_newick_tree, 'newick')
    for nd in tree:
        nd.label = None
    tree = tree.as_newick_string()
    return tree


def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            comma_back_paren_s = 1
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


def cleanPhyMLTree(tree, upper=True):
    tree = tree.strip()
    tree = remove_branch_lengths(tree)
    # removes support values (=total hack)
    tree = dendropy.Tree.get_from_string(tree, 'newick')
    if upper:
        for leaf in tree.leaf_nodes():
            leaf.taxon.label = leaf.taxon.label.upper()
    # lowercase names
    tree = tree.as_newick_string()
    # converts numbers in sci. notation
    # to decimals (e.g., 1e-22 = 0.000000)
    tree = branch_lengths_2_decimals(tree)
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

    def cleanPhybaseTree(self, tree):
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
        star_sptree = self.cleanPhybaseTree(str(star_sptree))
        steac_sptree = self.cleanPhybaseTree(str(steac_sptree))
        return (star_sptree, steac_sptree)
    
    def njst(self, trees, all_taxa):
        """ generate NJst trees from a list of trees. Requires Phybase and rpy2."""
        # get trees
        trees = self.obj.StrVector(trees)
        # get order of taxa based on last 2 chars of name
        ordered = sorted([[i, j[-2:]] for i, j in enumerate(all_taxa)], key=itemgetter(1))
        # re-order all_taxa
        reordered_all_taxa = [all_taxa[i[0]] for i in ordered]
        # conver to vector
        species_taxaname = self.obj.StrVector(reordered_all_taxa)
        # get species names based on taxon names
        spname = []
        for i in ordered:
            if i[1] not in spname:
                spname.append(i[1])
            else:
                pass
        # convert to string vector
        species_spname = self.obj.StrVector(spname)
        # get count of taxa associated with each species
        cnt = [0] * len(spname)
        for i in ordered:
            idx = spname.index(i[1])
            cnt[idx] += 1
        # get sequence representing count
        ss = self.obj.FloatVector(cnt)
        pdb.set_trace()
        species_structure = self.obj.r['spstructure'](ss)
        # run NJst
        njst_sptree = self.obj.r['NJst'](trees, species_spname, species_taxaname, species_structure)
        njst_sptree = self.cleanPhybaseTree(str(njst_sptree))
        return (njst_sptree)


def phyMLTrees(directory):
    """ get and format for phybase() PhyML 3.0 trees in a directory."""
    # this 'taxa_labels portion' corrects a bug where some trees are
    # missing expected taxa. This should be fixed such that it simply
    # checks for an expected number of taxa and no duplicates
    tree_list = []
    for tree_file in glob.glob(os.path.join(directory, '*tree.txt')):
        for tree in open(tree_file,'r'):
            tree = cleanTree(tree)
            tree_list.append(tree)
    return tree_list


def consensus(tree_list, min_freq=0.5):
    trees = dendropy.TreeList()
    for tree in tree_list:
        t = dendropy.Tree()
        t.read_from_string(tree, 'newick')
        trees.append(t)
    con_tree = trees.consensus(min_freq)
    return con_tree.as_string('newick')


def getTaxa(tree):
    t = dendropy.Tree()
    if "=" in tree:
        tree = tree.split('=')[-1].strip()
    t.read_from_string(tree, 'newick', preserve_underscores=True)
    taxa = t.taxon_set.labels()
    return taxa


def parseBootreps(args):  
    bootreps = {}
    fin = open(args.input_file, 'rU')
    taxa = []
    for count, line in enumerate(fin):
        key, tree = line.split("\t")
        tree = cleanPhyMLTree(tree)
        if count == 0:
            taxa = getTaxa(tree)
        if bootreps.has_key(key) != True:
            bootreps[key] = [tree]
        else:
            bootreps[key].append(tree)  
    star_fout = open(os.path.splitext(args.input_file)[0] + '.star.trees', 'w')
    steac_fout = open(os.path.splitext(args.input_file)[0] + '.steac.trees', 'w')
    steac_trees = []
    star_trees = []
    print 'starting bootreps'
    for count, key in enumerate(bootreps.keys()):
        trees = bootreps[key]
        star_tree, steac_tree = phybase(trees, args.outgroup, taxa)
        steac_trees.append(steac_tree)
        star_trees.append(star_tree)
        star_fout.write(star_tree)
        steac_fout.write(steac_tree)
        print 'processed', count
    steac_consensus = consensus(steac_trees)
    star_consensus = consensus(star_trees)
    template = """#NEXUS\n
                    begin trees;\n
                    tree 'STARConsensus' = %s\n
                    tree 'STEACConsensus' = %s\n
                    end;\n""" % (star_consensus, steac_consensus)
    template = textwrap.dedent(template)
    steac_star_cons_out = os.path.splitext(args.input_file)[0]
    steac_star_cons_out = os.path.join(steac_star_cons_out, 'steac_star.consensus.trees')
    steac_star_cons_out.write(template)


def get_file_chunks(args):
    # LOOP THROUGH SORTED FILE
    if os.path.splitext(args.input_file)[-1] == '.gz':
        f = gzip.open(args.input_file)
    else:
        f = open(args.input_file)
    # if gzip, skip a line
    if os.path.splitext(args.input_file)[-1] == '.gz':
        f.readline()
    taxa = None
    while 1:
        start = f.tell()
        line = f.readline()
        if not line:
            break
        bootrep, tree = line.split("\t")
        if not taxa:
            tree = cleanPhyMLTree(tree)
            taxa = getTaxa(tree)
        while line.split("\t")[0] == bootrep:
            line = f.readline()
        else:
            #pdb.set_trace()
            f.seek(-len(line), 1)
            yield start, f.tell() - start, args.input_file, args.outgroup, taxa, args.method
    f.close()


def get_chunk_data(chunk):
    f = open(chunk[2])
    f.seek(chunk[0])
    if chunk[1]:
        data = f.read(chunk[1])
    else:
        data = f.read()
    data = data.strip().split("\n")
    if data != ['']:
        return data


def worker(x):
    # instantiate Phybase
    phybase = Phybase()
    data = get_chunk_data(x)
    iden = data[0].split("\t")[0]
    trees = [cleanPhyMLTree(d.split("\t")[1]) for d in data]
    if x[5] == 'star':
        results = phybase.run(trees, x[3], x[4])
    elif x[5] == 'njst':
        results = phybase.njst(trees, x[4])
    print 'processed', len(trees), \
        'trees of bootstrap replicate', iden
    return results


def parseSortedBootreps(args):
    """docstring for parseBootreps"""
    chunks = get_file_chunks(args)
    if args.cores > 1:
        p = Pool(args.cores)
        results = p.map(worker, chunks)
	p.close()
    else:
        results = map(worker, chunks)
    if args.method == 'star':
        # SETUP OUTPUT FILES
        star_file = os.path.splitext(args.input_file)[0]
        star_file += '.star.trees'
        star_fout = open(star_file, 'w')
        steac_file = os.path.splitext(args.input_file)[0]
        steac_file += '.steac.trees'
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
        template = """#NEXUS\nbegin trees;\ntree 'STARConsensus' = {0}\ntree 'STEACConsensus' = {1}\nend;""".format(consensus(all_star), consensus(all_steac))
        steac_star_cons_out = os.path.splitext(args.input_file)[0]
        steac_star_cons_out += '.steac_star.consensus.trees'
        steac_star_cons_out = open(steac_star_cons_out, 'w')
        steac_star_cons_out.write(template)
        steac_star_cons_out.close()
    elif args.method == 'njst':
        # SETUP OUTPUT FILES
        njst_file = os.path.splitext(args.input_file)[0]
        njst_file += '.njst.trees'
        njst_fout = open(njst_file, 'w')
        all_njst = []
        pdb.set_trace()
        for tree in results:
            all_njst.append(tree)
            njst_fout.write("{}\n".format(tree))
        # Clean up files
        njst_fout.close()
        template = """#NEXUS\nbegin trees;\ntree 'NJstConsensus' = {0}\nend;""".format(consensus(all_njst))
        njst_cons_out = os.path.splitext(args.input_file)[0]
        njst_cons_out += '.njst.consensus.trees'
        njst_cons_out = open(njst_cons_out, 'w')
        njst_cons_out.write(template)
        njst_cons_out.close()

def get_genetree_chunks(args, is_nexus):
    """return generator from genetree input file where elements are trees"""
    if os.path.splitext(args.input_file)[-1] == '.gz':
        f = gzip.open(args.input_file)
    else:
        f = open(args.input_file)
    while 1:
        start = f.tell()
        line = f.readline()
        if not line:
            break
        yield start, f.tell() - start, args.input_file, is_nexus
    f.close()


def clean_genetree_worker(tree_offsets):
    """clean input genetree and return cleaned tree"""
    sys.stdout.write('.')
    sys.stdout.flush()
    line = get_chunk_data(tree_offsets)[0]
    is_nexus = tree_offsets[3]
    if is_nexus == True:
        if len(line.strip().split("=")) == 2:
            tree = line.strip().split("=")[-1]
            tree = tree.strip()
    else:
        tree = line.strip()
        if "=" in tree:
            tree = tree.split('=')[-1].strip()
    tree = tree.strip(";")
    tree = cleanPhyMLTree(tree)
    return tree


def parse_genetrees(args):
    """parse a set of genetrees in serial or parallel fashion and run through PHYBASE"""
    is_nexus = False
    if args.input_file.endswith('.nex') or args.input_file.endswith('.nexus'):
        is_nexus = True
    chunks = get_genetree_chunks(args, is_nexus)
    print "Cleaning genetrees"
    if args.cores > 1:
        p = Pool(args.cores)
        trees = p.map(clean_genetree_worker, chunks)
	p.close()
    else:
        trees = map(clean_genetree_worker, chunks)
    # get taxa from first tree
    taxa = getTaxa(trees[0])
    # instantiate Phybase instance and analyse trees
    phybase = Phybase()
    if args.method == 'star':
        star_tree, steac_tree = phybase.run(trees, args.outgroup, taxa)
        template = """\n#NEXUS\nbegin trees;\ntree 'STAR' = %s\ntree 'STEAC' = %s\nend;""" % (star_tree, steac_tree)
        print template
        star_steac_out = os.path.splitext(args.input_file)[0]
        star_steac_out += '.star_steac.trees'
        star_steac_out = open(star_steac_out, 'w')
        star_steac_out.write(template)
        star_steac_out.close()
    elif args.method == 'njst':
        njst_tree = phybase.njst(trees, taxa)
        template = """\n#NEXUS\nbegin trees;\ntree 'Njst' = %s\nend;""" % (njst_tree)
        print template
        njst_out = os.path.splitext(args.input_file)[0]
        njst_out += '.njst.trees'
        njst_out = open(njst_out, 'w')
        njst_out.write(template)
        njst_out.close()

def print_taxa(args):
    if os.path.splitext(args.input_file)[-1] == '.gz':
        fin = gzip.open(args.input_file, 'r')
    else:
        fin = open(args.input_file, 'rU')
    line = fin.readline()
    for count, line in enumerate(fin):
        if count == 1:
            tree = line.split('\t')[-1]
            taxa = getTaxa(tree)
            taxa.sort()
            # remove
            taxa.pop(0)
            for taxon_count, taxon in enumerate(taxa, 1):
                print taxon
            print "---------\n", taxon_count, 'total taxa.'
            break


def main():
    args = interface()
    if args.print_taxa == True:
        print_taxa(args)
        sys.exit()
    if args.bootstraps and not args.sorted:
        parseBootreps(args)
    elif args.genetrees:
        parse_genetrees(args)
    elif args.bootstraps and args.sorted:
        parseSortedBootreps(args)

if __name__ == '__main__':
    main()
