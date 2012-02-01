#!/usr/bin/env python
# encoding: utf-8
"""
nexus_to_concat.py

Created by Brant Faircloth on 2011-03-05.
Copyright (c) 2011 Brant Faircloth. All rights reserved.
"""

import os
import sys
import pdb
import glob
import optparse
import cPickle
from Bio.Nexus import Nexus
from collections import OrderedDict

def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--models', dest = 'models', action='store', 
type='string', default = None, help='The path to the models file.', 
metavar='FILE')

    p.add_option('--aligns', dest = 'aligns', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    p.add_option('--concat', dest = 'concat', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    p.add_option('--metadata', dest = 'metadata', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    p.add_option('--mr-bayes', dest = 'mrbayes', action='store_true', 
default=False, help='[Optional] Format output for MrBayes.')
    
    p.add_option('--fully', dest = 'fully', action='store_true', 
default=False, help='Fully partition the data (or partition by model).')

    p.add_option('--interleave', dest = 'interleave', action='store_true', 
default=False, help='[Optional] Interleave sequence in nexus.')

    p.add_option('--unlink', dest = 'unlink', action='store_true', 
default=False, help='[Optional] Unlink the models.')

    (options,arg) = p.parse_args()
    if not options.models or not options.aligns:
        p.print_help()
        sys.exit(2)
    return options, arg

def get_loci_and_models(infile):
    loci = OrderedDict()
    for line in open(infile, 'rU'):
        ls = line.strip().split('\t')
        group = loci.setdefault(ls[1], OrderedDict())
        group[ls[0]] = None
    return loci

def save_concat_metadata(metadata, outfile):
    o = open(outfile, 'w')
    cPickle.dump(metadata, o)
    o.close()

def add_mr_bayes_params(metadata, outfile, partition_fully, partition_name = 'fully', unlink = False):
    o = open(outfile, 'a')
    o.write('begin mrbayes;\n')
    #o.write('\tset autoclose=yes nowarn=yes;\n\texecute test.nex;\n')
    partitions = []
    params = OrderedDict()
    c = 1
    for model in metadata:
        short_name = model.split('-')[1]
        params[short_name] = []
        if partition_fully:
            for locus in metadata[model]:
                o.write('\tcharset {0} = {1} - {2};\n'.format(locus, metadata[model][locus][0], metadata[model][locus][1]))
                partitions.append(locus)
                params[short_name].append(str(c))
                c += 1
        else:
            #pdb.set_trace()
            o.write('\tcharset {0} = {1} - {2};\n'.format(short_name, metadata[model][0], metadata[model][1]))
            partitions.append(short_name)
            params[short_name].append(str(c))
            c += 1
    o.write('\tpartition {0} = {1}: {2};\n'.format(partition_name, len(partitions), ', '.join(partitions)))
    o.write('\tset partition = {0};\n'.format(partition_name))
    for model in params:
        m = SUBS[model]
        if m['rates']:
            #pdb.set_trace()
            lset = "\tlset applyto=({0}) nst={1} rates={2};\n".format(','.join(params[model]), m['nst'], m['rates'])
        else:
            lset = "\tlset applyto=({0}) nst={1};\n".format(','.join(params[model]), m['nst'])
        o.write(lset)
        if m['statefreqpr']:
            prset = "\tprset applyto=({0}) statefreqpr={1};\n".format(','.join(params[model]), m['statefreqpr'])
            o.write(prset)
    if unlink:
        o.write('\tunlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);\n')
    o.write('end;')
    o.close()

def fully_partition(metadata, aligns):
    to_combine = []
    start = 1
    for model in metadata:
        for locus in metadata[model]:
            nex = Nexus.Nexus(open(os.path.join(aligns, "{0}.nex".format(locus))))
            end = start + nex.nchar - 1
            metadata[model][locus] = (start, end)
            to_combine.append((locus, nex))
            start = end + 1
    combined = Nexus.combine(to_combine)
    #pdb.set_trace()
    return combined, metadata

def model_partition(metadata, aligns):
    to_combine = []
    start = 1
    end = 0
    new_metadata = OrderedDict()
    for model in metadata:
        for locus in metadata[model]:
            nex = Nexus.Nexus(open(os.path.join(aligns, "{0}.nex".format(locus))))
            #s = sum([1 if 'copy' in n else 0 for n in nex.get_original_taxon_order()])
            #if s > 0:
            #    pdb.set_trace()
            end += nex.nchar
            to_combine.append((locus, nex))
        new_metadata[model] = (start, end)
        start = end + 1
    combined = Nexus.combine(to_combine)
    #pdb.set_trace()
    return combined, new_metadata

def main():
    options, args = interface()
    metadata = get_loci_and_models(options.models)
    #concat, metadata = fully_partition(metadata, options.aligns)
    if options.fully:
        concat, metadata = fully_partition(metadata, options.aligns)
    else:
        concat, metadata = model_partition(metadata, options.aligns)
    concat.write_nexus_data(filename=options.concat, interleave=options.interleave, append_sets = False)
    if not options.mrbayes:
        save_concat_metadata(metadata, options.metadata)
    else:
        if options.fully:
            add_mr_bayes_params(metadata, options.concat, options.fully, partition_name = 'fully', unlink = options.unlink)
        else:
            add_mr_bayes_params(metadata, options.concat, options.fully, partition_name = 'partial', unlink = options.unlink)
            
    #pdb.set_trace()
    
SUBS = {
    'GTR':{'nst':6, 'rates':None, 'statefreqpr':None},
    'GTRI':{'nst':6, 'rates':'propinv', 'statefreqpr':None}, 
    'GTRG':{'nst':6, 'rates':'gamma', 'statefreqpr':None}, 
    'GTRIG':{'nst':6, 'rates':'invgamma', 'statefreqpr':None},
    'SYM':{'nst':6, 'rates': None, 'statefreqpr':'fixed(equal)'}, 
    'SYMI':{'nst':6, 'rates': 'propinv', 'statefreqpr':'fixed(equal)'},
    'SYMG':{'nst':6, 'rates': 'gamma', 'statefreqpr':'fixed(equal)'},
    'SYMIG':{'nst':6, 'rates': 'invgamma', 'statefreqpr':'fixed(equal)'},
    'HKY':{'nst':2, 'rates': None, 'statefreqpr':None},
    'HKYI':{'nst':2, 'rates': 'propinv', 'statefreqpr':None},
    'HKYG':{'nst':2, 'rates': 'gamma', 'statefreqpr':None},
    'HKYIG':{'nst':2, 'rates': 'invgamma', 'statefreqpr':None},
    'K2P':{'nst':2, 'rates': None, 'statefreqpr':'fixed(equal)'},
    'K2PI':{'nst':2, 'rates': 'propinv', 'statefreqpr':'fixed(equal)'},
    'K2PG':{'nst':2, 'rates': 'gamma', 'statefreqpr':'fixed(equal)'},
    'K2PIG':{'nst':2, 'rates': 'invgamma', 'statefreqpr':'fixed(equal)'},
    'F81':{'nst':1, 'rates': None, 'statefreqpr':None},
    'F81I':{'nst':1, 'rates': 'propinv', 'statefreqpr':None},
    'F81G':{'nst':1, 'rates': 'gamma', 'statefreqpr':None},
    'F81IG':{'nst':1, 'rates': 'invgamma', 'statefreqpr':None},
    'JC69':{'nst':1, 'rates': None, 'statefreqpr':'fixed(equal)'},
    'JC69I':{'nst':1, 'rates': 'propinv', 'statefreqpr':'fixed(equal)'},
    'JC69G':{'nst':1, 'rates': 'gamma', 'statefreqpr':'fixed(equal)'},
    'JC69IG':{'nst':1, 'rates': 'invgamma', 'statefreqpr':'fixed(equal)'}
}

if __name__ == '__main__':
    main()
