#!/usr/bin/env python
# encoding: utf-8
"""
run_lastz.py

Created by Brant Faircloth on 2010-02-24.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

# Description

A helper script to run lastz.

"""

import sys
import os
import time
import optparse
from phyluce.many_lastz import multi_lastz_runner

import pdb

def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"
    
    p = optparse.OptionParser(usage)
    
    p.add_option('--target', dest = 'target', action='store', \
type='string', default = None, help='The path to the target file (2bit)', \
metavar='FILE')
    p.add_option('--query', dest = 'query', action='store', \
type='string', default = None, help='The path to the query file (2bit)', \
metavar='FILE')
    p.add_option('--output', dest = 'output', action='store', \
type='string', default = None, help='The path to the output file', \
metavar='FILE')
    p.add_option('--nprocs', dest = 'nprocs', action='store', \
type='int', default = 1, help='The number of processors to use')
    p.add_option('--coverage', dest = 'coverage', action='store', \
type='float', default = 83, help='The fraction of bases in the \
entire input sequence (target or query, whichever is shorter) that are \
included in the alignment block, expressed as a percentage')
    p.add_option('--identity', dest = 'identity', action='store', \
type='float', default = 92.5, help='The fraction of aligned bases \
(excluding columns containing gaps or non-ACGT characters) that are \
matches, expressed as a percentage')
    p.add_option('--huge', dest = 'huge', action='store_true', default=False, \
help='Deal with poorly assembled (many scaffolds) genome sequences')
    p.add_option('--size', dest = 'size', action='store', \
type='int', default = 10000000, help='The chunk size (in bp) to stick in a \
file while using the --huge option')
    
    (options,arg) = p.parse_args()
    for f in [options.target, options.query, options.output]:
        if not f:
            p.print_help()
            sys.exit(2)
        if f != options.output and not os.path.isfile(f):
            print "You must provide a valid path to the query/target file."
            p.print_help()
            sys.exit(2)
    return options, arg


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    options, arg = interface()
    multi_lastz_runner(
            options.output,
            options.nprocs,
            options.target,
            options.query,
            options.huge,
            options.coverage,
            options.identity,
            options.size
        )
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time) / 60, 'minutes'

if __name__ == '__main__':
    main()

