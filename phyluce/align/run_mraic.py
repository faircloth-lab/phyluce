#!/usr/bin/env python
# encoding: utf-8

"""
run_mraic.py

Created by Brant Faircloth on 20 December 2010 14:40 PST (-0800).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import os
import sys
import glob
import shutil
import optparse
import subprocess
import multiprocessing


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to input directory.', 
metavar='FILE')

    p.add_option('--output', dest = 'output', action='store', 
type='string', default = None, help='The path to output directory.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.input or not options.output:
        p.print_help()
        sys.exit(2)
    if not os.path.isdir(options.input) or not os.path.isdir(options.output):
        print "You must provide a valid path to input and/or output directories"
        sys.exit(2)
    return options, arg

def run_mraic(files_to_process):
    sys.stdout.write(".")
    sys.stdout.flush()
    input, output = files_to_process
    cmd = "mraic_mod.pl --infile='{0}' --output_dir='{1}' >/dev/null 2>&1".format(input, output)
    p = subprocess.Popen(cmd, shell=True).communicate()

def main():
    options, arg = interface()
    phylip_files = glob.glob(os.path.join(options.input, '*.phylip'))
    output = [options.output] * len(phylip_files)
    files_to_process = zip(phylip_files, output)
    pool = multiprocessing.Pool(7)
    pool.map(run_mraic, files_to_process)
    #map(run_mraic, files_to_process)
    output_files = glob.glob(os.path.join(options.output, '*.phylip.AICc*'))
    results = {}
    for f in output_files:
        path, name = os.path.split(f)
        name_split = name.split('.')
        results[name_split[0]] = name_split[2]
    outf = open("output_models.txt",'w')
    for k,v in results.iteritems():
        outf.write("{0}\t{1}\n".format(k,v))
    outf.close()


if __name__ == '__main__':
    main()
