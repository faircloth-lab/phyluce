#!/usr/bin/env python
# encoding: utf-8

"""
run_mpest.py

Created by Brant Faircloth on 16 September 2010 09:37 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  Run an input file through our slightly customized version of mpest
either using a single core or multiple cores, and feeding mpest a random 
(integer) seed drawn from a uniform distribution.

USAGE:  python ../../run_mpest.py --control-file=408_loci_25_species.control \
            --iterations=1000 --cores=7 
            --output=408_loci_25_species_mpest.tree

"""

import pdb
import sys
import os
import time
import numpy
import optparse
import tempfile
import subprocess
import multiprocessing


def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"
    
    p = optparse.OptionParser(usage)
    
    p.add_option('--control-file', dest = 'input', action='store', \
type='string', default = None, help='The path to the input control file', \
metavar='FILE')
    p.add_option('--output', dest = 'output', action='store', \
type='string', default = None, help='The path to the output file', \
metavar='FILE')
    p.add_option('--iterations', dest = 'iterations', action='store', \
type='int', default = None, help='The number of iterations to run')
    p.add_option('--cores', dest = 'nprocs', action='store', \
type='int', default = None, help='The number of cores to use')
    (options,arg) = p.parse_args()
    return options, arg


def mpest_cli(input, seed, outfile):
    cli = 'mpest {0} {1} {2}'.format(input, seed, outfile)
    return cli


def single_mpest(input, iterations):
    '''docstring'''
    seeds = numpy.random.random_integers(100000, high=None, size=iterations)
    #pdb.set_trace()
    trees = []
    for seed in seeds:
        temp_fd, temp_out = tempfile.mkstemp(suffix='.mpestout', dir='/tmp')
        os.close(temp_fd)
        cli = mpest_cli(input, seed, temp_out)
        print cli
        #pdb.set_trace()
        mpest_out, mpest_stderr = subprocess.Popen(cli, shell=True, 
            stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        #print mpest_stderr
        #pdb.set_trace()
        for line in reversed(open(temp_out).readlines()):
            if 'mpestTree' in line:
                pre, post = line.strip('\n').split('=')
                post = post.strip(' ')
                post = post.rstrip(';')
                trees.append(post)
        os.remove(temp_out)
    return trees

def multi_mpest(input, output, control_file):
    '''docstring'''
    for seed in iter(input.get, 'STOP'):
        temp_fd, temp_out = tempfile.mkstemp(suffix='.mpestout', dir='/tmp')
        os.close(temp_fd)
        cli = mpest_cli(control_file, seed, temp_out)
        print cli
        #pdb.set_trace()
        mpest_out, mpest_stderr = subprocess.Popen(cli, shell=True, 
            stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        #pdb.set_trace()
        for line in reversed(open(temp_out).readlines()):
            if 'mpestTree' in line:
                pre, post = line.strip('\n').split('=')
                post = post.strip(' ')
                post = post.rstrip(';')
                output.put(post)
        os.remove(temp_out)


def q_runner(n_procs, list_item, function, *args):
    '''generic function used to start worker processes'''
    task_queue      = multiprocessing.Queue()
    results_queue   = multiprocessing.JoinableQueue()
    if args:
        arguments = (task_queue, results_queue,) + args
    else:
        arguments = (task_queue, results_queue,)
    results = []
    # reduce processer count if proc count > files
    if len(list_item) < n_procs:
        n_procs = len(list_item)
    for l in list_item:
        task_queue.put(l)
    for _ in range(n_procs):
        p = multiprocessing.Process(target = function, args = arguments).start()
        #print 'Starting %s' % function
    for _ in range(len(list_item)):
        # indicated done results processing
        results.append(results_queue.get())
        results_queue.task_done()
    #tell child processes to stop
    for _ in range(n_procs):
        task_queue.put('STOP')
    # join the queue until we're finished processing results
    results_queue.join()
    # not closing the Queues caused me untold heartache and suffering
    task_queue.close()
    results_queue.close()
    return results


def main():
    start_time      = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    options, arg = interface()
    if options.nprocs == 1:
        trees = single_mpest(options.input, options.iterations)
        outp = open(options.output, 'w')
        outp.write(';\n'.join(trees))
        outp.write(';')
        outp.close()
    else:
        seeds = map(None, numpy.random.random_integers(100000, high=None, size=options.iterations))
        #pdb.set_trace()
        trees = q_runner(options.nprocs, seeds, multi_mpest, options.input)
        outp = open(options.output, 'w')
        outp.write(';\n'.join(trees))
        outp.write(';')
        outp.close()
        #pdb.set_trace()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time)/60, 'minutes'



if __name__ == '__main__':
    main()
