#!/usr/bin/env python
# encoding: utf-8
"""
run_lastz.py

Created by Brant Faircloth on 2010-02-24.
Copyright (c) 2010 Brant Faircloth. All rights reserved.

# Description

A helper script to run lastz.

"""

import pdb
import sys
import os
import time
import optparse
import tempfile
import subprocess
import bx.seq.twobit
import multiprocessing



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


def lastzParams(chromo, probe, coverage, identity, temp_out):
    cli = \
    'lastz \
    %s \
    %s[nameparse=full]\
    --strand=both \
    --seed=12of19 \
    --transition \
    --nogfextend \
    --nochain \
    --gap=400,30 \
    --xdrop=910 \
    --ydrop=8370 \
    --hspthresh=3000 \
    --gappedthresh=3000 \
    --noentropy \
    --coverage=%s \
    --identity=%s \
    --output=%s \
    --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity,coverage' \
    % \
    (chromo, probe, coverage, identity, temp_out)
    return cli

def lastz(input, output, coverage, identity):
    '''docstring for worker2'''
    for chromo, probe in iter(input.get, 'STOP'):
        print '\t%s' % chromo
        temp_fd, temp_out = tempfile.mkstemp(suffix='.lastz')
        os.close(temp_fd)
        cli = lastzParams(chromo, probe, coverage, identity, temp_out)
        lzstdout, lztstderr = subprocess.Popen(cli, shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        if lztstderr:
            output.put(lztstderr)
        else:
            output.put(temp_out)

def SingleProcLastz(input, output, coverage, identity):
    '''docstring for worker2'''
    #pdb.set_trace()
    chromo, probe = input
    temp_fd, temp_out = tempfile.mkstemp(suffix='.lastz')
    os.close(temp_fd)
    cli = lastzParams(chromo, probe, coverage, identity, temp_out)
    #pdb.set_trace()
    lzstdout, lztstderr = subprocess.Popen(cli, shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
    if lztstderr:
        output.append(lztstderr)
    else:
        output.append(tmp_out)
    return output


def main():
    start_time      = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    options, arg = interface()
    if not options.huge:
        # get individual records from the 2bit file
        chromos =  [os.path.join(options.target, c) for c in bx.seq.twobit.TwoBitFile(file(options.target)).keys()]
    else:
        chromos = []
        # split target file into `options.size` (~10 Mbp) chunks
        temp_fd, temp_out = tempfile.mkstemp(suffix='.fasta')
        os.close(temp_fd)
        temp_out_handle = open(temp_out, 'w')
        tb = bx.seq.twobit.TwoBitFile(file(options.target))
        sequence_length = 0
        print 'Running with the --huge option.  Chunking files into {0} bp...'.format(options.size)
        for seq in tb.keys():
            sequence = tb[seq][0:]
            sequence_length += len(sequence)
            # write it to the outfile
            temp_out_handle.write('>{0}\n{1}\n'.format(seq, sequence))
            if sequence_length > options.size:
                temp_out_handle.close()
                # put tempfile name on stack
                chromos.append(temp_out+'[multiple]')
                # open a new temp file
                temp_fd, temp_out = tempfile.mkstemp(suffix='.fasta')
                os.close(temp_fd)
                temp_out_handle = open(temp_out, 'w')
                # reset sequence length
                sequence_length = 0
    probes = (options.query,) * len(chromos)
    cp = zip(chromos, probes)
    # put those record names on the stack
    print "Running the targets against %s queries..." % len(chromos)
    if options.nprocs == 1:
        results = []
        for each in cp:
            print each
            print results
            results = SingleProcLastz(each, results, options.coverage, options.identity)
    else:
        results = q_runner(options.nprocs, cp, lastz, options.coverage, options.identity)
    outp = open(options.output, 'wb')
    print "Writing the results file..."
    #pdb.set_trace()
    for f in results:
        print '\t%s' % f
        # read the file
        outp.write(open(f,'rb').read())
        # cleanup the lastz output files
        os.remove(f)
    outp.close()
    print 'Cleaning up the chunked files...'
    if options.huge:
        for f in chromos:
            # cleanup the chunked files
            os.remove(f.strip('[multiple]'))
    # stats
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()

