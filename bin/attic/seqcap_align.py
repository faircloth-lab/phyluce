#!/usr/bin/env python
# encoding: utf-8
"""
seqcap_align.py

Created by Nicholas Crawford and Brant C. Faircloth
Copyright (c) 2010 Nicholas Crawford and Brant C. Faircloth. All rights reserved.

This script aligns fasta sequence groups on a per locus basis (where the locus
name is in fasta header).  It takes a fasta file of reads, in arbitrary order,
groups reads by locus and uses MUSCLE (http://www.drive5.com/muscle/) to do
align reads by locus.  We use the Biopython Alignment class to hold and
output reads in various formats (fasta, nexus, etc).  We've also implemented
a class (ConcatenatedAlignment) which concatenates alignments by locus, and
outputs one large alignment in various formats.

We've also implemented options to trim aligned reads at specified distances
from an internal probe sequence and also remove the specified sequence
from within reads.

python seqcapAlign.py --input=all.fa --output=output/ \
    --probe-file=../Simulation/SureSelectProbes.fsa --trim-flank=100 \
    --multiprocessing --processors=6 --concatenate

"""

import pdb
import sys
import os
import shutil
import optparse
import tempfile
import multiprocessing
import phyluce.muscle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment


class ConcatenatedAlignment(Alignment):
    '''A child class of the biopython alignment class, created so that we could
    add the concat method to concatenate multiple alignments'''
    def __init__(self):
        Alignment.__init__(self, Gapped(IUPAC.unambiguous_dna, '-'))

    def concat(self, alignment):
        '''concatenate alignment objects (on a per SeqRecord basis)'''
        #pdb.set_trace()
        if not self._records:
            for seq in alignment:
                seq.id = seq.id.split('_')[-1]
                seq.name = seq.id
                self._records.append(seq)
        else:
            for seq in alignment:
                for pos, record in enumerate(self._records):
                    # this assumes that we will be using id as
                    # the joining attribute...
                    if seq.id.split('_')[-1] == record.id:
                        c_seq = SeqRecord(record.seq + seq.seq)
                        c_seq.name = record.id
                        c_seq.id = record.id
                        self._records[pos] = c_seq


class Locus(list):
    '''a subclass of list to hold a group of sequence reads on a per
    locus basis'''
    def __init__(self, *arg):
        list.__init__(self)
        if arg:
            self.append(arg[0])
        self.contents = ''
        self.formatted = False
        self.tempFile = None

    def formatSequences(self):
        '''create the contents of a fasta file on a per locus basis'''
        self.formatted = True
        for item in self:
            self.contents += item.format('fasta')

    def createTempFile(self):
        '''create a tempfile holding the contents of self.contents'''
        if self.contents:
            fd, self.tempFile = tempfile.mkstemp(suffix='.muscle')
            os.write(fd, self.contents)
            os.close(fd)

    #def cleanupTempFile(self):
    #    '''remove the tempfile'''
    #    os.remove(self.tempFile) 


def interfaceWarnings(p, message):
    '''generic warning function for interface options'''
    print message
    p.print_help()
    sys.exit(2)


def interface():
    '''Get the starting parameters from a configuration file'''
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    p.add_option('--input', dest='path', action='store', \
        type='string', default=None, \
        help='The path to the input file.', metavar='FILE')
    p.add_option('--output', dest='output', action='store', \
        type='string', default=None, \
        help='The path to the output directory.', metavar='FILE')
    p.add_option('--format', dest='format', action='store', \
        type='string', default='nexus', \
        help='Output format:  clustal, emboss, fasta, nexus, phylip,'\
            +'stockholm')
    p.add_option('--concatenate', dest='concatenate', action='store_true', \
        default=False, help='Concatenate alignments by locus')
    p.add_option('--probe-file', dest='probefile', action='store', \
        type='string', default=None, \
        help='Path to a fasta file containing the probe sequences.', \
        metavar='FILE')
    p.add_option('--trim-flank', dest='trimflank', action='store', \
        type='int', default=None, \
        help='The length to trim the flank.')
    p.add_option('--species', dest='species', action='store', \
        type='int', default=None, \
        help='Number of species expected in each alignment.')
    p.add_option('--no-ambiguous', dest='no_ambiguous', action='store_true', \
        default=False, help='Don\'t allow reads with ambiguous (N) bases')
    p.add_option('--trim-ambiguous', dest='trim_ambiguous', action='store_true', \
        default=False, help='Allow reads with ambiguous (N) bases, but trim' + \
                ' them out following alignment')
    p.add_option('--trim-running', dest='trimrunning', action='store_true', \
        default=False, help='Trim alignments using a running average')
    p.add_option('--remove-probe', dest='removeprobe', action='store_true', \
        default=False, help='Remove the probe sequence from the alignment.')
    p.add_option('--keep-probe', dest='keepprobe', action='store_true', \
        default=False, help='Add the probe sequence to the alignment.')
    p.add_option('--multiprocessing', dest='multiprocessing', action='store_true', \
        default=False, help='Use mutiple cores/processors.')
    p.add_option('--processors', dest='processors', action='store', \
        type='int', default=multiprocessing.cpu_count() - 2, \
        help='Number of cores to use.')

    p.add_option('--notstrict', dest='notstrict', action='store_true', \
        default=False, help='Allow alignments containing not all species.')

    p.add_option('--faircloth', dest='faircloth', action='store_true', \
        default=False, help='Take faircloth+stephens probe names')

    p.add_option('--verbose', dest='verbose', action='store_true', \
        default=False, help='List locus names while processing.')

    (options, arg) = p.parse_args()
    if not options.species:
        interfaceWarnings(p, "You must provide an expected species count per alignment")
    if options.processors < 1:
        options.processors = 1
    if not options.path:
        interfaceWarnings(p, None)
    if not os.path.isfile(options.path):
        interfaceWarnings(p, "You must provide a valid path to the input file or directory")
    if options.removeprobe and not options.probefile:
        interfaceWarnings(p, "You must provide a file of probe sequences if you are removing a probe")
    if os.path.isfile(options.output):
        interfaceWarnings("You must provide an output *directory*.")
    formats = {'clustal': '.clw',
                'emboss': '.emboss',
                'fasta': '.fa',
                'nexus': '.nex',
                'phylip': '.phylip',
                'stockholm': '.stockholm'}
    if options.format not in formats.keys():
        interfaceWarnings(p, "This is an unsupported output format")
    options.format_extension = formats[options.format]
    return options, arg


def outPFilename(path, fnn, extension):
    '''generalize the formatting of output filenames'''
    return (os.path.join(path, fnn)) + extension


def multiAlign(input, output):
    for locus_name, sequences, probe, options in iter(input.get, 'STOP'):
        results = singleAlign(locus_name, sequences, probe, options)
        output.put(results)
    return


def singleAlign(locus_name, sequences, probe, options):
    #pdb.set_trace()
    if options.verbose:
        print '\tLocus:  %s' % locus_name
    if options.trimflank or options.removeprobe:
        # give the probe a generic name so we can easily find it later
        probe.name = 'probe'
        probe.id = 'probe'
        probe.description = 'probe'
        # add the probe sequence to the locus
        sequences.append(probe)
        # format the locus reads to ('fasta')
        sequences.formatSequences()
        # create a tempfile to feed to muscle and feed it
        sequences.createTempFile()
        muscle = phyluce.muscle.Align(sequences.tempFile)
        muscle.run_alignment(consensus=False)
        muscle.get_probe_location()
        muscle.trim_alignment(method='trim', probe='remove')
    elif options.keepprobe:
        # give the probe a generic name so we can easily find it later
        probe.name = 'probe'
        probe.id = 'probe'
        probe.description = 'probe'
        # add the probe sequence to the locus
        sequences.append(probe)
        # format the locus reads to ('fasta')
        sequences.formatSequences()
        # create a tempfile to feed to muscle and feed it
        sequences.createTempFile()
        muscle = phyluce.muscle.Align(sequences.tempFile)
        muscle.run_alignment(consensus=False)
        muscle.trim_alignment(method='running-probe', window_size=20, threshold=0.5)
    elif options.trimrunning:
        #print sequences
        sequences.formatSequences()
        sequences.createTempFile()
        muscle = phyluce.muscle.Align(sequences.tempFile)
        muscle.run_alignment(consensus=False)
        muscle.trim_alignment(method='running', window_size=20, threshold=0.5)
    else:
        #print sequences
        sequences.formatSequences()
        sequences.createTempFile()
        muscle = phyluce.muscle.Align(sequences.tempFile)
        muscle.run_alignment(consensus=False)
        muscle.trim_alignment(method='notrim')
    if options.trim_ambiguous:
        muscle.trim_ambiguous_bases()
        return locus_name, muscle.perfect_trimmed_alignment
    else:
        return locus_name, muscle.trimmed_alignment


def q_runner(n_procs, loci, probes, options, function, *args):
    '''generic function used to start worker processes'''
    task_queue = multiprocessing.Queue()
    results_queue = multiprocessing.JoinableQueue()
    if args:
        arguments = (task_queue, results_queue,) + args
    else:
        arguments = (task_queue, results_queue,)
    results = []
    # reduce processer count if proc count > files
    if len(loci) < n_procs:
        n_procs = len(loci)
    for count, locus in enumerate(loci):
        if probes:
            task_queue.put([locus, loci[locus], probes[locus], options])
        else:
            task_queue.put([locus, loci[locus], None, options])
    for _ in range(n_procs):
        p = multiprocessing.Process(target=function, args=arguments).start()
        #print 'Starting %s' % function
    for _ in range(len(loci)):
        # indicated done results processing
        results.append(results_queue.get())
        results_queue.task_done()
    #tell child processes to stop
    for _ in range(n_procs):
        task_queue.put('STOP')
    # join the queue until we're finished processing results
    print 'Waiting for jobs to finish...'
    results_queue.join()
    # not closing the Queues caused me untold heartache and suffering
    task_queue.close()
    results_queue.close()
    return results


def makeOutPutDir(output):
    """create directory to store output"""
    if os.path.exists(output) == True:
        overwrite = raw_input('Path exists, overwrite? [Y/n]:')
        if overwrite.lower() == 'y' or 'yes':
            shutil.rmtree(output)
            os.mkdir(output)
        else:
            pass
    return output


def ambiguousBaseChecker(loci, locus, record):
    if not 'N' in record.seq:
        loci = buildLocusDict(loci, locus, record)
    else:
        print 'Skipping {0} because it contains ambiguous bases'.format(record.id)
    return loci


def buildLocusDict(loci, locus, record):
    if locus not in loci.keys():
        loci[locus] = Locus(record)
    else:
        loci[locus].append(record)
    return loci


def main():
    options, arg = interface()
    # create an instance of the ConcatenatedAlignment, if necc:
    if options.concatenate:
        concatenatedAlignment = ConcatenatedAlignment()
    else:
        alignmentBlocks = None
    if not options.output:
        options.output = options.path
    if options.trimflank or options.removeprobe or options.keepprobe:
        # create a dict of probe sequences
        #TODO: clean this up a bit
        probes = SeqIO.to_dict(SeqIO.parse(open(options.probefile, 'rU'), 'fasta'))
    else:
        probes = None
    # create a holder for all loci
    loci = {}
    # basically, cram all of the sequence records from a file into a
    # dictionary, indexed by the locus name.  This will allow us (soon)
    # to read from a single file of fastas, putting reads with their
    # correct loci then sending that group of locus records off for
    # processing.
    print 'Making output directory...'
    makeOutPutDir(options.output)
    print 'Building the locus dictionary...'
    if options.no_ambiguous:
        print 'Removing ALL sequences with ambiguous bases...'
    else:
        print 'NOT removing sequences with ambiguous bases...'
    for record in SeqIO.parse(open(options.path, 'rU'), 'fasta'):
        #pdb.set_trace()
        if not options.faircloth:
            locus = record.description.split('|')[1]
        else:
            locus = '_'.join([record.description.split('|')[0], \
                record.description.split('|')[1].split('_')[0]])
        #record.id = record.description.split('_')[0]           # added by NGC to fix seq IDs
        # skip records containing ambiguous bases
        if options.no_ambiguous:
            loci = ambiguousBaseChecker(loci, locus, record)
        else:
            loci = buildLocusDict(loci, locus, record)
    #pdb.set_trace()
    # iterate over loci to check for all species at a locus
    good_loci = {}
    for locus in loci:
        if options.notstrict:
            good_loci[locus] = loci[locus]
        elif len(loci[locus]) < options.species:
            #for critter in loci[locus]:
            #    #pdb.set_trace()
            #    print critter.id.split('_')[2]
            print 'Dropping Locus {0} because of missing species'.format(locus)
            #pdb.set_trace()
        else:
            good_loci[locus] = loci[locus]
    #pdb.set_trace()
    if options.multiprocessing:
        print 'Using %s cores to align DNA sequences...' % options.processors
        results = q_runner(options.processors, good_loci, probes, options, multiAlign)
    else:
        print 'Using 1 core to align DNA sequences...'
        results = []
        for locus in good_loci:
            if probes:
                probe = probes[locus]
            else:
                probe = None
            results.append(singleAlign(locus, loci[locus], probe, options))
    print 'Writing output files...'
    for locus, alignment in results:
        if alignment and alignment.get_alignment_length() > 50:
            # renamed from nexus to (generic) outp because function is generic now...
            print outPFilename(options.output, locus, options.format_extension)
            outp = open(outPFilename(options.output, locus, options.format_extension), 'w')
            #
            # NOTE:  this is a stop-gap measure to deal with different probe names
            #
            if options.format == "phylip":
                for r in alignment:
                    r.id = r.id.split('_')[-1]
            outp.write(alignment.format(options.format))
            outp.close()
        else:
            print '\tLocus %s not successfully aligned due to trimming errors - skipped from writing' % locus
            pass
    if options.concatenate:
        for locus, alignment in results:
            if alignment and alignment.get_alignment_length() > 50:
                concatenatedAlignment.concat(alignment)
            else:
                print '\tLocus %s not successfully aligned due to trimming errors - skipped from writing' % locus
                pass
        outp = open(outPFilename(options.output, 'concat', options.format_extension), 'w')
        try:
            outp.write(concatenatedAlignment.format(options.format))
        except:
            pdb.set_trace()
        outp.close()
        #pdb.set_trace()

if __name__ == '__main__':
    main()
