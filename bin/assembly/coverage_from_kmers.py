import pdb

import os
import sys
import math
import numpy
import argparse
from seqtools.sequence import fasta

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument('--input', nargs='?', default=sys.stdin)
    parser.add_argument('--output', nargs='?', default=sys.stdout)
    parser.add_argument('--csv', action='store_true', default = False)
    return parser.parse_args()

def get_average_read_length(input):
    try:
        stat_name = os.path.basename(input).split('.')[0].capitalize() + ".fastq.n-less.lengths"
        tld = os.path.split(os.path.dirname(input))[0]
        stat_path = os.path.join(tld, 'stats', stat_name)
        contents = open(stat_path, 'rU').read().split('\n')
        s_contents = contents[0].split(',')
        assert s_contents[0] == 'length'
        return float(s_contents[1])
    except:
        l = raw_input("What was n-less length? ")
        return float(l)

def main():
    args = get_args()
    avg_read_length = get_average_read_length(args.input)
    kmer = raw_input("What was kmer length? ")
    kmer = int(kmer)
    avg_c = []
    for read in fasta.FastaReader(args.input):
        s_read = read.identifier.split('_')
        ck = float(s_read[-1])
        c = ck * avg_read_length/(avg_read_length - kmer + 1)
        avg_c.append(c)
    avg_c = numpy.array(avg_c)
    if not args.csv:
        print "mean:\t", numpy.mean(avg_c)
        print "95ci:\t", 1.96 * (numpy.std(avg_c, ddof=1)/math.sqrt(len(avg_c)))
        print "min:\t", min(avg_c)
        print "max:\t", max(avg_c)
        print "median:\t", numpy.median(avg_c)
        print "<10x:\t", sum(avg_c < 10)
        print "<25x:\t", sum(avg_c < 25)
        print "<50x:\t", sum(avg_c < 50)
        print "<100x:\t", sum(avg_c < 100)
    else:
        print "{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(numpy.mean(avg_c), 
            1.96 * (numpy.std(avg_c, ddof=1)/math.sqrt(len(avg_c))), 
            min(avg_c), max(avg_c), numpy.median(avg_c), sum(avg_c < 10), 
            sum(avg_c < 25), sum(avg_c < 50), sum(avg_c < 100))
    
if __name__ == '__main__':
    main()
