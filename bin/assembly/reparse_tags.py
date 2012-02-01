import pdb
import sys
import argparse
from Levenshtein import hamming
from seqtools.sequence import fastq
from seqtools.sequence import transform
from collections import Counter

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastqs from input based on tag sequence')
    parser.add_argument('--input', nargs='?', default=sys.stdin)
    parser.add_argument('--tag', required = True, dest = 'tag')
    parser.add_argument('--output', nargs='?', default=sys.stdout)
    return parser.parse_args()

def main():
    args = get_args()
    #pdb.set_trace()
    fastqs = fastq.FastqReader(args.input)
    out_fastq = fastq.FastqWriter(args.output)
    for read in fastqs:
        read_tag = read.identifier.split('#')[-1].split('/')[0]
        if hamming(read_tag, args.tag) <= 1:
            out_fastq.write(read)
        else:
            pass
    out_fastq.close()

if __name__ == '__main__':
    main()
