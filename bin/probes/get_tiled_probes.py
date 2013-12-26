#!/usr/bin/env python
# encoding: utf-8

"""
py_tiler.py

Created by Brant Faircloth on 2010-05-27.
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

Updated:  2011-07-15

This program "designs" selects probes based on a tile length, tile density
and tiling scheme (strict or flexibly overlap).

USAGE:  python py_tiler.py --probe-length=120 --tiling-density=1.2 \
            --overlap="flush-left" --input=chicken-mhc-cd1c-to-bg.fa \
            --masking=0.25 --output=test.fa --remove-ambiguous \
            --remove-gc --coords="chr16:1-250,600" \
            --bed=chicken-mhc.bed

Masking option gives the frequency of masked bases to accept in a probe.
Remove-gc removes probes 30 < x < 70 GC %.  Remove-ambiguous removes any
probes having "N" bases - since these are not easy to synthesize and could
cause problems.

"""

import os
import re
import sys
import argparse
import tempfile

from Bio import SeqIO
from phyluce import lastz
from phyluce.helpers import get_dupes

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Tile sequence capture probes across fastas.')
    parser.add_argument(
        '--input',
        required=True,
        dest='input',
        type=str,
        help='The path to the input file',
    )
    parser.add_argument(
        '--output',
        required=True,
        type=str,
        default=None,
        help='The path to the output file'
    )
    parser.add_argument(
        '--probe-prefix',
        required=True,
        type=str,
        default=None,
        help='The prefix (e.g. "uce-") to add to all probes designed'
    )
    parser.add_argument(
        '--designer',
        required=True,
        type=str,
        default=None,
        help='Your last name (to indicate who designed the probes)'
    )
    parser.add_argument(
        '--probe-length',
        dest='length',
        type=int,
        default=120,
        help='The length of the probes sequence to design'
    )
    parser.add_argument(
        '--tiling-density',
        dest='density',
        type=float,
        default=2,
        help='The tiling density'
    )
    parser.add_argument(
        '--overlap',
        type=str,
        choices=['middle', "flush-left"],
        default='middle',
        help='The method of tiling'
    )
    parser.add_argument(
        '--probe-bed',
        type=str,
        default=None,
        help='The path to an output file for outputting the probe coordinates in BED format'
    )
    parser.add_argument(
        '--locus-bed',
        type=str,
        default=None,
        help='The path to an output file for outputting the locus coordinates in BED format'
    )
    parser.add_argument(
        '--masking',
        dest='mask',
        type=float,
        default=None,
        help='The maximum frequency of per-probe masking allowed containing the sequence'
    )
    parser.add_argument(
        '--do-not-remove-ambiguous',
        dest='amb',
        action='store_true',
        default=True,
        help='Do not remove loci with probes containing ambiguous bases'
    )
    parser.add_argument(
        '--remove-gc',
        dest='gc',
        action='store_true',
        default=False,
        help='Remove loci with GC content outside 30 <= GC <= 70'
    )
    return parser.parse_args()


def middle_overlapper(region, args):
    '''
    for the middle class, the idea is that you start designing probes at
    the middle of the sequence and work outwards towards the ends.  At the
    middle, the typical overlap is split to each side of the middle position,
    like so (space in sequence denotes middle)

    probe-top-2     TTGATCAGCGGCCC
    probe-top-1              CGGCCCCTTCCGA G ATTA
    sequence    TGGATTGATCAGCGGCCCCTTCCCGA G ATTAAACTTGTAGCAGCTGATACACTTGGC
    probe-bottom-1                    CCGA G ATTAAACTTGTAG
    probe-bottom-2                                ACTTGTAGCAGCTGAT

    '''
    seq_len = len(region.seq)
    if seq_len == 179:
        extra = 1
    else:
        extra = 0
    # determine the degree of overlap between tiles
    tile_overlap = args.length - (args.length/args.density)
    tile_non_overlap = args.length - tile_overlap
    coords = []
    #pdb.set_trace()
    middle = seq_len/2
    halfsies = tile_overlap/2
    r_prb_strt = middle - halfsies - extra
    l_prb_strt = middle + halfsies + 2*extra
    end = 0
    while r_prb_strt + args.length <= seq_len:
        end = r_prb_strt + args.length
        coords.append((int(r_prb_strt), int(end)))
        r_prb_strt += tile_non_overlap
    start = l_prb_strt
    while l_prb_strt - args.length >= 0:
        start = l_prb_strt - args.length
        coords.append((int(start), int(l_prb_strt)))
        l_prb_strt -= tile_non_overlap
    return coords


def left_flush_overlapper(region, args):
    """
    for the left-flush, the idea is that you start designing probes at
    the 5' end of the sequence you are targetting and work outwards towards the
    3'.  The resulting probes will the flush on the left and, depending on
    the sequence, ragged on the right.

    probe-top-2     TTGATCAGCGGCCC
    probe-top-1              CGGCCCCTTCCGA G ATTA
    sequence        TTGATCAGCGGCCCCTTCCCGA G ATTAAACTTGTAGCAGCTGATACACTTGGC
    probe-bottom-1                    CCGA G ATTAAACTTGTAG
    probe-bottom-2                                ACTTGTAGCAGCTGAT
    """
    # determine the degree of overlap between tiles
    #pdb.set_trace()
    tile_overlap = args.length - (args.length/args.density)
    if tile_overlap == 0:
        step = 0
    else:
        step = int(round(tile_overlap))
    starts = range(0, len(region.seq), args.length - step)
    coords = [(start, start + args.length) for start in starts]
    return coords


def dots(letter=None):
    """flush, to stdout, some indicator of progress when screening"""
    if not letter:
        sys.stdout.write(".")
    else:
        sys.stdout.write("{}".format(letter))
    sys.stdout.flush()


def check_for_dupes(probe_set):
    """create some temp files and search a newly-designed probe-set for dupes"""
    # write to a tempfile
    f = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for ps in probe_set:
        SeqIO.write(ps, f, 'fasta')
    f.close()
    # align f to itself
    lz = lastz.Align(f.name, f.name, 70, 80)
    lz.run()
    dupes = get_dupes(lz.output, pos=2)
    os.remove(f.name)
    os.remove(lz.output)
    return dupes


def prune_probe_set(probe_set, dupes):
    """prune a probe set of duplicates"""
    final = []
    for ps in probe_set:
        good = []
        for probe in ps:
            name = '_'.join(probe.id.split("|")[:2])
            if name not in dupes:
                good.append(probe)
        final.append(good)
        print "Removed {0} probes that were likely duplicates".format(len(ps) - len(good))
    return final


def meta_to_dict(s):
    return dict([i.split(':') for i in s.split(',')])


def main():
    """
    main loop
    TODO: clean
    """
    args = get_args()
    probe_set = []
    print "Probes removed for masking (.) / low GC % (G) / ambiguous bases (N):"
    for locus_count, locus in enumerate(SeqIO.parse(open(args.input, 'rU'), 'fasta')):
        #pdb.set_trace()
        global_coords = locus.description.split('|')[1]
        global_chromo, global_chromo_positions = global_coords.split(':')
        global_chromo_start, global_chromo_end = [int(i) for i in global_chromo_positions.split('-')]
        if args.overlap == 'middle':
            coords = middle_overlapper(locus, args)
        elif args.overlap == 'flush-left':
            coords = left_flush_overlapper(locus, args)
        coords.sort()
        probes = []
        for k, coord in enumerate(coords):
            # need to get global starts and ends for probe
            global_probe_start = global_chromo_start + coord[0]
            global_probe_end = global_chromo_start + coord[1]
            #
            probe_id = '{0}{1}_p{2}'.format(
                args.probe_prefix,
                locus_count,
                k + 1
            )
            probe_description = " |designer:{0},probes-locus:{1},probes-probe:{2},probes-global-chromo:{3},probes-global-start:{4},probes-global-end:{5},probes-local-start:{6},probes-local-end:{7}".format(
                args.designer,
                locus_count,
                k + 1,
                global_chromo,
                global_probe_start,
                global_probe_end,
                coord[0],
                coord[1]
            )
            # slice the sequence
            probe = locus[coord[0]:coord[1]]
            # set slice to uppercase
            probe.seq = probe.seq.upper()
            # set the id and description
            probe.id, probe.name = probe_id, probe_id
            probe.description = probe_description
            masked = sum([1 for base in probe.seq if base.islower()]) / len(probe.seq)
            gc = sum([1. for base in probe.seq if base.upper() == 'C' or base.upper() == 'G']) / len(probe.seq)
            if args.mask and masked >= args.mask:
                dots('M')
            elif args.amb and ('N' in probe.seq or 'n' in probe.seq):
                dots('N')
            elif args.gc and (gc > 0.7 or gc < 0.3):
                dots('G')
            elif len(probe.seq) < args.length:
                dots('L')
            else:
                probes.append(probe)
        if probes:
            probe_set.append(probes)
    cons_count = len(probe_set)
    probe_count = 0
    probe_count = sum([len(cons) for cons in probe_set])
    print '\n\n'
    print 'Conserved locus count = {0}'.format(cons_count)
    print 'Probe Count = {0}'.format(probe_count)
    outp = open(args.output, 'w')
    if args.probe_bed:
        outpb = open(args.probe_bed, 'w')
        outpb.write('track name=get_tiled_probes description="get_tiled_probes designed probes" useScore=1 useScore=1 itemRgb="On"\n')
    if args.locus_bed:
        outlb = open(args.locus_bed, 'w')
        outlb.write('track name=get_tiled_probes_loci description="get_tiled_probes loci" useScore=1 useScore=1 itemRgb="On"\n')
    for ps in probe_set:
        lb_coords = []
        for probe in ps:
            meta_dict = meta_to_dict(probe.description)
            lb_coords.extend([int(meta_dict['probes-global-start']), int(meta_dict['probes-global-end'])])
            outp.write(">{0}{1}\n{2}\n".format(
                probe.id,
                probe.description,
                probe.seq
            ))
            if args.probe_bed:
                outpb.write("{}\t{}\t{}\t{}\t450\t+\t0\t0\t0,0,205\n".format(
                    meta_dict['probes-global-chromo'],
                    meta_dict['probes-global-start'],
                    meta_dict['probes-global-end'],
                    probe.id
                ))
        lb_coords.sort()
        mn = min(lb_coords)
        mx = max(lb_coords)
        try:
            assert int(mx) > int(mn)
        except:
            pdb.set_trace()
        if args.locus_bed:
            outlb.write("{0}\t{1}\t{2}\t{3}{4}\t450\t+\t0\t0\t0,0,205\n".format(
                meta_dict['probes-global-chromo'],
                mn,
                mx,
                args.probe_prefix,
                meta_dict['probes-locus']
            ))

    if args.probe_bed:
        outpb.close()
    if args.locus_bed:
        outlb.close()
    outp.close()

if __name__ == '__main__':
    main()
