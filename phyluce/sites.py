# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 22 May 2015 15:32 CDT (-0500)
"""

from collections import Counter


def get_informative_sites(count):
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    sufficient_sites = len(count)
    if sufficient_sites >= 2:
        sufficient_sequences = sum([1 for i in count.values() if i >= 2])
        if sufficient_sequences >= 2:
            return True
    return False


def get_differences(count):
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    # remove X
    del count['X']
    sufficient_sites = len(count)
    # counted, different = (1,1)
    if sufficient_sites >= 2:
        return (1, 1)
    # counted, not different = (1,0)
    elif sufficient_sites >= 1 and count.most_common()[0][1] > 1:
        return (1, 0)
    # not counted, not different = (0,0)
    else:
        return (0, 0)


def compute_informative_sites(align):
    informative_sites = []
    differences = []
    counted_sites =  []
    for idx in xrange(align.get_alignment_length()):
        col = align[:, idx].upper()
        count = Counter(col)
        if get_informative_sites(count):
            informative_sites.append(1)
        else:
            informative_sites.append(0)
        diff = get_differences(count)
        if diff == (1, 1):
            counted_sites.append(1)
            differences.append(1)
        elif diff == (1, 0):
            differences.append(0)
            counted_sites.append(1)
        else:
            differences.append(0)
            counted_sites.append(0)
    return sum(informative_sites), sum(differences), sum(counted_sites)
