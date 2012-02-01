import os
import sys
import argparse
from seqcap.lib import lastz
from operator import itemgetter
from collections import defaultdict

def get_name(header, splitchar = "_", items = 2):
    """use own function vs. import from match_contigs_to_probes - we don't want lowercase"""
    if splitchar:
        return "_".join(header.split(splitchar)[:items]).lstrip('>')
    else:
        return header.lstrip('>')

def get_dupes(lastz_file):
    matches = defaultdict(list)
    dupes = set()
    for lz in lastz.Reader(lastz_file):
        target_name = get_name(lz.name1, "|", 1)
        query_name = get_name(lz.name2, "|", 1)
        matches[target_name].append(query_name)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k,v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    return dupes

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_names_from_config(config, group):
    try:
        return [i[0] for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None