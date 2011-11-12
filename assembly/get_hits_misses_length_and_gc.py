"""
File: get_hits_misses_length_and_gc.py
Author: Brant Faircloth

Created by Brant Faircloth on 10 November 2011 14:11 PST (-0800)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

to plot:

    ggplot(d, aes(y = tm, x = length)) + \
        geom_point(aes(colour = present), alpha=0.5, size = 2) + \
        facet_wrap(~taxon) + scale_colour_hue(h=c(180, 360))

"""

import os
import sys
import sqlite3
from collections import defaultdict
from seqtools.sequence import fasta

import pdb

def main():
    uces = []
    # get all ids of probes in 2560 set
    for seq in fasta.FastaReader('../archive/probe-subset-2560-synthesized.fasta'):
        name_split = seq.identifier.split('_')
        if name_split[0] not in ['>chrE22C19W28','>chrUn']:
            iden = '_'.join(name_split[:2]).strip('>')
        else:
            iden = '_'.join(name_split[:3]).strip('>')
        uces.append(iden)
    # get names, lengths, and GC content of loci in dbase
    conn = sqlite3.connect('/Users/bcf/Git/brant/seqcap/Non-repo/probe.sqlite')
    cur = conn.cursor()
    metadata = defaultdict(dict)
    for uce in uces:
        cur.execute("SELECT cons, cons_len FROM cons WHERE seq = ?", (uce,))
        data = cur.fetchall()
        # ensure we only get one record back
        assert len(data) == 1, "More than one record"
        read, length = data[0]
        gc = round((read.count('C') + read.count('G')) / float(len(read)), 3)
        cur.execute('''SELECT count(*) FROM sureselect WHERE seq = ? AND
                selected = 1''', (uce,))
        data = cur.fetchall()
        assert len(data) == 1, "More than one record"
        count = data[0][0]
        if count > 1:
            cur.execute('''SELECT avg(tm), avg(masked_bases),
            avg(added_bases) from sureselect where seq = ? 
            group by seq''', (uce,))
        else:
            cur.execute('''SELECT tm, masked_bases,
            added_bases from sureselect where seq = ? 
            group by seq''', (uce,))
        tm, masked, added = cur.fetchall()[0]
        metadata[uce] = {
            'gc':gc,
            'length':length,
            'count':count,
            'tm':tm,
            'masked':masked,
            'added':added
            }
    cur.close()
    conn.close()
    conn = sqlite3.connect('../archive/birds-probe-matches.sqlite')
    cur = conn.cursor()
    taxa = [
            'anser_erythropus',
            'gallus_gallus',
            'pitta_guajana',
            'dromaius_novaehollandiae',
            'megalaima_virens',
            'struthio_camelus',
            'eudromia_elegans',
            'phalacrocorax_carbo',
            'urocolius_indicus',
        ]
    query = "SELECT {} FROM matches WHERE uce = ?".format(', '.join(taxa))
    for uce in metadata.keys():
        cur.execute(query, (uce.lower(),))
        data = cur.fetchall()
        for k,v in enumerate(data[0]):
            metadata[uce][taxa[k]] = v
        #pdb.set_trace()
    outfile = open('gc-length-species-matches.csv', 'w')
    outfile.write('uce,gc,length,count,tm,masked,added,present,taxon\n')
    for uce in sorted(metadata.keys()):
        for taxon in taxa:
            outfile.write('{},{},{},{},{},{},{},{},{}\n'.format(
                uce,
                metadata[uce]['gc'],
                metadata[uce]['length'],
                metadata[uce]['count'],
                metadata[uce]['tm'],
                metadata[uce]['masked'],
                metadata[uce]['added'],
                metadata[uce][taxon],
                taxon.replace('_',' ').capitalize()
            ))
    outfile.close()
        
if __name__ == '__main__':
    main()

