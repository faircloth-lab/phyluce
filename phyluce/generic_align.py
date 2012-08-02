#!/usr/bin/env python
# encoding: utf-8
"""
File: generic_align.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 March 2012 12:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import numpy
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
#from Bio.Align.Generic import Alignment
from Bio.Align import MultipleSeqAlignment

import pdb


class GenericAlign(object):
    """docstring for Align"""
    def __init__(self, input):
        self.input = input
        self.alignment = None
        self.trimmed_alignment = None
        self.perfect_trimmed_alignment = None

    def _clean(self, outtemp):
        if type(outtemp) is list:
            for f in outtemp:
                os.remove(f)
        else:
            os.remove(outtemp)
        # cleanup temp file
        try:
            os.remove(self.input)
        except:
            pass

    def _find_ends(self, forward=True):
        """determine the first (or last) position where all reads in an 
        alignment start/stop matching"""
        if forward:
            theRange = xrange(self.alignment.get_alignment_length())
        else:
            theRange = reversed(xrange(self.alignment.get_alignment_length()))
        for col in theRange:
            if '-' in self.alignment.get_column(col):
                pass
            else:
                break
        return col

    def _base_checker(self, bases, sequence, loc):
        """ensure that any trimming that occurs does not start beyong the
        end of the sequence being trimmed"""
        # deal with the case where we just want to measure out from the
        # middle of a particular sequence
        if len(loc) == 1:
            loc = (loc, loc)
        if not bases > len(sequence.seq[:loc[0]]) and \
            not bases > len(sequence.seq[loc[1]:]):
            return True

    def _record_formatter(self, temp):
        """return a string formatted as a biopython sequence record"""
        temp_record = SeqRecord(temp)
        return temp_record

    def _alignment_summary(self, alignment):
        """return summary data for an alignment object using the AlignInfo
        class from BioPython"""
        summary = AlignInfo.SummaryInfo(alignment)
        consensus = summary.dumb_consensus()
        return summary, consensus

    def _read(self, format):
        """read an alignment from the CLI - largely for testing purposes"""
        self.alignment = AlignIO.read(open(self.input, 'rU'), format)

    def get_probe_location(self):
        '''Pull the probe sequence from an alignment object and determine its position
        within the read'''
        # probe at bottom => reverse order
        for record in self.alignment[::-1]:
            if record.id == 'probe':
                start = re.search('^-*', str(record.seq))
                end = re.search('-*$', str(record.seq))
                # should be first record
                break
        # ooh, this seems so very backwards
        self.ploc = (start.end(), end.start(),)

    def running_average(self, window_size, threshold, proportion=0.3, k=None, running_probe=False):
        # iterate across the columns of the alignment and determine presence
        # or absence of base-identity in the column
        differences = []
        members = len(self.alignment)
        if not running_probe:
            for column in xrange(self.alignment.get_alignment_length()):
                column_values = self.alignment[:, column]
                # get the count of different bases in a column (converting
                # it to a set gets only the unique values)
                column_list = list(column_values)
                # use proportional removal of gaps
                if column_list.count('-') <= int(round(proportion * members, 0)):
                    column_list = [i for i in column_list if i != '-']
                #pdb.set_trace()
                if len(set(column_list)) > 1:
                    differences.append(0)
                else:
                    differences.append(1)
        else:
            for column in xrange(self.alignment.get_alignment_length()):
                column_values = list(self.alignment[:, column])
                # drop the index of the probe from the column_values
                del column_values[k]
                # get the count of different bases in a column (converting
                # it to a set gets only the unique values).
                #
                # no need to convert to a list here because it is already one
                if len(set(column_values)) > 1:
                    differences.append(0)
                else:
                    differences.append(1)
        differences = numpy.array(differences)
        weight = numpy.repeat(1.0, window_size) / window_size
        running_average = numpy.convolve(differences, weight)[window_size - 1:-(window_size - 1)]
        good = numpy.where(running_average >= threshold)[0]
        # remember to add window size onto end of trim
        try:
            start_clip, end_clip = good[0], good[-1] + window_size
        except IndexError:
            start_clip, end_clip = None, None
        return start_clip, end_clip, good

    def trim_alignment(self, method='edges', remove_probe=None, bases=None, consensus=True, window_size=20, threshold=0.5, proportion=0.3):
        """Trim the alignment"""
        if method == 'edges':
            # find edges of the alignment
            start = self._find_ends(forward=True)
            end = self._find_ends(forward=False)
        elif method == 'running':
            start, end, good = self.running_average(window_size, threshold, proportion=proportion)
        elif method == 'running-probe':
            # get position of probe
            for k, v in enumerate(self.alignment):
                if v.name == 'probe':
                    break
                else:
                    pass
            start, end = self.running_average(window_size, threshold, proportion, k, True)
        #pdb.set_trace()
        if method == 'notrim':
            self.trimmed_alignment = self.alignment
        else:
            # create a new alignment object to hold our alignment
            self.trimmed_alignment = MultipleSeqAlignment([], Gapped(IUPAC.ambiguous_dna, "-"))
            for sequence in self.alignment:
                # ignore the probe sequence we added
                if (method == 'edges' or method == 'running' or method == 'running-probe') and not remove_probe:
                    # it is totally retarded that biopython only gives us the option to
                    # pass the Alignment object a name and str(sequence).  Given this
                    # level of retardation, we'll fudge and use their private method
                    if start >= 0 and end:
                        sl = []
                        for k, v in enumerate(sequence.seq):
                            if k in good:
                                sl.append(v)
                            else:
                                sl.append('.')
                        pdb.set_trace()
                        # sequence to array
                        seq_array = numpy.array(list(str(sequence.seq)))
                        # reindex by good bases
                        seq_array = seq_array[good]
                        # convert to sequence object
                        new_seq = Seq(seq_array.tostring(), IUPAC.ambiguous_dna)
                        new_seq_record = SeqRecord(new_seq, id=sequence.id, name=sequence.name, description=sequence.description)
                        self.trimmed_alignment.append(new_seq_record)
                    else:
                        self.trimmed_alignment = None
                        break
                elif method == 'static' and not remove_probe and bases:
                    # get middle of alignment and trim out from that - there's a
                    # weakness here in that we are not actually locating the probe
                    # region, we're just locating the middle of the alignment
                    mid_point = len(sequence) / 2
                    if self._base_checker(bases, sequence, mid_point):
                        self.trimmed_alignment._records.append(
                                sequence[mid_point - bases:mid_point + bases]
                            )
                    else:
                        self.trimmed_alignment = None
                elif method == 'static' and not remove_probe and bases and self.ploc:
                    # get middle of alignment and trim out from that - there's a
                    # weakness here in that we are not actually locating the probe
                    # region, we're just locating the middle of the alignment
                    if self._base_checker(bases, sequence, self.ploc):
                        self.trimmed_alignment._records.append(
                                sequence[self.ploc[0] - bases:self.ploc[1] + bases]
                            )
                    else:
                        self.trimmed_alignment = None
                elif remove_probe and self.ploc:
                    # we have to drop to sequence level to add sequence slices
                    # where we basically slice around the probes location
                    temp = sequence.seq[:self.ploc[0]] + sequence.seq[self.ploc[1]:]
                    self.trimmed_alignment._records.append( \
                            self._record_formatter(temp)
                        )
                elif method == 'static' and remove_probe and bases and self.ploc:
                    if self._base_checker(bases, sequence, self.ploc):
                        temp = sequence.seq[self.ploc[0] - bases:self.ploc[0]] + \
                                sequence.seq[self.ploc[1]:self.ploc[1] + bases]
                        self.trimmed_alignment._records.append( \
                                self._record_formatter(temp)
                            )
                    else:
                        self.trimmed_alignment = None
        pdb.set_trace()
        # build a dumb consensus
        if consensus and self.trimmed_alignment:
            self.trimmed_alignment_summary, self.trimmed_alignment_consensus = \
                self._alignment_summary(self.trimmed_alignment)
        if not self.trimmed_alignment:
            print "\tAlignment {0} dropped due to trimming".format(self.alignment._records[0].description)

    def trim_ambiguous_bases(self):
        """snip ambiguous bases from a trimmed_alignment"""
        ambiguous_bases = []
        # do this by finding all ambiguous bases and then snipping the largest
        # chunk with no ambiguous bases from the entire alignment
        if not self.trimmed_alignment:
            self.perfect_trimmed_alignment = self.trimmed_alignment
        else:
            for column in xrange(0, self.trimmed_alignment.get_alignment_length()):
                if 'N' in self.trimmed_alignment[:,column]:
                    ambiguous_bases.append(column)
            maximum = 0
            maximum_pos = None
            #pdb.set_trace()
            if not ambiguous_bases:
                self.perfect_trimmed_alignment = self.trimmed_alignment
            if ambiguous_bases:
                # prepend and append the start and end of the sequence so consider
                # those chunks outside the stop and start of ambiguous base runs.
                ambiguous_bases.insert(0, 0)
                ambiguous_bases.append(self.trimmed_alignment.get_alignment_length() - 1)
                # create a new alignment object to hold our alignment
                self.perfect_trimmed_alignment = \
                    MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
                for pos in xrange(len(ambiguous_bases)):
                    if pos + 1 < len(ambiguous_bases):
                        difference = ambiguous_bases[pos + 1] - \
                            ambiguous_bases[pos]
                        if difference > maximum:
                            maximum = difference
                            maximum_pos = (pos, pos + 1)
                    else:
                        pass
                # make sure we catch cases where there is not best block
                if maximum_pos:
                    for sequence in self.trimmed_alignment:
                        self.perfect_trimmed_alignment.append(
                                sequence[ambiguous_bases[maximum_pos[0]] + 1:ambiguous_bases[maximum_pos[1]]]
                            )
                else:
                    self.perfect_trimmed_alignment = None

if __name__ == '__main__':
    pass
