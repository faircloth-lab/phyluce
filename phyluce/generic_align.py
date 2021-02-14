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
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

import pdb


class GenericAlign(object):
    """docstring for Align"""

    def __init__(self, input):
        self.input = input
        self.alignment = None
        self.trimmed = None

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

    def _get_ends(self, seq):
        """find the start and end of sequence data for a given alignment row"""
        f = re.compile("^([-]+)")
        result = f.search(str(seq.seq))
        if result:
            start_gap = len(result.groups()[0])
        else:
            start_gap = 0
        r = re.compile("([-]+)$")
        result = r.search(str(seq.seq))
        if result:
            end_gap = len(result.groups()[0])
        else:
            end_gap = 0
        return start_gap, len(seq.seq) - end_gap

    def _gap_replacement(self, match, r="?"):
        """function called by replace_ends to add group of (missing data) characters"""
        if match.groups():
            return r * len(match.groups()[0])
        else:
            pass

    def _replace_ends(self, seq):
        """replace the ends of a given alignment with a character (usually gap/missing data)"""
        seq = re.sub("^([-]+)", self._gap_replacement, seq)
        seq = re.sub("([-]+)$", self._gap_replacement, seq)
        return seq

    def _alignment_consensus(self, alignment):
        """return consensus for an alignment object using BioPython"""
        consensus = []
        for pos in range(alignment.get_alignment_length()):
            col = alignment[:, pos].upper()
            cnt = Counter(col)
            # get item with max count
            base, occurrence = cnt.most_common(1)[0]
            # get proportion
            proportion = float(occurrence) / len(col)
            if proportion >= 0.5:
                consensus.append(base)
            else:
                consensus.append("N")
        return "".join(consensus)

    def _read(self, format):
        """read an alignment from the CLI - largely for testing purposes"""
        self.alignment = AlignIO.read(open(self.input, "rU"), format)

    def _record_formatter(self, trim, name):
        """return a string formatted as a biopython sequence record"""
        return SeqRecord(
            Seq(trim),
            id=name,
            name=name,
            description=name,
            annotations={"molecule_type": "DNA"},
        )

    def running_average(self, alignment, window_size, proportion, threshold):
        """
        trim an alignment (assuming default parameters) such that `proportion`
        of taxa have `threshold` of residues at each position across the first
        and last `window_size` column slice of the alignment
        """
        # iterate across the columns of the alignment and determine presence
        # or absence of base-identity in the column
        good_alignment = []
        # get count of taxa in alignment
        taxa = len(alignment)
        # get what constitutes the count of characters we need to
        # make a "majority", based on values of `proportion` passed
        # to function.
        majority_of_characters = int(round(proportion * taxa, 0))
        for column in range(alignment.get_alignment_length()):
            # get the count of different bases in a column
            column_count = Counter(alignment[:, column])
            # don't start considering base differences until we
            # have data from > required_characters (meaning we've
            # past the gappy parts of a given alignment)
            if column_count["-"] <= majority_of_characters:
                # we don't want to count gaps
                del column_count["-"]
                # alignment is "good" where the count of identities
                # at a given base is >= 50% across all taxa
                if column_count.most_common(1)[0][1] >= majority_of_characters:
                    good_alignment.append(True)
                # alignment is "bad" when we have < majority_of_characters
                # identities at a given base.
                else:
                    good_alignment.append(False)
            # alignment is also "bad" when we have > majority_of_characters
            # gaps in a column
            else:
                good_alignment.append(False)
        # convert good_alignment to array
        good_alignment = numpy.array(good_alignment)
        # replace numpy convolve with trimming in from left then in from right
        for start_clip in range(good_alignment.size):
            # make sure we start on a "good" base
            if good_alignment[start_clip] != False:
                # get successive window-sized slices
                window = good_alignment[
                    start_clip : (start_clip + window_size)
                ]
                proportion = float(sum(window)) / len(window)
                # stop if we hit a point where porportion of good bases > threshold
                # this is he "good" start of an alignment
                if proportion > threshold:
                    break
        reverse_good_alignment = good_alignment[::-1]
        for end_clip in range(reverse_good_alignment.size):
            # make sure we start on a "good" base
            if reverse_good_alignment[end_clip] != False:
                window = reverse_good_alignment[
                    end_clip : (end_clip + window_size)
                ]
                proportion = float(sum(window)) / len(window)
                # stop if we hit a point where porportion of good bases > threshold
                # this is he "good" start of an alignment
                if proportion >= threshold:
                    end_clip = reverse_good_alignment.size - end_clip
                    break
        return start_clip, end_clip

    def stage_one_trimming(
        self,
        alignment,
        window_size,
        proportion,
        threshold,
        min_len,
        replace_ends=False,
    ):
        """
        First stage (of 3) alignment trimming to find and trim edges of a given
        alignment.  Calls running_average function above to determine reasonable
        alignment start and end trimming for the entire alignment block.
        """
        # get the trim positions that we determine begin and end "good"
        # alignments
        start, end = self.running_average(
            alignment, window_size, proportion, threshold
        )
        # create a new alignment object to hold our alignment
        s1_trimmed = MultipleSeqAlignment([])
        for sequence in alignment:
            # ensure correct sequence alphabet or we'll get a conflict when
            # we try to generate a consensus
            # sequence.seq.alphabet = IUPAC.IUPACAmbiguousDNA()
            sequence.seq.annotations = {"molecule_type": "DNA"}
            if start >= 0 and end:
                trim = sequence[start:end]
                # ensure we don't just add a taxon with only gaps/missing
                # data and that alignments are >= min_len
                if (
                    set(trim) != set(["-"])
                    and set(trim) != (["?"])
                    and len(trim) >= min_len
                ):
                    if not replace_ends:
                        s1_trimmed.append(sequence[start:end])
                    else:
                        # replace end gaps with missing data character ?
                        # called on third iteration of trimming
                        repl = self._replace_ends(str(sequence[start:end].seq))
                        s1_trimmed.append(
                            self._record_formatter(repl, sequence.id)
                        )
                else:
                    s1_trimmed = None
                    break
            else:
                s1_trimmed = None
                break
        return s1_trimmed

    def stage_two_trimming(
        self, s1_trimmed, window_size, max_divergence, min_len
    ):
        """
        Alignment row-by-row trimming.  After stage one trimming, iterate
        over rows of alignment to find differences between the alignment
        consensus and the row (taxon) of data.  Trim those ends that differ
        from the consensus with > `divergence` across a `window_size` window.
        Goes to third round of filtering to remove edges that end up with only '----'
        characters to start or end alignment block.
        """
        # create new alignment object to hold trimmed alignment
        s2_trimmed = MultipleSeqAlignment([])
        # get consensus of alignment in array form
        consensus_array = numpy.array(
            list(self._alignment_consensus(s1_trimmed))
        )
        # iterate over each alignment sequence
        for sequence in s1_trimmed:
            sequence.seq.annotations = {"molecule_type": "DNA"}
            # ensure sequence is uppercase - consensus will be, too
            sequence = sequence.upper()
            # get the true ends of the sequence by walking in until we hit some bases
            start, end = self._get_ends(sequence)
            # convert sequence to array
            orig_seq_array = numpy.array(list(sequence))
            # trim down gaps at edges so they do not exert undue influence
            # on trimming the sequence row
            seq_array = orig_seq_array[start:end]
            # set default values for trim to `start` and `end`, just for safety
            # this ensure we don't carry anything over from previous iteration
            # (we shouldn't)
            bad_start = 0
            bad_end = len(sequence)
            # =============================================================
            # get first 5' => 3' positions that start a `window_size` block
            # of sequence having a divergence of less than `max_divergence`
            # from the consensus sequence of all alignments
            # =============================================================
            # compare the sequence to the consensus, returns an array of
            # boolean values representing equality relative to the consensus
            compare = seq_array != consensus_array[start:end]
            # begin working from 5' => 3' across `compare` array
            for bad_start in range(compare.size):
                # get successive window-sized slices
                window = compare[bad_start : bad_start + window_size]
                divergence = float(sum(window)) / window.size
                # stop if we hit a point where divergence < max_divergence
                if divergence < max_divergence:
                    break
            # reverse the `compare` array and begin working 3' => 5'
            reversed_compare = compare[::-1]
            for bad_end in range(reversed_compare.size):
                window = reversed_compare[bad_end : bad_end + window_size]
                divergence = float(sum(window)) / window.size
                # get 5 value slices
                if divergence < max_divergence:
                    bad_end = reversed_compare.size - bad_end
                    break
            # given original edge trimming and `bad_start`/`bad_end` values,
            # set the starting values of the sequece array to '-'
            orig_seq_array[: start + bad_start] = "-"
            orig_seq_array[start + bad_end :] = "-"
            trim = "".join(orig_seq_array)
            # ensure alignment consists of something other than '-' or '?'
            # and that alignments are >= min_len
            if (
                set(trim) != set(["-"])
                and set(trim) != (["?"])
                and len(trim) >= min_len
            ):
                s2_trimmed.append(self._record_formatter(trim, sequence.id))
            # if they're not, return None
            else:
                s2_trimmed = None
                break
        return s2_trimmed

    def trim_alignment(
        self,
        method="running",
        window_size=20,
        proportion=0.65,
        threshold=0.65,
        max_divergence=0.20,
        min_len=100,
    ):
        """
        Trim a given alignment from one of the alignment engines.  Uses three-pass
        approach - one to trim alignment block ends, one to trim row data, and a
        third to re-trim alignment block ends after trimming row data.  Drop alignments
        shorter than 100 bp.

        Returns self.trimmed
        """
        if method == "notrim":
            self.trimmed = self.alignment
        else:
            try:
                s1_trimmed = self.stage_one_trimming(
                    self.alignment, window_size, proportion, threshold, min_len
                )
                # pdb.set_trace()
                s2_trimmed = self.stage_two_trimming(
                    s1_trimmed, window_size, max_divergence, min_len
                )
                # cleanup any edges on which we've masked the data
                self.trimmed = self.stage_one_trimming(
                    s2_trimmed,
                    window_size,
                    proportion,
                    threshold,
                    min_len,
                    replace_ends=True,
                )
            # very generic exception statement here - self.alignment will have None
            # value, "X" will write to stdout, and we'll report failed loci at end
            # of run.
            except:
                pass


if __name__ == "__main__":
    aln = GenericAlign("../test-data/uce-3.nexus")
    aln._read("nexus")
    aln.trim_alignment()
    outf = open("alignment-test.nex", "w")
    outf.write(aln.trimmed.format("nexus"))
    outf.close()
    # pdb.set_trace()
