.. include:: ../global.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

.. _List of Programs:

************************
List of Phyluce Programs
************************

Assembly
========

phyluce_assembly_assemblo_abyss
-------------------------------

Assemble fastq data for phyluce_ using abyss_.


phyluce_assembly_assemblo_spades
--------------------------------

Assemble fastq data for phyluce_ using spades_.


phyluce_assembly_assemblo_velvet
--------------------------------

Assemble fastq data for phyluce_ using velvet_.


phyluce_assembly_explode_get_fastas_file
----------------------------------------

Given an input "monolithic" fasta file of UCE contigs (from phyluce_), break that file up into locus- or taxon-specific, individual files.


phyluce_assembly_extract_contigs_to_barcodes
--------------------------------------------

Takes as input the LOG file created during ``phyluce_assembly_match_contigs_to_barcodes`` (below), and outputs a more nicely-formatted table of results.


phyluce_assembly_get_bed_from_lastz
-----------------------------------

Given a lastz_ file produced by phyluce_, convert those results to BED format.


phyluce_assembly_get_fasta_lengths
----------------------------------

Given an input FASTA-formatted file, summarize the info on contigs within that file and output summary statistics.


phyluce_assembly_get_fastas_from_match_counts
---------------------------------------------

Given a match-count file (produced from ``phyluce_assembly_get_match_counts``, output a monolithic FASTA-formatted file of UCE loci.


phyluce_assembly_get_fastq_lengths
----------------------------------

Given some input FASTQ data, output summary stastistics about those reads.


phyluce_assembly_get_match_counts
---------------------------------

Given results from ``phyluce_assembly_match_contigs_to_probes`` (below) and a configuration file, output those taxa and loci for which matches exist in the UCE database. 


phyluce_assembly_match_contigs_to_barcodes
------------------------------------------

Given a directory of assembled contigs and a file containing an organismal barcode, check the contigs for presence of the barcode sequence, extract that region of the contig, and run the result against the BOLD database.  Useful for checking species ID and also searching for potential contamination.


phyluce_assembly_match_contigs_to_probes
----------------------------------------

Given a directory of assembled contigs and a file of UCE baits/probes, search the contigs for those that match part/all of a given bait/probe at some level of stringency.


phyluce_assembly_screen_probes_for_dupes
----------------------------------------

Check a probe/bait file for potential duplicate baits/probes.


Alignment
=========





Genetree
NCBI
Probes
Utilities
Workflow
