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

Given an input "monolithic" fasta file of UCE contigs (from phyluce_), break that file up into locus- or taxon-specific individual files.


phyluce_assembly_extract_contigs_to_barcodes
--------------------------------------------

Takes as input the LOG file created during ``phyluce_assembly_match_contigs_to_barcodes`` (below) and outputs a more nicely-formatted table of results.


phyluce_assembly_get_bed_from_lastz
-----------------------------------

Given a lastz_ file produced by phyluce_, convert those results to BED format.


phyluce_assembly_get_fasta_lengths
----------------------------------

Given an input FASTA-formatted file, summarize the info on contigs within that file and output summary statistics.


phyluce_assembly_get_fastas_from_match_counts
---------------------------------------------

Given a match-count file (produced from ``phyluce_assembly_get_match_counts``), output a monolithic FASTA-formatted file of UCE loci.


phyluce_assembly_get_fastq_lengths
----------------------------------

Given some input FASTQ data, output summary stastistics about those reads.


phyluce_assembly_get_match_counts
---------------------------------

Given results from ``phyluce_assembly_match_contigs_to_probes`` (below) and a configuration file, output those taxa and loci for which matches exist in the UCE database. The config file looks like:

.. code-block:: bash

    [all]
    alligator_mississippiensis
    anolis_carolinensis
    gallus_gallus
    mus_musculus


phyluce_assembly_match_contigs_to_barcodes
------------------------------------------

Given a directory of assembled contigs and a file containing an organismal barcode in FASTA format, check the contigs of all taxa in the directory for presence of the barcode sequence, extract that region of each contig for each taxon, and run the result against the BOLD database (for each taxon).  Useful for checking species ID and also searching for potential contamination.


phyluce_assembly_match_contigs_to_probes
----------------------------------------

Given a directory of assembled contigs and a file of UCE baits/probes, search the contigs for those that match part/all of a given bait/probe at some level of stringency.


phyluce_assembly_screen_probes_for_dupes
----------------------------------------

Check a probe/bait file for potential duplicate baits/probes.


Alignment
=========

phyluce_align_add_missing_data_designators
------------------------------------------

Sometimes alignments do not contain the same taxa as other alignments, and those "missing taxa" need to be added in.  This program allows you to add those missing entries for the missing taxa, although this is often not needed when using the phyluce_ concatenation tools (see below), which automatically deal with this problem.


phyluce_align_concatenate_alignments
------------------------------------

Given an input file of alignments, concatenate those alignments together and output either a ``--nexus`` or a ``--phylip`` formatted file, along with charset information (either as an extra file, for phylip, or within the ``--nexus`` formatted file).


phyluce_align_convert_degen_bases
---------------------------------

If there are IUPAC degenerate base codes within an alignment, convert those to "N".


phyluce_align_convert_one_align_to_another
------------------------------------------

Convert alignments between formats.  Can convert freely between FASTA, Nexus, Phylip, Phylip-relaxed, Clustal, Emboss, and Stockholm.


phyluce_align_explode_alignments
--------------------------------

Given an input directry of alignments, "explode" those files into taxon- or locus-specific sequence files.


phyluce_align_extract_taxa_from_alignments
------------------------------------------

Given a set of alignments and a list of taxa to keep or remove from the alignments, make a new directory of alignments with those taxa kept or removed.


phyluce_align_extract_taxon_fasta_from_alignments
-------------------------------------------------

Given a set of alignments and a taxon to extract from them, extract the data for the taxon and format those data as a FASTA file.


phyluce_align_filter_alignments
-------------------------------

Filter alignments having certain taxa or certain lengths and make a new directory without those alignments.


phyluce_align_format_concatenated_phylip_for_paml
-------------------------------------------------

This will convert a Phylip-formatted concatenated alignment for PAML's weird, internal format.  Not sure if PAML needs this format any longer.


phyluce_align_get_align_summary_data
------------------------------------

Given a directory of alignments, output summary statistics for those alignments quickly.


phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed
-----------------------------------------------------------

Given a directory of alignments, use gblocks_ to trim alignment edges and output a new folder of trimmed alignments.


phyluce_align_get_incomplete_matrix_estimates
---------------------------------------------

Given a directory of alignments, estimate the number of taxa present in various incomplete matrix scenarios.


phyluce_align_get_informative_sites
-----------------------------------

Given a directory of alignments, compute the number of informative sites.  Can output a list of these, which can be used for additional filtering by # of sites.


phyluce_align_get_only_loci_with_min_taxa
-----------------------------------------

Given a directory of alignments, filter those alignments for a minimum number of taxa, and output the filtered alignments to a new directory.


phyluce_align_get_ry_recoded_alignments
---------------------------------------

Given a directory of alignments, recode those as either "R/Y" or "0/1" and output the converted alignments to a new directory.


phyluce_align_get_smilogram_from_alignments
-------------------------------------------

Given a directory of alignments, ouput a CSV-formatted file showing the number of sites in the alignment from the center to the edges.  Can be input to R to make "smilogram" figures of UCE variation.


phyluce_align_get_taxon_locus_counts_in_alignments
--------------------------------------------------

Given a directory of alignments, get a count of taxa in each alignment.  Can be used to filter alignments based on occupancy.


phyluce_align_get_trimal_trimmed_alignments_from_untrimmed
----------------------------------------------------------

Given a directory of alignments, use trimAL_ to trim the alignments and write the trimmed alignments to a new directory.


phyluce_align_get_trimmed_alignments_from_untrimmed
---------------------------------------------------

Given a directory of alignments, use the phyluce_ edge-trimming algorithm to trim alignment edges and output a new folder of trimmed alignments.


phyluce_align_move_align_by_conf_file
-------------------------------------

Given an input configuration file, copy alignments present in the configuration file from one directory to a new directory.  Useful for filtering alignments. The format of the config file looks like:

.. code-block:: bash

    [all]
    locus-name-1.nexus
    locus-name-2.nexus
    locus-name-3.nexus

And this will copy all three loci to the specified (new) directory.


phyluce_align_randomly_sample_and_concatenate
---------------------------------------------

Given an input directory of alignments, randomly sample those alignments and output a file of those random sequences concatenated together.


phyluce_align_reduce_alignments_with_raxml
------------------------------------------

Given an input directory of alignments, use raxml_ to create "reduced" files of each alignment (where missing and low-info site patterns have been removed). Was needed to work around a pargenes_ bug but should not be necessary any longer.


phyluce_align_remove_empty_taxa
-------------------------------

Given an input directory of alignments, remove those taxa within each alignment having no data.


phyluce_align_remove_locus_name_from_files
------------------------------------------

Given an input directory of alignments that still contain the names of each UCE loci along with each taxon, strip the name of the UCE locus from each taxon and output the result into a new directory.


phyluce_align_screen_alignments_for_problems
--------------------------------------------

Given an input directory of alignments, screen those for problems such as weird nucleotide codes ("X") or runs of ambigious bases ("N").


phyluce_align_seqcap_align
--------------------------

Given a monolithic fasta file, align fasta sequences by locus and output the resulting alignments into a new directory.


phyluce_align_split_concat_nexus_to_loci
----------------------------------------

Given an input NEXUS-formatted file of a concatenated alignment (with charset info in the file), split the concatenated alignment back into component parts.


Genetrees
=========

phyluce_genetrees_generate_multilocus_bootstrap_count
-----------------------------------------------------

This is used with site and locus-based boostrap resampling.  Not particularly recommended any longer.


phyluce_genetrees_get_mean_bootrep_support
------------------------------------------

Given a set of input genetrees, compute the mean bootrep support among those trees.


phyluce_genetrees_get_tree_counts
---------------------------------

Given an input directory of alignments, uses the symmetric difference to count the number of similar and difference gene tree topologies.


phyluce_genetrees_rename_tree_leaves
------------------------------------

Given an input tree and a config file detailing the mapping of old names to new names, convert the old leaf names in a tree to the new names in the config file.  The config file format is similar to:

.. code-block:: bash

    [standard]
    acanthisitta_chloris_APP_002:Passeriformes__Acanthisittidae__acanthisitta_chloris_APP_002__4763__1051
    acrocephalus_arundinaceus_OUT_0054:Passeriformes__Acrocephalidae__acrocephalus_arundinaceus_OUT_0054__4729__1062
    aegithalos_caudatus_B10K_DU_002_10:Passeriformes__Aegithalidae__aegithalos_caudatus_B10K_DU_002_10__4572__1038
    aegotheles_bennettii_B10K_DU_029_76:Caprimulgiformes__Aegothelidae__aegotheles_bennettii_B10K_DU_029_76__4806__1058
    agapornis_roseicollis_OUT_0001:Psittaciformes__Psittaculidae__agapornis_roseicollis_OUT_0001__4760__1040
    agelaius_phoeniceus_OUT_0050:Passeriformes__Icteridae__agelaius_phoeniceus_OUT_0050__4786__1060

Where the old name is on the left of the colon and the new name is on the right of the colon.

phyluce_genetrees_sort_multilocus_bootstraps
--------------------------------------------

This is used with site and locus-based boostrap resampling.  Not particularly recommended any longer.

NCBI
====

phyluce_ncbi_chunk_fasta_for_ncbi
---------------------------------

Splits an input fasta file into chunks of 10,000 sequences because Sequin files should not contain more than 10,000 records.

phyluce_ncbi_prep_uce_align_files_for_ncbi
------------------------------------------

Given an input file of alignments, prep those for input to tbl2asn, which formats them for submission to NCBI.  Also takes a config file with a format similar to what follows.

.. code-block:: bash

    [exclude taxa]
    ichthyopis kohtaoensis
    plethodon chlorobryonis
    phalacrocorax carbo

    [exclude loci]
    uce-1809

    [metadata]
    molecule:DNA
    moltype:genomic
    location:genomic
    note:ultra conserved element locus {}
    specimen_voucher:{}

    [vouchers]
    Ardeotis kori:FLMNH 44254
    Balaeniceps rex:LSUMZ B13372
    Cathartes aura:LSUMZ B17242

    [remap]
    pterocles:pterocles exustus
    zanclostomus javanicus:sphyrapicus varius

Probes
======

phyluce_probe_easy_lastz
------------------------

Run an "easy" lastz_ search of one file against another file.


phyluce_probe_get_genome_sequences_from_bed
-------------------------------------------

Given an input BED file, extracts fasta information matching the coordinates in the BED file.


phyluce_probe_get_locus_bed_from_lastz_files
--------------------------------------------

Given a lastz_ results file where baits/probes were searched against a genome, output a BED-formatted file of the **locus** coodinates for each match of bait/baits to the genome.


phyluce_probe_get_multi_fasta_table
-----------------------------------

Make a table containing multi-way fasta information.


phyluce_probe_get_multi_merge_table
-----------------------------------

Make a table containing multi-way fasta information.


phyluce_probe_get_probe_bed_from_lastz_files
--------------------------------------------

Given a lastz_ results file where baits/probes were searched against a genome, output a BED-formatted file of the **bait** coodinates for each match of bait/baits to the genome.

phyluce_probe_get_screened_loci_by_proximity
--------------------------------------------

Given a FASTA file of properly formatted baits, keep only 1 locus (randomly) of those falling within a set distance from one another.


phyluce_probe_get_subsets_of_tiled_probes
-----------------------------------------

Given a bait file that contains baits designed from multiple organisms, prune that bait file to contain only those baits from a desired subset of organisms.


phyluce_probe_get_tiled_probe_from_multiple_inputs
--------------------------------------------------

Design baits from multiple input genomes.


phyluce_probe_get_tiled_probes
------------------------------

Design baits from a single input genome.


phyluce_probe_query_multi_fasta_table
-------------------------------------

Query a multifasta table.


phyluce_probe_query_multi_merge_table
-------------------------------------

Query a multimerge table.


phyluce_probe_reconstruct_uce_from_probe
----------------------------------------

From a UCE bait set, reconstruct the UCE locus used for design.


phyluce_probe_remove_duplicate_hits_from_probes_using_lastz
-----------------------------------------------------------

Given a bait set, and some lastz_ results of matching that bait set to itself, screen those probes from the bait set that match other probes (these are putative duplicates).


phyluce_probe_remove_overlapping_probes_given_config
----------------------------------------------------

Given a config file, filter baits from a set that are in the config file.


phyluce_probe_run_multiple_lastzs_sqlite
----------------------------------------

Use phyluce_ to run multiple lastz searches across multiple input genomes.


phyluce_probe_slice_sequence_from_genomes
-----------------------------------------

Given results from ``phyluce_probe_run_multiple_lastzs_sqlite``, slice fasta sequences from the genomes where there were matches.


phyluce_probe_strip_masked_loci_from_set
----------------------------------------

Remove baits from a putative bait set design where bait sequences have a high degree of masking.


Utilities
=========

phyluce_utilities_combine_reads
-------------------------------

Combine groups of reads based on an input file in config format.


phyluce_utilities_filter_bed_by_fasta
-------------------------------------

Filter a BED file of UCEs given a FASTA file of UCEs.


phyluce_utilities_get_bed_from_fasta
------------------------------------

Given an input fasta file of baits, prepared a BED-formatted file of their locations.


phyluce_utilities_merge_multiple_gzip_files
-------------------------------------------

Merge together multiple gzip files from the same sample.


phyluce_utilities_merge_next_seq_gzip_files
-------------------------------------------

Merge together multiple fastq gzip files from the next seq (these sometimes come as 4 files per sample).


phyluce_utilities_replace_many_links
------------------------------------

Use a config file to reformat many symlinks all at once.


phyluce_utilities_sample_reads_from_files
-----------------------------------------

Automatically randomly sample a fraction of reads from a fastq file.


phyluce_utilities_unmix_fasta_reads
-----------------------------------

Given an interleaved fastq file, convert that file into R1, R2, and singleton reads.


Workflow
========

phyluce_workflow
----------------

A single program to run a variety of Snakemake_ workflows.
