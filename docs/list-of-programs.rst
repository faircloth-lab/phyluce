.. include:: global.rst

****************
List of Programs
****************

Below is a list, by function (directory), of the various programs included in
the * phyluce_ package.  This list includes ALL programs in phyluce_, and some may
be lightly tested.  This list of programs **does not** include programs
currently in the ``ATTIC`` directory of the repository, because these programs
are deprecated.

assembly
========

* phyluce_assembly_assemblo_abyss
* phyluce_assembly_assemblo_trinity
* phyluce_assembly_assemblo_velvet
* phyluce_assembly_copy_trinity_symlinks
* phyluce_assembly_explode_get_fastas_file
* phyluce_assembly_extract_contigs_to_barcodes
* phyluce_assembly_get_fasta_lengths
* phyluce_assembly_get_fastas_from_match_counts
* phyluce_assembly_get_fastq_lengths
* phyluce_assembly_get_match_counts
* phyluce_assembly_get_trinity_coverage
* phyluce_assembly_get_trinity_coverage_for_uce_loci
* phyluce_assembly_match_contigs_to_barcodes
* phyluce_assembly_match_contigs_to_probes
* phyluce_assembly_parse_trinity_coverage_for_uce_loci_log
* phyluce_assembly_parse_trinity_coverage_log
* phyluce_assembly_screen_probes_for_dupes


align
=====

* phyluce_align_add_missing_data_designators
* phyluce_align_convert_one_align_to_another
* phyluce_align_explode_alignments
* phyluce_align_extract_taxa_from_alignments
* phyluce_align_extract_taxon_fasta_from_alignments
* phyluce_align_filter_alignments
* phyluce_align_filter_characters_from_alignments
* phyluce_align_format_concatenated_phylip_for_paml
* phyluce_align_format_nexus_files_for_mrbayes
* phyluce_align_format_nexus_files_for_raxml
* phyluce_align_get_align_summary_data
* phyluce_align_get_bed_from_lastz
* phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed
* phyluce_align_get_incomplete_matrix_estimates
* phyluce_align_get_indels_from_alignments
* phyluce_align_get_informative_sites
* phyluce_align_get_only_loci_with_min_taxa
* phyluce_align_get_smilogram_from_alignments
* phyluce_align_get_taxon_locus_counts_in_alignments
* phyluce_align_get_trimmed_alignments_from_untrimmed
* phyluce_align_move_align_by_conf_file
* phyluce_align_output_list_of_taxon_counts
* phyluce_align_parallel_sate
* phyluce_align_randomly_sample_and_concatenate
* phyluce_align_remove_empty_taxa
* phyluce_align_remove_locus_name_from_nexus_lines
* phyluce_align_screen_alignments_for_problems
* phyluce_align_seqcap_align
* phyluce_align_split_concat_nexus_to_loci

genetrees
=========

* phyluce_genetrees_phybase
* phyluce_genetrees_reformat_raxml_output
* phyluce_genetrees_rename_tree_leaves
* phyluce_genetrees_run_raxml_bootstraps
* phyluce_genetrees_run_raxml_genetrees
* phyluce_genetrees_run_raxml_multilocus_bootstraps
* phyluce_genetrees_split_models_from_genetrees


ncbi
====

* phyluce_ncbi_chunk_fasta_for_ncbi
* phyluce_ncbi_example-prep.conf
* phyluce_ncbi_prep_uce_align_files_for_ncbi
* phyluce_ncbi_prep_uce_fasta_files_for_ncbi

snps
====

* phyluce_snp_bwa_align
* phyluce_snp_convert_vcf_to_snapp
* phyluce_snp_convert_vcf_to_structure
* phyluce_snp_find_snps_in_bed_interval
* phyluce_snp_get_dbsnp_freq_stats
* phyluce_snp_get_dbsnp_variability_for_all_uces
* phyluce_snp_prep_interval_list_file_for_picard
* phyluce_snp_screen_vcf_files
* phyluce_snp_summarize_vcf_file

utilities
=========

* phyluce_utilities_combine_reads
* phyluce_utilities_filter_bed_by_fasta
* phyluce_utilities_merge_multiple_gzip_files
* phyluce_utilities_merge_next_seq_gzip_files
* phyluce_utilities_replace_many_links
* phyluce_utilities_sample_reads_from_files
* phyluce_utilities_unmix_fasta_reads

probes
======

* phyluce_probe_easy_lastz
* phyluce_probe_get_clusters_from_bed
* phyluce_probe_get_clusters_from_taxa
* phyluce_probe_get_genome_sequences_from_bed
* phyluce_probe_get_locus_bed_from_lastz_files
* phyluce_probe_get_multi_fasta_table
* phyluce_probe_get_multi_merge_table
* phyluce_probe_get_probe_bed_from_lastz_files
* phyluce_probe_get_screened_loci_by_proximity
* phyluce_probe_get_subsets_of_tiled_probes
* phyluce_probe_get_tiled_probe_from_multiple_inputs
* phyluce_probe_get_tiled_probes
* phyluce_probe_query_multi_fasta_table
* phyluce_probe_query_multi_merge_table
* phyluce_probe_remove_duplicate_hits_from_probes_using_lastz
* phyluce_probe_remove_overlapping_probes_given_config
* phyluce_probe_run_multiple_lastzs_sqlite
* phyluce_probe_slice_sequence_from_genomes
* phyluce_probe_strip_masked_loci_from_set
