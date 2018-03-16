.. include:: global.rst

.. _Tutorial II:

*********************************
Tutorial II: UCE Phasing UCE data
*********************************

The following workflow derives from Andermann et al. 2018 (https://doi.org/10.1101/255752) and focuses on phasing SNPs in UCE data.

To phase your UCE data, you need to have individual-specific "reference" contigs against which to align your raw reads.  Generally speaking, you can create these individual-specific reference contigs at several stages of the phyluce_ pipeline, and the stage at which you choose to do this may depend on the analyses that you are running.  Some of the factors that come into play when making this decision are the sequence divergence between your organisms of interest, the quality of your alignments, and the types of locus trimming you are employing.

The following give you two different options for creating these individual-specific FASTA reference sequences: (1) uses the UCE contigs identified just after running :ref:`UceExtraction` and (2) which assumes your organisms are closely-related, and that you have aligned your sequence data using mafft_ with edge-trimming (this means that you have reached the end of the :ref:`EdgeTrimming` section).

.. attention::  We have not yet fully-implemented code that you can
    use if you are trimming your alignment data with some other
    approach (e.g. gblocks_ or trimal_).

The first thing we will do is separate the alignments that we've generated into individual-specific FASTA files, then we will align raw reads to these FASTA files using bwa_, and we will call SNPs and phase the SNPs using samtools_.

Approach 1: Exploding the UCE FASTA library
-------------------------------------------

As outlined in :ref:`Tutorial I`, once you extract your contigs that are UCE loci, you can create individual-specific reference FASTA files by "exploding" the monolithic fasta file.  You can do that with the following:

.. code-block:: bash

    # explode the monolithic FASTA by taxon (you can also do by locus)
    phyluce_assembly_explode_get_fastas_file \
        --input all-taxa-incomplete.fasta \
        --output-dir exploded-fastas \
        --by-taxon

You may want to get stats on these exploded-fastas by running something like the following:

.. code-block:: bash

    # get summary stats on the FASTAS
    for i in exploded-fastas/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done

    # samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
    alligator-mississippiensis.unaligned.fasta,4315,3465679,803.170104287,3.80363492428,224,1794,823.0,980
    anolis-carolinensis.unaligned.fasta,703,400214,569.294452347,9.16433421241,224,1061,546.0,7
    gallus-gallus.unaligned.fasta,3923,3273674,834.482283966,4.26048496461,231,1864,852.0,1149
    mus-musculus.unaligned.fasta,825,594352,720.426666667,9.85933217965,225,1178,823.0,139


Approach 2: Exploding aligned and trimmed UCE sequences
-------------------------------------------------------

You can also choose loci that have already been aligned and trimmed as the basis for SNP calling and haplotype phasing.  The benefits of this approach is that the individual-specific reference contigs you are inputting to the process will be somewhat normalized across all of your individuals because you have already generated alignments from all of your UCE loci and trimmed the edges of these loci.

To follow this approach, first proceed through the :ref:`EdgeTrimming` section of :ref:`TutorialOne`.  Then, you can "explode" the directory of alignments you have generated to create separate FASTA files for each individual using the following (this assumes your alignments are in `mafft-nexus-edge-trimmed` as in the tutorial).

.. code-block:: bash

    # explode the alignment files in mafft-nexus-edge-trimmed by taxon create a taxon-specific FASTA
    phyluce_align_explode_alignments \
        --input mafft-nexus-edge-trimmed \
        --output-dir exploded-fastas \
        --by-taxon

You may want to get stats on these exploded-fastas by running something like the following:

.. code-block:: bash

    # get summary stats on the FASTAS
    for i in exploded-fastas/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done


Creating a configuration file
-----------------------------

Before you run the script, you have to create a configuration file, telling the program where the cleaned and trimmed fastq reads are stored for each sample and where to find the contig FASTA library for each sample.
The configuration file should look like in the following example and should be saved as e.g. ``phasing.conf``::

    [references]
    genus_species1:/path/to/uce/taxon-set1/exploded-fastas/genus_species1_contigs.fasta
    genus_species2:/path/to/uce/taxon-set1/exploded-fastas/genus_species2_contigs.fasta

    [individuals]
    genus_species1:/path/to/clean-fastq/genus_species1
    genus_species2:/path/to/clean-fastq/genus_species2

    [flowcell]
    genus_species1:XXYYZZ
    genus_species2:XXYYZZ



[references]
^^^^^^^^^^^^

In this section you simply state the sample ID (``genus_species1``) followed by a colon (``:``) and the full path to the sample-specific FASTA library which was generated in the previous step.

[individuals]
^^^^^^^^^^^^^

In this section you give the complete path to the cleaned and trimmed reads folder for each sample.

.. attention:: The cleaned reads used by this program should be generated by illumiprocessor_ because the folder structure of the cleaned reads files is assumed to be that of illumiprocessor_ . This means that the zipped fastq files (fastq.gz) have to be located in a subfolder with the name ``split-adapter-quality-trimmed`` within each sample-specific folder.

[flowcell]
^^^^^^^^^^

The flowcell section is meant to add flowcell information from the Illumina run to the header to the BAM file that is created.  This can be helpful for later identication of sample and run information.  If you do not know the flowcell information for the data you are processing, you can look inside of the `fastq.gz` file using a program like `less`.  The flowcell identifier is the set of digits and numbers after the 2nd colon.

.. code-block:: bash

  @J00138:149:HT23LBBXX:8:1101:5589:1015 1:N:0:ATAAGGCG+CATACCAC
              ^^^^^^^^^

Alternatively, you can enter any string of information here (no spaces) that you would like to help identify a given sample.


Mapping reads against contigs
-----------------------------

To map the fastq read files against the contig reference database for each sample, run:


.. code-block:: bash

    phyluce_snp_bwa_multiple_align \
        --config /path/to/phasing.conf \
        --output /path/to/mapping_results \
        --subfolder split-adapter-quality-trimmed


This will produce an output along these lines::

  2016-03-09 16:40:22,628 - phyluce_snp_bwa_multiple_align - INFO - ============ Starting phyluce_snp_bwa_multiple_align ============
  2016-03-09 16:40:22,628 - phyluce_snp_bwa_multiple_align - INFO - Version: 1.5.0
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --config: /path/to/phasing.conf
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --cores: 1
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --log_path: None
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --mem: False
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --no_remove_duplicates: False
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --output: /path/to/mapping_results
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --subfolder: split-adapter-quality-trimmed
  2016-03-09 16:40:22,629 - phyluce_snp_bwa_multiple_align - INFO - Argument --verbosity: INFO
  2016-03-09 16:40:22,630 - phyluce_snp_bwa_multiple_align - INFO - ============ Starting phyluce_snp_bwa_multiple_align ============
  2016-03-09 16:40:22,631 - phyluce_snp_bwa_multiple_align - INFO - Getting input filenames and creating output directories
  2016-03-09 16:40:22,633 - phyluce_snp_bwa_multiple_align - INFO - ---------------------- Processing genus_species1 ----------------------
  2016-03-09 16:40:22,633 - phyluce_snp_bwa_multiple_align - INFO - Finding fastq/fasta files
  2016-03-09 16:40:22,636 - phyluce_snp_bwa_multiple_align - INFO - File type is fastq
  2016-03-09 16:40:22,637 - phyluce_snp_bwa_multiple_align - INFO - Creating read index file for genus_species1-READ1.fastq.gz
  2016-03-09 16:40:33,999 - phyluce_snp_bwa_multiple_align - INFO - Creating read index file for genus_species1-READ2.fastq.gz
  2016-03-09 16:40:45,142 - phyluce_snp_bwa_multiple_align - INFO - Building BAM for genus_species1
  2016-03-09 16:41:33,195 - phyluce_snp_bwa_multiple_align - INFO - Cleaning BAM for genus_species1
  2016-03-09 16:42:03,410 - phyluce_snp_bwa_multiple_align - INFO - Adding RG header to BAM for genus_species1
  2016-03-09 16:42:49,518 - phyluce_snp_bwa_multiple_align - INFO - Marking read duplicates from BAM for genus_species1
  2016-03-09 16:43:26,917 - phyluce_snp_bwa_multiple_align - INFO - Creating read index file for genus_species1-READ-singleton.fastq.gz
  2016-03-09 16:43:27,066 - phyluce_snp_bwa_multiple_align - INFO - Building BAM for genus_species1
  2016-03-09 16:43:27,293 - phyluce_snp_bwa_multiple_align - INFO - Cleaning BAM for genus_species1
  2016-03-09 16:43:27,748 - phyluce_snp_bwa_multiple_align - INFO - Adding RG header to BAM for genus_species1
  2016-03-09 16:43:28,390 - phyluce_snp_bwa_multiple_align - INFO - Marking read duplicates from BAM for genus_species1
  2016-03-09 16:43:30,633 - phyluce_snp_bwa_multiple_align - INFO - Merging BAMs for genus_species1
  2016-03-09 16:44:05,811 - phyluce_snp_bwa_multiple_align - INFO - Indexing BAM for genus_species1
  2016-03-09 16:44:08,047 - phyluce_snp_bwa_multiple_align - INFO - ---------------------- Processing genus_species2 ----------------------
  ...


Phasing mapped reads
--------------------

In the previous step you mapped the reads against the contig FASTA file for each sample. The results are stored in the output folder in bam-format. Now you can start the actual phasing of the reads. This will sort the reads within each bam file into two separate bam files (``genus_species1.0.bam`` and ``genus_species1.1.bam``).
The program is very easy to run and just requires the path to the bam files (output folder from previous mapping program, ``/path/to/mapping_results``) and the path to the configuration file, which is the same file as used in the previous step (``/path/to/phasing.conf``):

.. code-block:: bash

    phyluce_snp_phase_uces \
        --config /path/to/phasing.conf \
        --bams /path/to/mapping_results/ \
        --output /path/to/phased_reads


The output is supposed to look like this::

  2016-03-09 17:31:43,790 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
  2016-03-09 17:31:43,790 - phyluce_snp_phase_uces - INFO - Version: 1.5.0
  2016-03-09 17:31:43,790 - phyluce_snp_phase_uces - INFO - Argument --bams: /path/to/mapping_results/
  2016-03-09 17:31:43,790 - phyluce_snp_phase_uces - INFO - Argument --config: /path/to/phasing.conf
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - Argument --conservative: False
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - Argument --cores: 1
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - Argument --log_path: None
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - Argument --output: /path/to/phased_reads
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - Argument --verbosity: INFO
  2016-03-09 17:31:43,791 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
  2016-03-09 17:31:43,793 - phyluce_snp_phase_uces - INFO - Getting input filenames and creating output directories
  2016-03-09 17:41:32,196 - phyluce_snp_phase_uces - INFO - ----------------------- Processing genus_species1 ----------------------
  2016-03-09 17:41:32,196 - phyluce_snp_phase_uces - INFO - Phasing BAM file for genus_species1
  2016-03-09 17:41:42,787 - phyluce_snp_phase_uces - INFO - Sorting BAM for genus_species1
  2016-03-09 17:41:44,239 - phyluce_snp_phase_uces - INFO - Sorting BAM for genus_species1
  2016-03-09 17:41:45,705 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 0
  2016-03-09 17:42:02,203 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 1
  2016-03-09 17:42:18,776 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file unphased
  2016-03-09 17:42:58,258 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 0 from FASTQ 0
  2016-03-09 17:42:58,273 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 1 from FASTQ 1
  2016-03-09 17:42:58,286 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file unphased from FASTQ unphased
  2016-03-09 17:42:58,298 - phyluce_snp_phase_uces - INFO - Checking for correct FASTA files
  2016-03-09 17:42:58,298 - phyluce_snp_phase_uces - INFO - Cleaning FASTA files
  2016-03-09 17:42:58,475 - phyluce_snp_phase_uces - INFO - Balancing FASTA files
  2016-03-09 17:42:58,627 - phyluce_snp_phase_uces - INFO - Symlinking FASTA files
  2016-03-09 17:42:58,627 - phyluce_snp_phase_uces - INFO - ---------------------- Processing genus_species2 ---------------------
  2016-03-09 17:42:58,628 - phyluce_snp_phase_uces - INFO - Phasing BAM file for genus_species2
  2016-03-09 17:43:02,459 - phyluce_snp_phase_uces - INFO - Sorting BAM for genus_species2
  2016-03-09 17:43:03,012 - phyluce_snp_phase_uces - INFO - Sorting BAM for genus_species2
  2016-03-09 17:43:03,565 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 0
  2016-03-09 17:43:11,131 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 1
  2016-03-09 17:43:18,723 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file unphased
  2016-03-09 17:43:37,441 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 0 from FASTQ 0
  2016-03-09 17:43:37,454 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 1 from FASTQ 1
  2016-03-09 17:43:37,464 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file unphased from FASTQ unphased
  2016-03-09 17:43:37,472 - phyluce_snp_phase_uces - INFO - Checking for correct FASTA files
  2016-03-09 17:43:37,473 - phyluce_snp_phase_uces - INFO - Cleaning FASTA files
  2016-03-09 17:43:37,633 - phyluce_snp_phase_uces - INFO - Balancing FASTA files
  2016-03-09 17:43:37,776 - phyluce_snp_phase_uces - INFO - Symlinking FASTA files
  2016-03-09 17:43:37,779 - phyluce_snp_phase_uces - INFO - ------------------ Merging alleles from all loci-----------------
  2016-03-09 17:43:38,577 - phyluce_snp_phase_uces - INFO - Wrote 819 loci for genus_species1
  2016-03-09 17:43:38,669 - phyluce_snp_phase_uces - INFO - Wrote 812 loci for genus_species2
  2016-03-09 17:43:38,675 - phyluce_snp_phase_uces - INFO - ================ Completed phyluce_snp_phase_uces ===============



The program automatically produces a consensus sequence for each of these phased bam files (= allele sequence) and stores these allele sequences of all samples in a joined FASTA file (``joined_allele_sequences_all_samples.fasta``). This allele FASTA is deposited in the subfolder ``fastas`` within your output folder (e.g. ``/path/to/phased_reads``) and can be used as input for the following alignment steps.
