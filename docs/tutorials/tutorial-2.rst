.. include:: ../global.rst

.. _TutorialTwo:

*********************************
Tutorial II: Phasing UCE data
*********************************

The following workflow derives from Andermann et al. 2018 (https://doi.org/10.1093/sysbio/syy039) and focuses on phasing SNPs in UCE data.

To phase your UCE data, you need to have individual-specific "reference" contigs against which to align your raw reads.  Generally speaking, you can create these individual-specific reference contigs at several stages of the phyluce_ pipeline, and the stage at which you choose to do this may depend on the analyses that you are running.  That said, I think that the best way to proceed uses edge-trimmed exploded alignments as your reference contigs, aligns raw reads to those, and uses the exploded alignments and raw reads to phase your data.

.. attention::  We have not implemented code that you can
    use if you are trimming your alignment data with some other
    approach (e.g. gblocks_ or trimal_).

Exploding aligned and trimmed UCE sequences
-------------------------------------------

Probably the best way to proceed (you can come up with other ways to do this) is to choose loci that have already been aligned and edge-trimmed as the basis for SNP calling and haplotype phasing.  The benefit of this approach is that the individual-specific reference contigs you are inputting to the process will be somewhat normalized across all of your individuals because you have already generated alignments from all of your UCE loci and trimmed the edges of these loci.

To follow this approach, first proceed through the :ref:`EdgeTrimming` section of :ref:`TutorialOne`.  Then, you can "explode" the directory of alignments you have generated to create separate FASTA files for each individual using the following (this assumes your alignments are in `mafft-nexus-edge-trimmed` as in the tutorial).

.. code-block:: bash

    # explode the alignment files in mafft-nexus-edge-trimmed by taxon create a taxon-specific FASTA
    phyluce_align_explode_alignments \
        --alignments mafft-nexus-edge-trimmed \
        --input-format nexus \
        --output mafft-nexus-edge-trimmed-exploded \
        --by-taxon

The current directory structure should look like (I've collapsed a number of
branches in the tree):

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ...
    ├── taxon-sets
    │       └── all
    │           ├── all-taxa-incomplete.conf
    │           ├── all-taxa-incomplete.fasta
    │           ├── all-taxa-incomplete.incomplete
    │           ├── exploded-fastas
    │           ├── log
    │           ├── mafft-nexus-edge-trimmed
    |           └── mafft-nexus-edge-trimmed-exploded
    |               ├── alligator_mississippiensis.fasta
    |               ├── anolis_carolinensis.fasta
    |               ├── gallus_gallus.fasta
    |               └── mus_musculus.fasta
    ...
    └── uce-search-results

You may want to get stats on these exploded-fastas by running something like the following:

.. code-block:: bash

    # get summary stats on the FASTAS
    for i in mafft-nexus-edge-trimmed-exploded/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done


Creating a re-alignment configuration file
------------------------------------------

Before aligning raw reads back to these reference contigs using ``bwa``, you have to create a configuration file, which tells the program where the cleaned and trimmed fastq reads are stored for each sample and where to find the reference FASTA file for each sample. The configuration file should look like in the following example and should be saved as e.g. ``phasing.conf``

.. code-block:: text

    [references]
    alligator_mississippiensis:/Users/bcf/tmp/phyluce/mafft-nexus-edge-trimmed-exploded/alligator_mississippiensis.fasta
    anolis_carolinensis:/Users/bcf/tmp/phyluce/mafft-nexus-edge-trimmed-exploded/anolis_carolinensis.fasta
    gallus_gallus:/Users/bcf/tmp/phyluce/mafft-nexus-edge-trimmed-exploded/gallus_gallus.fasta
    mus_musculus:/Users/bcf/tmp/phyluce/mafft-nexus-edge-trimmed-exploded/mus_musculus.fasta

    [individuals]
    alligator_mississippiensis:/Users/bcf/tmp/phyluce/clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
    anolis_carolinensis:/Users/bcf/tmp/phyluce/clean-fastq/anolis_carolinensis/split-adapter-quality-trimmed/
    gallus_gallus:/Users/bcf/tmp/phyluce/clean-fastq/gallus_gallus/split-adapter-quality-trimmed/
    mus_musculus:/Users/bcf/tmp/phyluce/clean-fastq/mus_musculus/split-adapter-quality-trimmed/

    [flowcell]
    alligator_mississippiensis:D1HTMACXX
    anolis_carolinensis:C0DBPACXX
    gallus_gallus:A8E3E
    mus_musculus:C0DBPACXX


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

Alternatively, you can enter any string of information here (no spaces) that you would like to help identify a given sample (e.g. ``XXYYZZ``).


Mapping reads against contigs
-----------------------------

To map the fastq read files against the contig reference database for each sample, run the folliwing.  This will use ``bwa mem`` to map the raw reads to the "reference" contigs:

.. code-block:: bash

    phyluce_snp_bwa_multiple_align \
        --config phasing.conf \
        --output multialign-bams \
        --cores 12 \
        --log-path log \
        --mem

This will produce an output along these lines

.. code-block:: text

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

In the previous step you aligned your sequence reads against the reference FASTA file for each sample. The results are stored in the output folder in ``bam`` format. Now you can start the actual phasing of the reads. This will analyze and sort the reads within each bam file into two separate bam files (``genus_species1.0.bam`` and ``genus_species1.1.bam``).

The program is very easy to run and just requires the path to the bam files (output folder from previous mapping program, ``/path/to/mapping_results``) and the path to the configuration file, which is the same file as used in the previous step (``/path/to/phasing.conf``).  Then, run:

.. code-block:: bash

    phyluce_snp_phase_uces \
        --config phasing.conf \
        --bams multialign-bams \
        --output multialign-bams-phased-reads

The output will look something like the following

.. code-block:: text

    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Version: git 0babc1a
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Argument --bams: /Users/bcf/tmp/phyluce/multialign-bams
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Argument --config: /Users/bcf/tmp/phyluce/phasing.conf
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Argument --conservative: False
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Argument --cores: 1
    2018-07-27 09:58:35,963 - phyluce_snp_phase_uces - INFO - Argument --log_path: None
    2018-07-27 09:58:35,964 - phyluce_snp_phase_uces - INFO - Argument --output: /Users/bcf/tmp/phyluce/multialign-bams-phased-reads
    2018-07-27 09:58:35,964 - phyluce_snp_phase_uces - INFO - Argument --verbosity: INFO
    2018-07-27 09:58:35,964 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
    2018-07-27 09:58:35,964 - phyluce_snp_phase_uces - INFO - Getting input filenames and creating output directories
    2018-07-27 10:02:10,526 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Version: git 0babc1a
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --bams: /Users/bcf/tmp/phyluce/multialign-bams
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --config: /Users/bcf/tmp/phyluce/phasing.conf
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --conservative: False
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --cores: 1
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --log_path: None
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --output: /Users/bcf/tmp/phyluce/multialign-bams-phased-reads
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - Argument --verbosity: INFO
    2018-07-27 10:02:10,527 - phyluce_snp_phase_uces - INFO - ================ Starting phyluce_snp_phase_uces ================
    2018-07-27 10:02:10,528 - phyluce_snp_phase_uces - INFO - Getting input filenames and creating output directories
    2018-07-27 10:02:10,528 - phyluce_snp_phase_uces - INFO - ------------- Processing alligator_mississippiensis -------------
    2018-07-27 10:02:10,528 - phyluce_snp_phase_uces - INFO - Phasing BAM file for alligator_mississippiensis
    2018-07-27 10:02:21,695 - phyluce_snp_phase_uces - INFO - Sorting BAM for alligator_mississippiensis
    2018-07-27 10:02:23,115 - phyluce_snp_phase_uces - INFO - Sorting BAM for alligator_mississippiensis
    2018-07-27 10:02:24,533 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 0
    2018-07-27 10:02:24,583 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file 1
    2018-07-27 10:02:24,613 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTQ file unphased
    2018-07-27 10:02:24,643 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 0 from FASTQ 0
    2018-07-27 10:02:24,654 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file 1 from FASTQ 1
    2018-07-27 10:02:24,662 - phyluce_snp_phase_uces - INFO - Creating REF/ALT allele FASTA file unphased from FASTQ unphased
    2018-07-27 10:02:24,673 - phyluce_snp_phase_uces - INFO - Checking for correct FASTA files
    2018-07-27 10:02:24,673 - phyluce_snp_phase_uces - INFO - Cleaning FASTA files
    2018-07-27 10:02:24,681 - phyluce_snp_phase_uces - INFO - Balancing FASTA files
    2018-07-27 10:02:24,682 - phyluce_snp_phase_uces - INFO - Symlinking FASTA files


The program automatically produces a consensus sequence for each of these phased bam files (= allele sequence) and stores these allele sequences of all samples in a joined FASTA file (``joined_allele_sequences_all_samples.fasta``). This allele FASTA is deposited in the subfolder ``fastas`` within your output folder (e.g. ``/path/to/multialign-bams-phased-reads/fastas``).

You can directly input that file (``joined_allele_sequences_all_samples.fasta``) back into the alignment pipeline, like so:

.. code-block:: bash

    phyluce_align_seqcap_align \
        --fasta /path/to/multialign-bams-phased-reads/fastas/joined_allele_sequences_all_samples.fasta \
        --output PHASED-DATA mafft-nexus-edge-trimmed \
        --taxa 4 \
        --aligner mafft \
        --cores 12 \
        --incomplete-matrix \
        --log-path log

Following alignment, you can choose how you'd like to treat these data (e.g. internally trim, analyze, etc).
