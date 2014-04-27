.. include:: global.rst

***************
Quality control
***************

When you receive your data from the sequencer, typically they are already
demultiplexed (split by sequence tags) for you. If your data are from the
MiSeq_, they may also have been trimmed of adapter contamination.

.. note:: If your data are not demultiplexed, they can come in a vast array
   of different types, although the common type is as so-called *BCL* files.
   Regardless, if your data are not demultiplexed fastq files, you will need to
   talk to your sequencing provider about how to accomplish demultiplexing.  I
   provide an **unsupported** guide to demultiplexing *BCL* files using Casava_
   here: https://gist.github.com/brantfaircloth/3125885.

Regardless, you need to do a fair bit of quality control on your read data.
**At a minimum**, this includes:

* getting some idea of how much data you have
* trimming adapter contamination of reads
* trimming low quality bases from reads

Although the MiSeq may trim some adapter contamination, running your reads
through an additional round of trimming won't hurt.  There is also `lots of
evidence <http://scholar.google.com/scholar?q=sequence+quality+affects+short+rea
d+assembly&btnG=&hl=en&as_sdt=0%2C5>`_ showing that quality control of your read
data has a **HUGE** effect on your overall success.  Bottom line is:  *garbage
in, garbage out*.

Read Counts
===========

The first thing to do once you have your read data in hand is to to get an idea
of the split of reads among your samples.  This does two things:

#. Allows you to determine how well you split the run among your sequence tags
#. Shows you which samples may be suboptimal (have very few reads)

Really unequal read counts mean that you may want to switch up your library
quantification process (or your pooling steps, prior to enrichment). Suboptimal
read counts for particular libraries may or may not mean that the enrichments
of those samples "failed", but it's generally a reasonable indication.

You can get read counts in one of two way - the first is very simple, the
second uses phyluce_.

Count reads using shell tools
-----------------------------

You can get a quick and dirty idea of the number of reads you have for each
sample using simple shell or "terminal" tools - counting lines in lots of files
is a task that's really suited to UNIX-like operating systems.

Because FASTQ files contain 4 lines for each read, we can count all the lines
in every file and divide each by 4 to get the count of reads in a given file.

Assuming your reads are in a directory structure like::

    some-directory-name/
        sample1_R1.fastq.gz
        sample1_R2.fastq.gz
        sample2_R1.fastq.gz
        sample2_R2.fastq.gz

you can count lines across all files in all directories, using:

.. code-block:: bash

    for i in some-directory-name/*_R1.fastq.gz; do echo $i; gunzip -c $i | wc -l; done

This will output a list of results like::

    bfidt-000_R1.fastq.gz
     10000000
    bfidt-000_R2.fastq.gz
     16000000

If you **divide these numbers by 4**, then that is the count of R1 reads.  The
number of reads in the R2 files, if you have paired-end data, should **always**
be equal.

Count reads using phyluce_
--------------------------

phyluce_ has slightly more sophisticated options to summarize reads in FASTQ
files. In addition to counting the number of reads in each file, phyluce_ will
give you summary statistics about each read, and it will also summarize files
in parallel (which you can do in the shell using the `xargs` command).  To
summarize the same directory of reads using phyluce, you can run::

    $ phyluce misc fastq --by-file --fastqs some-directory-name

and you should see output similar to::

    $ phyluce misc fastq --by-file --fastqs some-directory-name
    2013-12-31 18:03:22,915 - phyluce - INFO - ======================== Starting phyluce =======================
    2013-12-31 18:03:22,916 - phyluce - INFO - Version: xxxxxx
    2013-12-31 18:03:22,916 - phyluce - INFO - Argument --by_file: True
    2013-12-31 18:03:22,916 - phyluce - INFO - Argument --cmd: fastq
    2013-12-31 18:03:22,916 - phyluce - INFO - Argument --cores: 1
    2013-12-31 18:03:22,916 - phyluce - INFO - Argument --csv: False
    2013-12-31 18:03:22,916 - phyluce - INFO - Argument --exclude: None
    2013-12-31 18:03:22,917 - phyluce - INFO - Argument --fastqs: some-directory-name
    2013-12-31 18:03:22,917 - phyluce - INFO - Argument --func: <function fastq_lengths at 0x1019d5410>
    2013-12-31 18:03:22,917 - phyluce - INFO - Argument --log_path: None
    2013-12-31 18:03:22,917 - phyluce - INFO - Argument --output: None
    2013-12-31 18:03:22,917 - phyluce - INFO - Argument --verbosity: INFO
    2013-12-31 18:03:22,917 - phyluce - INFO - ........................ Running phyluce ........................
    2013-12-31 18:03:22,918 - phyluce - INFO - Getting FASTQ files
    2013-12-31 18:03:22,919 - phyluce - INFO - Getting FASTA stats using 1 cores
    ...
    2013-12-31 18:04:03,045 - phyluce - INFO - file,reads,bp,avg_len,stderr_len,min_len,max_len,median_len,contigs>1kb
    2013-12-31 18:04:03,528 - phyluce - INFO - ecoli-S14-R2.fastq.gz,603108,138387212,229.456767279,0.0517896770301,40,251,250.0
    2013-12-31 18:04:04,042 - phyluce - INFO - ecoli-S14-R1.fastq.gz,603108,143836439,238.492009723,0.0382828666402,40,251,251.0
    2013-12-31 18:04:04,046 - phyluce - INFO - ecoli-S14-Singleton.fastq.gz,4678,1007251,215.316588286,0.823631334677,40,251,250.0
    2013-12-31 18:04:04,046 - phyluce - INFO - ======================= Completed phyluce =======================

.. attention:: **LOG files**
   A quick word on the output... what gets printed to the screen is a
   **LOG** of what happened when you ran the program.  This LOG
   is both written to the screen and simultaneously saved in a file named
   (by default) `phyluce.log`. Note that it records the date and time you ran
   the command, the details of the command you ran, and output from
   the progam that you ran. The LOG file created is always *appended to* as you
   run subsequent analyses, and in this manner, it helps you track
   **specifically** what you did to your data.

This LOG output above contains the summary data in csv-format with the header
printed directly above.  If you pass the `--output OUTPUT-file.csv`, to phyluce_
using::

    $ phyluce misc fastq --by-file --fastqs some-directory-name --output
    my-output-file

then the results (without the log output) are printed to a csv file.  As you can
see, phyluce_ gives you the more detailed stats about your data, including the:

* file name
* read count
* average read length
* stderr of read length
* min read length
* max read length
* median read length

You can pass the `--cores <integer>` option to the program to use the number of
compute cores you have available to you to speed up the counting process across
many files.

.. attention:: You should always pass phyluce the number of **physical**
   compute cores on your machine.  Some technologies (like Hyper-threading)
   make it look like you have 2X as many physicaly cores as you do.  On
   laptops, the number of phyiscal cores typical is 2, on desktops 2-4, and on
   workstations and server-class computers, it ranges from 4-32.

phyluce_ will also recursively search a file path for FASTQ files, so if you
have a directory structure like::

    dir1/
        dir2/
            sample1_L001_R1_001.fastq.gz
            sample1_L001_R2_001.fastq.gz
            dir3/
                sample1_L001_R1_001.fastq.gz
                sample1_L001_R2_001.fastq.gz
        dir4/
            sample1_L001_R1_001.fastq.gz
            sample1_L001_R2_001.fastq.gz

It *should* find all the files and report back on their lengths.

.. warning:: Be careful because this recursive option can cause problems when
   folder hierarchies are very deep or very complex.

Directory and file organization
===============================

Your demultiplexed

Adapter- and quality- trimming
==============================


    Project_name/
        Sample_BFIDT-000
            BFIDT-000_AGTATAGCGC_L001_R1_001.fastq.gz
            BFIDT-000_AGTATAGCGC_L001_R2_001.fastq.gz
        Sample_BFIDT-001
            BFIDT-001_TTGTTGGCGG_L001_R1_001.fastq.gz
            BFIDT-001_TTGTTGGCGG_L001_R2_001.fastq.gz

Now, what you want to do is to clean these reads of adapter contamination
and trim low-quality bases from all reads (and probably also drop reads
containing "N" (ambiguous) bases.  Then you want to interleave the resulting
data, where read pairs are maintained, and also have an output file of singleton
data, where read pairs are not.

You can do this however you like. For the steps that follow, you want a
resulting directory structure that looks like (replace genus_species with your
taxon names)::

    genus_species1/
        interleaved-adapter-quality-trimmed/
            genus_species1-READ1and2-interleaved.fastq.gz
            genus_species1-READ-singleton.fastq.gz
    genus_species2/
        interleaved-adapter-quality-trimmed/
            genus_species2-READ1and2-interleaved.fastq.gz
            genus_species2-READ-singleton.fastq.gz

The "READ1and2" file should be in interleaved file containig both reads kept,
and the READ file should contain those data where one read of the pair are kept
(AKA singleton reads).

Trimming with illumiprocessor
-----------------------------

I created a program that I use for adapter and quality trimming named
illumiprocessor_. It generally automates these processes and produces output
if the format we want downstream.

If you used Casava, it will be easiest to place all of the demuliplexed reads
into a single directory.  If you used `splitaake`_, then things should be all
set.

You need to generate a configuration file `your-illumiprocessor.conf`, that
gives details of your reads, how you want them processed, and what renaming
options to use. This file is an extension of the INI file used for `splitaake`_
and it looks like this::


    [adapters]
    truseq1:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTGAAAAA
    truseq2:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA

    [indexes]
    bfidt-000:AACCGAGTTA
    bfidt-001:AATACTTCCG
    bfidt-002:AATTAAGGCC
    bfidt-003:AAGCTTATCC

    [params]
    separate reads:True
    read1:{name}_R1.fastq.gz
    read2:{name}_R2.fastq.gz

    [combos]
    genus_species1:bfidt-000
    genus_species2:bfidt-001
    genus_species3:bfidt-002
    genus_species4:bfidt-003

    [remap]
    bfidt-000:genus_species1
    bfidt-001:genus_species2
    bfidt-002:genus_species3
    bfidt-003:genus_species4

This gives the adapter sequences to trim (with the index indicated by an
asterisk) in `[adapters]`, the indexes we used `[indexes]`, whether reads are
separate and a name-formatting convention `[params]`, a mapping of species
names to index `[combos]`, and a mapping of index names to species `[remap]`.

You can run illumiprocessor against your data (in `demultiplexed`) with the
following.  If you do not have a multicore machine, you may with to run with
`--cores=1`.  Additionally, multicore operations require a fair amount of RAM,
so if you're low on RAM, run with fewer cores:

.. code-block:: bash

    mkdir uce-clean
    python ~/git/illumiprocessor/illumiprocessor.py \
        demultiplexed \
        uce-clean \
        your-illumiprocesser.conf \
        --remap \
        --clean \
        --cores 12 \
        --complex

The clean data will appear in `uce-clean` with the following structure::

    uce-clean/
        genus_species1/
            interleaved-adapter-quality-trimmed/
                genus_species-READ1and2-interleaved.fastq.gz
                genus_species-READ-singleton.fastq.gz
            stats/
                genus_species-READ1.fastq.gz-adapter-contam.txt
                genus_species--READ2.fastq.gz-adapter-contam.txt
                sickle-trim.txt
            untrimmed/
                genus_species-READ1.fastq.gz (symlink)
                genus_species-READ1.fastq.gz (symlink)
        genus_species2/
            interleaved-adapter-quality-trimmed/
                genus_species-READ1and2-interleaved.fastq.gz
                genus_species-READ-singleton.fastq.gz
            stats/
                genus_species-READ1.fastq.gz-adapter-contam.txt
                genus_species--READ2.fastq.gz-adapter-contam.txt
                sickle-trim.txt
            untrimmed/
                genus_species-READ1.fastq.gz (symlink)
                genus_species-READ1.fastq.gz (symlink)

`interleaved-adapter-quality-trimmed` contains the cleaned read data in
interleaved format, with one file containing "READ1and2" (both reads kept) and
another file containing "READ" data, where singleton reads are kept (singletons
sometimes result from the QC routines).
