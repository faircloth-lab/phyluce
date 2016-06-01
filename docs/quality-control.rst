.. include:: global.rst

.. _Quality Control:

***************
Quality Control
***************

When you receive your data from the sequencer, typically they are already
demultiplexed (split by sequence tags) for you. If your data are from the
MiSeq_, they may also have been trimmed of adapter contamination.

.. note:: If your data are not demultiplexed, they can come in a vast array
   of different types, although one common type are so-called BCL-files.
   Regardless, if your data are not demultiplexed fastq files, you will need to
   talk to your sequencing provider about how to accomplish demultiplexing.  I
   provide an **unsupported** guide to demultiplexing *BCL* files using Casava_
   or bcl2fastq here: https://gist.github.com/brantfaircloth/3125885.

Regardless, you need to do a fair bit of quality control on your read data.
**At a minimum**, this includes:

* getting some idea of how much data you have
* trimming adapter contamination of reads
* trimming low quality bases from reads

Although the MiSeq may trim some adapter contamination, running your reads
through an additional round of trimming won't hurt.  There is also `lots of
evidence <http://scholar.google.com/scholar?q=sequence+quality+affects+short+rea
d+assembly&btnG=&hl=en&as_sdt=0%2C5>`_ showing that quality control of your read
data has a **large** effect on your overall success.  Bottom line is:  *garbage
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

Get read counts using phyluce_
------------------------------

There is a bit of code in phyluce that lets you count reads.  To use it, you
pass the code the directory containing reads you want to summarize and run:

.. code-block:: bash

    phyluce_assembly_get_fastq_lengths --input /directory/containing/reads/ --csv; done

You can run this across many directories of reads as described in
:ref:`Tutorial I`.


Adapter- and quality-trimming
==============================

Generally speaking, Casava_ and bcl2fastq_, the two programs used to demultiplex
sequence data from Illumina_ platforms, output fastq files in a format similar
to the following::

    Project_name/
        Sample_BFIDT-000
            BFIDT-000_AGTATAGCGC_L001_R1_001.fastq.gz
            BFIDT-000_AGTATAGCGC_L001_R2_001.fastq.gz
        Sample_BFIDT-001
            BFIDT-001_TTGTTGGCGG_L001_R1_001.fastq.gz
            BFIDT-001_TTGTTGGCGG_L001_R2_001.fastq.gz

What you want to do is to clean these reads of adapter contamination
and trim low-quality bases from all reads (and probably also drop reads
containing "N" (ambiguous) bases.  Then you want to interleave the resulting
data, where read pairs are maintained, and also have an output file of singleton
data, where read pairs are not.

You can do this however you like. phyluce_ assumes that the results structure of
your data after trimming will look like the following (replace genus_species
with your taxon names)::

    genus_species1/
        split-adapter-quality-trimmed/
            genus_species1-READ1.fastq.gz
            genus_species1-READ2.fastq.gz
            genus_species1-READ-singleton.fastq.gz
    genus_species2/
        split-adapter-quality-trimmed/
            genus_species2-READ1.fastq.gz
            genus_species2-READ2.fastq.gz
            genus_species2-READ-singleton.fastq.gz

.. note:: We not longer create `interleaved` fastq files because most programs
    assume that PE (paired-end) reads are maintained as read pairs.  Rather than
    create an interleaved set of data and a non-interleaved set of data, we
    have decided to keep files for each of R1, R2, and singleton reads (those
    that have lost their mates) and never produce an interleaved file.

Trimming with illumiprocessor
-----------------------------

You can run your adapter and quality trimming and output the files in the
correct format using a program I wrote called illumiprocessor_. It automates
these processes over hundred of files and produces output in the format we want
downstream.

.. attention:: Prior to running illumiprocessor_, if you used Casava, it will
    be easiest to place all of the demuliplexed reads into a single directory.

You need to generate a configuration file ``your-illumiprocessor.conf``, that
gives details of your reads, how you want them processed, and what renaming
options to use. There are **several** variations in formatting required
depending on the library preparation method that you used.

.. attention:: See the docuementation for illumiprocessor_ for configuration
    information.

You can run illumiprocessor against your data (in `demultiplexed`) with the
following.  If you do not have a multicore machine, you may with to run with
``--cores=1``.  Additionally, multicore operations require a fair amount of RAM,
so if you're low on RAM, run with fewer cores:

.. code-block:: bash

    illumiprocessor \
        --input demultiplexed \
        --output uce-clean \
        --config your-illumiprocesser.conf \
        --cores 12

The clean data will appear in ``uce-clean`` with the following structure::

    uce-clean/
        genus_species1/
            adapters.fasta
            raw-reads/
                genus_species1-READ1.fastq.gz (symlink)
                genus_species1-READ2.fastq.gz (symlink)
            split-adapter-quality-trimmed/
                genus_species1-READ1.fastq.gz
                genus_species1-READ2.fastq.gz
                genus_species1-READ-singleton.fastq.gz
            stats/
                genus_species1-adapter-contam.txt

        genus_species2/
            adapters.fasta
            raw-reads/
                genus_species2-READ1.fastq.gz (symlink)
                genus_species2-READ2.fastq.gz (symlink)
            split-adapter-quality-trimmed/
                genus_species2-READ1.fastq.gz
                genus_species2-READ2.fastq.gz
                genus_species2-READ-singleton.fastq.gz
            stats/
                genus_species2-adapter-contam.txt

You are now ready to move onto processing the cleaned read data.
