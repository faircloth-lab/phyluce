.. include:: global.rst

.. _Tutorial I:

*****************************
Tutorial I: UCE Phylogenomics
*****************************

In the following example, we are going to process raw read data from UCE
enrichments performed against several divergent taxa so that you can get a feel
for how a typical analysis goes.  I'm also going to use several tricks that I
did not cover in the :ref:`UCE Processing` section.

The taxa we are working with will be:

* Mus musculus (PE100)
* Anolis carolinensis (PE100)
* Alligator mississippiensis (PE150)
* Gallus gallus (PE250)

Download the data
=================

You can download the data from figshare
(http://dx.doi.org/10.6084/m9.figshare.1284521).  If you want to use the command
line, you can use something like:

.. code-block:: bash

    # create a project directory
    mkdir uce-tutorial

    # change to that directory
    cd uce-tutorial

    # download the data into a file names fastq.zip
    wget -O fastq.zip https://ndownloader.figshare.com/articles/1284521/versions/1

    # make a directory to hold the data
    mkdir raw-fastq

    # move the zip file into that directory
    mv fastq.zip raw-fastq

    # move into the directory we just created
    cd raw-fastq

    # unzip the fastq data
    unzip fastq.zip

    # delete the zip file
    rm fastq.zip

    # you should see 6 files in this directory now
    ls -l

    -rw-r--r--. 1 bcf data 152M Apr 11 14:03 Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf data 147M Apr 11 14:03 Alligator_mississippiensis_GGAGCTATGG_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf data 173M Apr 11 14:03 Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf data 174M Apr 11 14:03 Anolis_carolinensis_GGCGAAGGTT_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf data  67M Apr 11 14:03 Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf data  73M Apr 11 14:03 Gallus_gallus_TTCTCCTTCA_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf data 122M Apr 11 14:03 Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf data 121M Apr 11 14:03 Mus_musculus_CTACAACGGC_L001_R2_001.fastq.gz

Alternatively, if you think of the filesystem as a tree-like structure, the
directory in which we are working (`uce-tutorial`) would look like:

.. code-block:: bash

    uce-tutorial
    └── raw-fastq
        ├── Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
        ├── Alligator_mississippiensis_GGAGCTATGG_L001_R2_001.fastq.gz
        ├── Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
        ├── Anolis_carolinensis_GGCGAAGGTT_L001_R2_001.fastq.gz
        ├── Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
        ├── Gallus_gallus_TTCTCCTTCA_L001_R2_001.fastq.gz
        ├── Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
        └── Mus_musculus_CTACAACGGC_L001_R2_001.fastq.gz


If you do not want to use the command line, you can download the data using the
figshare interface or by clicking:

http://downloads.figshare.com/article/public/1284521

Count the read data
===================

Usually, we want a count of the actual number of reads in a given sequence file
for a given species.  As mentioned in the :ref:`UCE Processing` section, we can
do this several ways.  We'll use tools from unix, because they are fast. The
next line of code will count the lines in each R1 file (which should be equal to
the reads in the R2 file) and divide that number by 4 to get the number of
sequence reads.

.. code-block:: bash

    for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done

You should see:

.. code-block:: bash

    Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
    1750000
    Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
    1874362
    Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
    376559
    Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
    1298196

Clean the read data
===================

The data you just downloaded are actual, raw, untrimmed fastq data.  This means
they contain adapter contamination and low quality bases.  We need to remove
these - which you can do several ways.  We'll use another program that I
wrote (illumiprocessor_) because it allows us to trim many different indexed
adapters from individual-specific fastq files - something that is a pain to do
by hand.  That said, you can certainly trim your reads however you would like.
See the illumiprocessor_ website for instructions on installing the program.

To use this program, we will create a configuration file that we will use to
inform the program about which adapters are in which READ1 and READ2 files.  The
data we are trimming, here, are from TruSeq v3 libraries, but the indexes are 10
nucleotides long.  We will set up the trimming file with these parameters, but
please see the illumiprocessor_ documentation for other options.

.. code-block:: bash

    # this is the section where you list the adapters you used.  the asterisk
    # will be replaced with the appropriate index for the sample.
    [adapters]
    i7:AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
    i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

    # this is the list of indexes we used
    [tag sequences]
    BFIDT-166:GGAGCTATGG
    BFIDT-016:GGCGAAGGTT
    BFIDT-045:TTCTCCTTCA
    BFIDT-011:CTACAACGGC

    # this is how each index maps to each set of reads
    [tag map]
    Alligator_mississippiensis_GGAGCTATGG:BFIDT-166
    Anolis_carolinensis_GGCGAAGGTT:BFIDT-016
    Gallus_gallus_TTCTCCTTCA:BFIDT-045
    Mus_musculus_CTACAACGGC:BFIDT-011

    # we want to rename our read files something a bit more nice - so we will
    # rename Alligator_mississippiensis_GGAGCTATGG to alligator_mississippiensis
    [names]
    Alligator_mississippiensis_GGAGCTATGG:alligator_mississippiensis
    Anolis_carolinensis_GGCGAAGGTT:anolis_carolinensis
    Gallus_gallus_TTCTCCTTCA:gallus_gallus
    Mus_musculus_CTACAACGGC:mus_musculus

I create this file in a directory **above** the one holding my reads, so the
structure looks like:

.. code-block:: bash

    uce-tutorial
    ├── illumiprocessor.conf
    └── raw-fastq
        ├── Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
        ├── Alligator_mississippiensis_GGAGCTATGG_L001_R2_001.fastq.gz
        ├── Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
        ├── Anolis_carolinensis_GGCGAAGGTT_L001_R2_001.fastq.gz
        ├── Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
        ├── Gallus_gallus_TTCTCCTTCA_L001_R2_001.fastq.gz
        ├── Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
        └── Mus_musculus_CTACAACGGC_L001_R2_001.fastq.gz

Now I run illumiprocessor_ against the data.  Note that I am using **4 physical
CPU cores** to do this work.  You need to use the number of physical cores
available on your machine.

.. code-block:: bash

    # go to the directory containing our config file and data
    cd my-analysis

    # run illumiprocessor

    illumiprocessor \
        --input raw-fastq/ \
        --output clean-fastq \
        --config illumiprocessor.conf \
        --cores 4

The output should look like the following:

.. code-block:: bash

    2015-04-11 14:23:57,912 - illumiprocessor - INFO - ==================== Starting illumiprocessor ===================
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Version: 2.0.6
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Argument --config: illumiprocessor.conf
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Argument --cores: 4
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Argument --input: /scratch/bfaircloth-uce-tutorial/raw-fastq
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Argument --log_path: None
    2015-04-11 14:23:57,913 - illumiprocessor - INFO - Argument --min_len: 40
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --no_merge: False
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --output: /scratch/bfaircloth-uce-tutorial/clean-fastq
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --phred: phred33
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --r1_pattern: None
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --r2_pattern: None
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --se: False
    2015-04-11 14:23:57,914 - illumiprocessor - INFO - Argument --trimmomatic: /home/bcf/anaconda/jar/trimmomatic.jar
    2015-04-11 14:23:57,915 - illumiprocessor - INFO - Argument --verbosity: INFO
    2015-04-11 14:23:57,942 - illumiprocessor - INFO - Trimming samples with Trimmomatic
    Running....
    2015-04-11 14:25:17,714 - illumiprocessor - INFO - =================== Completed illumiprocessor ===================

Notice that the program has created a `log` file showing what it did, and it has
also created a new directory holding the clean data that has the name `clean-
fastq` (what you told it to name the directory). Within that new directory,
there are taxon-specific folder for the cleaned reads. More specifically, your
directory structure should look similar to the following (I've collapsed the
list of raw-reads):

.. code-block:: bash

    uce-tutorial
    ├── clean-fastq
    │   ├── alligator_mississippiensis
    │   ├── anolis_carolinensis
    │   ├── gallus_gallus
    │   └── mus_musculus
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    └── raw-fastq

Within each organism specific directory, there are more files and folders:

.. code-block:: bash

    uce-tutorial
    ├── clean-fastq
    │   ├── alligator_mississippiensis
    │   │   ├── adapters.fasta
    │   │   ├── raw-reads
    │   │   ├── split-adapter-quality-trimmed
    │   │   └── stats
    │   ├── anolis_carolinensis
    │   │   ├── adapters.fasta
    │   │   ├── raw-reads
    │   │   ├── split-adapter-quality-trimmed
    │   │   └── stats
    │   ├── gallus_gallus
    │   │   ├── adapters.fasta
    │   │   ├── raw-reads
    │   │   ├── split-adapter-quality-trimmed
    │   │   └── stats
    │   └── mus_musculus
    │       ├── adapters.fasta
    │       ├── raw-reads
    │       ├── split-adapter-quality-trimmed
    │       └── stats
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    └── raw-fastq

And, within each of those directories nested within the species-specific
directory, there are additional files or links to files:

.. code-block:: bash

    uce-tutorial
    ├── clean-fastq
    │   ├── alligator_mississippiensis
    │   │   ├── adapters.fasta
    │   │   ├── raw-reads
    │   │   │   ├── alligator_mississippiensis-READ1.fastq.gz -> <PATH>
    │   │   │   └── alligator_mississippiensis-READ2.fastq.gz -> <PATH>
    │   │   ├── split-adapter-quality-trimmed
    │   │   │   ├── alligator_mississippiensis-READ1.fastq.gz
    │   │   │   ├── alligator_mississippiensis-READ2.fastq.gz
    │   │   │   └── alligator_mississippiensis-READ-singleton.fastq.gz
    │   │   └── stats
    │   │       └── alligator_mississippiensis-adapter-contam.txt
    │   ├── anolis_carolinensis
    │   ├── gallus_gallus
    │   └── mus_musculus
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    └── raw-fastq

I have collapsed the listing to show only the first taxon.

The `->` in the `raw-reads` directory above means there are symlinks_ to the
files. I have removed the file paths and replaced them with `<PATH>` so that the
figure will fit on a page.

The really important information is in the `split-adapter-quality-trimmed`
directory - which now holds our reads that have had adapter-contamination and
low-quality bases removed. Within this `split-adapter-quality-trimmed`
directory,  the `READ1` and `READ2` files hold reads that remain in a pair (the
reads are in the same consecutive order in each file).  The `READ-singleton`
file holds READ1 reads **OR** READ2 reads that lost their "mate" or "paired-
read" because of trimming or removal.

Quality control
---------------

You might want to get some idea of what effect the trimming has on read counts
and overall read lengths. There are certainly other (better) tools out there to
do this (like FastQC_), but you can get a reasonable idea of how good your reads
are by running the following, which will output a CSV listing of read stats by
sample:

.. code-block:: bash

    # move to the directory holding our cleaned reads
    cd clean-fastq/

    # run this script against all directories of reads

    for i in *;
    do
        phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
    done

The output you see should look like this:

.. code-block:: bash

    # sample,reads,total bp,mean length, 95 CI length,min,max,median
    All files in dir with alligator_mississippiensis-READ1.fastq.gz,3279362,294890805,89.9232243955,0.01003996813,40,100,100.0
    All files in dir with anolis_carolinensis-READ2.fastq.gz,3456457,314839345,91.0873026917,0.00799863728974,40,100,100.0
    All files in dir with gallus_gallus-READ2.fastq.gz,749026,159690692,213.197795537,0.0588973605567,40,251,250.0
    All files in dir with mus_musculus-READ-singleton.fastq.gz,2332785,211828511,90.8049867433,0.0102813002698,40,100,100.0

Now, we're ready to assemble our reads.

.. _tutorial-assembly:

Assemble the data
=================

phyluce_ has a lot of options for assembly - you can use velvet_, abyss_, or
trinity_.  For this tutorial, we are going to use trinity_, because I believe it
works best for most purposes.  The helper programs for the other assemblers use
the same config file, too, so you can easily experiment with all of the
assemblers.

To run an assembly, we need to create a another configuration file.  The
assembly configuration file looks like the following, assuming we want to
assemble all of our data from the organisms above:

.. code-block:: bash

    [samples]
    alligator_mississippiensis:/scratch/bfaircloth-uce-tutorial/clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
    anolis_carolinensis:/scratch/bfaircloth-uce-tutorial/clean-fastq/anolis_carolinensis/split-adapter-quality-trimmed/
    gallus_gallus:/scratch/bfaircloth-uce-tutorial/clean-fastq/gallus_gallus/split-adapter-quality-trimmed/
    mus_musculus:/scratch/bfaircloth-uce-tutorial/clean-fastq/mus_musculus/split-adapter-quality-trimmed/

We will save this into a file named `assembly.conf` at the top of our
`uce-tutorial` directory:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    └── trinity-assemblies

If you want to change the names on the left hand side of the colon in the config
file, you can do so, but the `PATHs` on the right hand side need to point to our
"clean" UCE raw reads.  If you have files in multiple locations, you can use
different `PATHs` on the right-hand side.

.. attention:: Although you can easily input new PATHs in this file, the
    **structure** of the data below the PATH you use must be the same - meaning
    that the structure and naming scheme for READ1, READ2, and READ-singleton
    must be the same.  Or, put another way, it must look like the following:

    .. code-block:: bash

        some-random-data
        ├── clean-fastq
            ├── alligator_mississippiensis
            │    └── split-adapter-quality-trimmed
            │       ├── alligator_mississippiensis-READ1.fastq.gz
            │       ├── alligator_mississippiensis-READ2.fastq.gz
            │       └── alligator_mississippiensis-READ-singleton.fastq.gz
            └── anolis_carolinensis
                └── split-adapter-quality-trimmed
                    ├── anolis_carolinensis-READ1.fastq.gz
                    ├── anolis_carolinensis-READ2.fastq.gz
                    └── anolis_carolinensis-READ-singleton.fastq.gz

Now that we have that file created, copy it to our working directory, and run
the `phyluce_assembly_assemblo_trinity` program:

.. code-block:: bash

    # make sure we are at the top-level of our uce tutorial directory
    cd uce-tutorial

    # run the assembly
    phyluce_assembly_assemblo_trinity \
        --conf assembly.conf \
        --output trinity-assemblies \
        --clean \
        --cores 12

.. warning:: Note that I am using 12 physical CPU cores to do this work.  You
    need to use the number of physical cores available on *your* machine. The
    `phyluce.conf` file also assumes you have **at least** 8 GB of RAM on your
    system, and it is better to have much more.  If you use more CPU cores than
    you have or you specify more RAM than you have, bad things will happen.

As the assembly proceeds, you should see output similar to the following:

.. code-block:: bash

    2015-04-11 18:22:41,128 - phyluce_assembly_assemblo_trinity - INFO - -------------------- Processing gallus_gallus ----
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --dir: None
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --log_path: None
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --min_kmer_coverage: 2
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --output: /scratch/bfaircloth-uce-tutorial/trinity-assemblies
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --subfolder:
    2015-04-11 15:30:54,183 - phyluce_assembly_assemblo_trinity - INFO - Argument --verbosity: INFO
    2015-04-11 15:30:54,184 - phyluce_assembly_assemblo_trinity - INFO - Getting input filenames and creating output directories
    2015-04-11 15:30:54,186 - phyluce_assembly_assemblo_trinity - INFO - ------------- Processing alligator_mississippiensis -------------
    2015-04-11 15:30:54,186 - phyluce_assembly_assemblo_trinity - INFO - Finding fastq/fasta files
    2015-04-11 15:30:54,189 - phyluce_assembly_assemblo_trinity - INFO - File type is fastq
    2015-04-11 15:30:54,190 - phyluce_assembly_assemblo_trinity - INFO - Copying raw read data to /scratch/bfaircloth-uce-tutorial/trinity-assemblies/alligator_mississippiensis_trinity
    2015-04-11 15:30:54,561 - phyluce_assembly_assemblo_trinity - INFO - Combining singleton reads with R1 data
    2015-04-11 15:30:54,576 - phyluce_assembly_assemblo_trinity - INFO - Running Trinity.pl for PE data
    2015-04-11 16:29:19,127 - phyluce_assembly_assemblo_trinity - INFO - Removing extraneous Trinity files
    2015-04-11 16:29:20,957 - phyluce_assembly_assemblo_trinity - INFO - Symlinking assembled contigs into /scratch/bfaircloth-uce-tutorial/trinity-assemblies/contigs
    2015-04-11 16:29:20,957 - phyluce_assembly_assemblo_trinity - INFO - ----------------- Processing anolis_carolinensis ----------------
    ...[continued]...
    2015-04-11 22:30:19,558 - phyluce_assembly_assemblo_trinity - INFO - ========== Completed phyluce_assembly_assemblo_trinity ==========

.. attention:: This is not a toy tutorial - the data you downloaded contain
    roughly 100 MB for each of READ1 and READ2, which means it's going to take
    some time for these samples to assemble.

One the assembly is finished, have a look at the directory structure:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    └── trinity-assemblies
        ├── alligator_mississippiensis_trinity
        │   ├── contigs.fasta -> Trinity.fasta
        │   ├── Trinity.fasta
        │   └── trinity.log
        ├── anolis_carolinensis_trinity
        │   ├── contigs.fasta -> Trinity.fasta
        │   ├── Trinity.fasta
        │   └── trinity.log
        ├── contigs
        │   ├── alligator_mississippiensis.contigs.fasta -> ../alligator_mississippiensis_trinity/Trinity.fasta
        │   ├── anolis_carolinensis.contigs.fasta -> ../anolis_carolinensis_trinity/Trinity.fasta
        │   ├── gallus_gallus.contigs.fasta -> ../gallus_gallus_trinity/Trinity.fasta
        │   └── mus_musculus.contigs.fasta -> ../mus_musculus_trinity/Trinity.fasta
        ├── gallus_gallus_trinity
        │   ├── contigs.fasta -> Trinity.fasta
        │   ├── Trinity.fasta
        │   └── trinity.log
        └── mus_musculus_trinity
            ├── contigs.fasta -> Trinity.fasta
            ├── Trinity.fasta
            └── trinity.log

Your species-specific assembly files are in the `trinity-assemblies` directory
nested within species-specific directories that correspond to the name you used
in the `assembly.conf` file (to the left of the colon).  Each name is appended
with `_trinity` because that's what trinity_ requires.  There is a symlink_
within each species-specific folder so that you now the "contigs" you assembled
are in `Trinity.fasta`.  `trinity.log` holds the log output from trinity_.

There is also a `contigs` directory within this folder.  The `contigs` directory
is the important one, because it contains symlinks_ to all of the species-
specific contigs.  This means that you can treat this single folder as if it
contains all of your assembled contigs.

Assembly QC
-----------

We can get a sense of how well the assembly worked by running the following from
the top of our working directory:

.. code-block:: bash

    # run this script against all directories of reads

    for i in trinity-assemblies/contigs/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done

This should output something similar to the following.  I've added the header as
a comment:

.. code-block:: bash

    # samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
    alligator_mississippiensis.contigs.fasta,10587,5820479,549.776046094,3.5939422934,224,11285,413.0,1182
    anolis_carolinensis.contigs.fasta,2458,1067208,434.177379984,5.72662897806,224,4359,319.0,34
    gallus_gallus.contigs.fasta,19905,8841661,444.192966591,2.06136172068,224,9883,306.0,1530
    mus_musculus.contigs.fasta,2162,1126231,520.920906568,7.75103292163,224,6542,358.0,186

.. admonition:: Question: Why are my numbers slightly different than your numbers?
    :class: admonition tip

    The process of read assembly often differs by operating system and sometimes
    by OS version, and some of these differences are due to libraries that
    underlie many of the assembly programs.  Expect to see differences.  You
    should not expect for them to be huge.

.. attention:: If you see max-contig sizes around 16KB (for vertebrates), that
    is commonly the entire or almost-entire mtDNA genome.  You do not tend to
    see entire mtDNA assemblies when the input DNA was extracted from a source
    having few mitochondria (e.g. blood).

There are many, many other assembly QC steps you can run other than simply
looking at the stats of the assembled contigs.  We will not go into those here.

Finding UCE loci
================

Now that we've assembled our contigs from raw reads, it's time to find those
contigs which are UCE loci and move aside those that are not.  The directory
structure before we do this should look like the following:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    └── trinity-assemblies

Before we locate UCE loci, you need to get the probe set used for the
enrichments:

.. code-block:: bash

    wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta

Now, our directory structure looks like:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    ├── trinity-assemblies
    └── uce-5k-probes.fasta

Now, run the `phyluce_assembly_match_contigs_to_probes` program:

.. code-block:: bash

    phyluce_assembly_match_contigs_to_probes \
        --contigs trinity-assemblies/contigs \
        --probes uce-5k-probes.fasta \
        --output uce-search-results

You should see output similar to the following (also stored in
`phyluce_assembly_assemblo_trinity.log`):

.. code-block:: bash

    2015-04-12 12:44:58,951 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Starting phyluce_assembly_match_contigs_to_probes =======
    2015-04-12 12:44:58,951 - phyluce_assembly_match_contigs_to_probes - INFO - Version: git 2a9c49d
    2015-04-12 12:44:58,951 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --contigs: /scratch/uce-tutorial/trinity-assemblies/contigs
    2015-04-12 12:44:58,951 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --dupefile: None
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --keep_duplicates: None
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --log_path: None
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_coverage: 80
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_identity: 80
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --output: /scratch/uce-tutorial/uce-search-results
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --probes: /scratch/uce-tutorial/uce-5k-probes.fasta
    2015-04-12 12:44:58,952 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --regex: ^(uce-\d+)(?:_p\d+.*)
    2015-04-12 12:44:58,953 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --verbosity: INFO
    2015-04-12 12:44:59,111 - phyluce_assembly_match_contigs_to_probes - INFO - Creating the UCE-match database
    2015-04-12 12:44:59,405 - phyluce_assembly_match_contigs_to_probes - INFO - Processing contig data
    2015-04-12 12:44:59,405 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2015-04-12 12:45:12,735 - phyluce_assembly_match_contigs_to_probes - INFO - alligator_mississippiensis: 4315 (40.76%) uniques of 10587 contigs, 0 dupe probe matches, 230 UCE loci removed for matching multiple contigs, 40 contigs removed for matching multiple UCE loci
    2015-04-12 12:45:15,919 - phyluce_assembly_match_contigs_to_probes - INFO - anolis_carolinensis: 703 (28.60%) uniques of 2458 contigs, 0 dupe probe matches, 138 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    2015-04-12 12:45:33,479 - phyluce_assembly_match_contigs_to_probes - INFO - gallus_gallus: 3923 (19.71%) uniques of 19905 contigs, 0 dupe probe matches, 625 UCE loci removed for matching multiple contigs, 47 contigs removed for matching multiple UCE loci
    2015-04-12 12:45:36,604 - phyluce_assembly_match_contigs_to_probes - INFO - mus_musculus: 825 (38.16%) uniques of 2162 contigs, 0 dupe probe matches, 93 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci
    2015-04-12 12:45:36,604 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2015-04-12 12:45:36,605 - phyluce_assembly_match_contigs_to_probes - INFO - The LASTZ alignments are in /scratch/uce-tutorial/uce-search-results
    2015-04-12 12:45:36,605 - phyluce_assembly_match_contigs_to_probes - INFO - The UCE match database is in /scratch/uce-tutorial/uce-search-results/probe.matches.sqlite
    2015-04-12 12:45:36,605 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Completed phyluce_assembly_match_contigs_to_probes ======

The header info at the top tells us exactly what version of the code we are
running and keeps track of our options.  The important output is:

.. code-block:: bash

    alligator_mississippiensis: 4315 (40.76%) uniques of 10587 contigs, 0 dupe probe matches, 230 UCE loci removed for matching multiple contigs, 40 contigs removed for matching multiple UCE loci
    anolis_carolinensis: 703 (28.60%) uniques of 2458 contigs, 0 dupe probe matches, 138 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    gallus_gallus: 3923 (19.71%) uniques of 19905 contigs, 0 dupe probe matches, 625 UCE loci removed for matching multiple contigs, 47 contigs removed for matching multiple UCE loci
    mus_musculus: 825 (38.16%) uniques of 2162 contigs, 0 dupe probe matches, 93 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci

Which we can break down to the following (for `alligator_mississippiensis`):

.. code-block:: bash

    alligator_mississippiensis:
        4315 (40.76%) uniques of 10587 contigs
        0 dupe probe matches
        230 UCE loci removed for matching multiple contigs
        40 contigs removed for matching multiple UCE loci

These are the capture data for the `alligator_mississippiensis` sample.  We
targeted 5k UCE loci in this sample and recovered roughly 4484 of those loci.
Before reaching that total of 4484 loci, we removed 45 and 44 loci from the data
set because they looked like duplicates (probes supposedly targeting different
loci hit the same contig or two supposedly different contigs hit probes
designed for a single UCE locus).

.. admonition:: Question: Why is the count of UCE loci different by sample?
    :class: admonition tip

    For these example data, we enriched some samples (alligator_mississippiensis
    and gallus_gallus) for 5k UCE loci, while we enriched others
    (anolis_carolinensis and mus_musculus) for 2.5k UCE loci.  Additionally, the
    2.5k UCE enrichments did not work very well (operator error).

The directory structure now looks like the following (everything collapsed but
the `uce-search-results directory`):

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── phyluce_assembly_match_contigs_to_probes.log
    ├── raw-fastq
    ├── trinity-assemblies
    ├── uce-5k-probes.fasta
    └── uce-search-results
        ├── alligator_mississippiensis.contigs.lastz
        ├── anolis_carolinensis.contigs.lastz
        ├── gallus_gallus.contigs.lastz
        ├── mus_musculus.contigs.lastz
        └── probe.matches.sqlite

The search we just ran created lastz_ search result files for each taxon, and
stored summary results of these searches in the `probe.matches.sqlite` database
(see :ref:`Database` for more information on this database and its structure).

.. _tutorial-uce-extraction:

Extracting UCE loci
===================

Now that we have located UCE loci, we need to determine which taxa we want in
our analysis, create a list of those taxa, and then generated a list of which
UCE loci we enriched in each taxon (the "data matrix configuration file").  We
will then use this list to extract FASTA data for each taxon for each UCE locus.

First, we need to decide which taxa we want in our "taxon set".  So, we create a
configuration file like so:

.. code-block:: bash

    [all]
    alligator_mississippiensis
    anolis_carolinensis
    gallus_gallus
    mus_musculus

These names need to match the assembly names we used.  Here, we have just put
all 4 taxa in a list that we named `all`.  However, we can adjust this list in
many ways (see :ref:`locus-counts`).

Save this file as `taxon-set.conf` at the top level of our `uce-tutorial`
directory.  The directory should look like this, now:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── phyluce_assembly_match_contigs_to_probes.log
    ├── raw-fastq
    ├── taxon-set.conf
    ├── trinity-assemblies
    ├── uce-5k-probes.fasta
    └── uce-search-results

Now that we have this file created, we do the following to create the initial list of loci for each
taxon:

.. code-block:: bash

    # create an output directory for this taxon set - this just keeps
    # things neat
    cd uce-tutorial
    mkdir -p taxon-sets/all

    # create the data matrix configuration file
    phyluce_assembly_get_match_counts \
        --locus-db uce-search-results/probe.matches.sqlite \
        --taxon-list-config taxon-set.conf \
        --taxon-group 'all' \
        --incomplete-matrix \
        --output taxon-sets/all/all-taxa-incomplete.conf

The output should look like the following:

.. code-block:: bash

    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - =========== Starting phyluce_assembly_get_match_counts ==========
    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - Version: git 2a9c49d
    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - Argument --extend_locus_db: None
    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - Argument --incomplete_matrix: True
    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - Argument --keep_counts: False
    2015-04-12 12:56:11,835 - phyluce_assembly_get_match_counts - INFO - Argument --locus_db: /scratch/uce-tutorial/uce-search-results/probe.matches.sqlite
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --log_path: None
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --optimize: False
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.conf
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --random: False
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --sample_size: 10
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --samples: 10
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --silent: False
    2015-04-12 12:56:11,836 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_group: all
    2015-04-12 12:56:11,837 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_list_config: /scratch/uce-tutorial/taxon-set.conf
    2015-04-12 12:56:11,837 - phyluce_assembly_get_match_counts - INFO - Argument --verbosity: INFO
    2015-04-12 12:56:11,837 - phyluce_assembly_get_match_counts - INFO - There are 4 taxa in the taxon-group '[all]' in the config file taxon-set.conf
    2015-04-12 12:56:11,838 - phyluce_assembly_get_match_counts - INFO - Getting UCE names from database
    2015-04-12 12:56:11,844 - phyluce_assembly_get_match_counts - INFO - There are 5041 total UCE loci in the database
    2015-04-12 12:56:11,946 - phyluce_assembly_get_match_counts - INFO - Getting UCE matches by organism to generate a INCOMPLETE matrix
    2015-04-12 12:56:11,948 - phyluce_assembly_get_match_counts - INFO - There are 4710 UCE loci in an INCOMPLETE matrix
    2015-04-12 12:56:11,949 - phyluce_assembly_get_match_counts - INFO - Writing the taxa and loci in the data matrix to /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.conf
    2015-04-12 12:56:11,952 - phyluce_assembly_get_match_counts - INFO - ========== Completed phyluce_assembly_get_match_counts ==========

And, our directory structure should now look like this (collapsing all but
`taxon-sets`):

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── phyluce_assembly_get_match_counts.log
    ├── phyluce_assembly_match_contigs_to_probes.log
    ├── taxon-set.conf
    ├── taxon-sets
    │   └── all
    │       └── all-taxa-incomplete.conf
    ├── trinity-assemblies
    ├── uce-5k-probes.fasta
    └── uce-search-results

Now, we need to extract FASTA data that correspond to the loci in
`all-taxa-incomplete.conf`:

.. code-block:: bash

    # change to the taxon-sets/all directory
    cd taxon-sets/all

    # make a log directory to hold our log files - this keeps things neat
    mkdir log

    # get FASTA data for taxa in our taxon set
    phyluce_assembly_get_fastas_from_match_counts \
        --contigs ../../trinity-assemblies/contigs \
        --locus-db ../../uce-search-results/probe.matches.sqlite \
        --match-count-output all-taxa-incomplete.conf \
        --output all-taxa-incomplete.fasta \
        --incomplete-matrix all-taxa-incomplete.incomplete \
        --log-path log

The output should look something like the following:

.. code-block:: bash

    2015-04-12 12:58:40,646 - phyluce_assembly_get_fastas_from_match_counts - INFO - ===== Starting phyluce_assembly_get_fastas_from_match_counts ====
    2015-04-12 12:58:40,646 - phyluce_assembly_get_fastas_from_match_counts - INFO - Version: git a6a957a
    2015-04-12 12:58:40,646 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --contigs: /scratch/uce-tutorial/trinity-assemblies/contigs
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --extend_locus_contigs: None
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --extend_locus_db: None
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --incomplete_matrix: /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.incomplete
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --locus_db: /scratch/uce-tutorial/uce-search-results/probe.matches.sqlite
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --match_count_output: /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.conf
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.fasta
    2015-04-12 12:58:40,647 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --verbosity: INFO
    2015-04-12 12:58:40,702 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 4 taxa in the match-count-config file named all-taxa-incomplete.conf
    2015-04-12 12:58:40,716 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 4710 UCE loci in an INCOMPLETE matrix
    2015-04-12 12:58:40,717 - phyluce_assembly_get_fastas_from_match_counts - INFO - ---------Getting UCE loci for alligator_mississippiensis---------
    2015-04-12 12:58:40,773 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 4315 UCE loci for alligator_mississippiensis
    2015-04-12 12:58:40,774 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for alligator_mississippiensis
    2015-04-12 12:58:51,084 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.incomplete
    2015-04-12 12:58:51,086 - phyluce_assembly_get_fastas_from_match_counts - INFO - -------------Getting UCE loci for anolis_carolinensis------------
    2015-04-12 12:58:51,120 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 703 UCE loci for anolis_carolinensis
    2015-04-12 12:58:51,121 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for anolis_carolinensis
    2015-04-12 12:58:51,693 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.incomplete
    2015-04-12 12:58:51,701 - phyluce_assembly_get_fastas_from_match_counts - INFO - ----------------Getting UCE loci for gallus_gallus---------------
    2015-04-12 12:58:51,753 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 3923 UCE loci for gallus_gallus
    2015-04-12 12:58:51,753 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for gallus_gallus
    2015-04-12 12:59:10,959 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.incomplete
    2015-04-12 12:59:10,962 - phyluce_assembly_get_fastas_from_match_counts - INFO - ----------------Getting UCE loci for mus_musculus----------------
    2015-04-12 12:59:10,997 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 825 UCE loci for mus_musculus
    2015-04-12 12:59:10,997 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for mus_musculus
    2015-04-12 12:59:11,610 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.incomplete
    2015-04-12 12:59:11,618 - phyluce_assembly_get_fastas_from_match_counts - INFO - ==== Completed phyluce_assembly_get_fastas_from_match_counts ====

And, our directory structure should now look like this (collapsing all but
`taxon-sets`):

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── phyluce_assembly_get_match_counts.log
    ├── phyluce_assembly_match_contigs_to_probes.log
    ├── taxon-set.conf
    ├── taxon-sets
    │   └── all
    │       ├── all-taxa-incomplete.conf
    │       ├── all-taxa-incomplete.fasta
    │       ├── all-taxa-incomplete.incomplete
    │       └── log
    │           └── phyluce_assembly_get_fastas_from_match_counts.log
    └── phyluce_assembly_get_fastas_from_match_counts.log
    ├── trinity-assemblies
    ├── uce-5k-probes.fasta
    └── uce-search-results

The extracted FASTA data are in a monolithic FASTA file (all data for all
organisms) named `all-taxa-incomplete.fasta`.

Exploding the monolithic FASTA file
-----------------------------------

Lots of times we want to know individual statistics on UCE assemblies for a
given taxon.  We can do that by exploding the monolithic fasta file into a file
of UCE loci that we have enriched by taxon, then running stats on those exploded
files.  To do that, run the following:

.. code-block:: bash

    # explode the monolithic FASTA by taxon (you can also do by locus)
    phyluce_assembly_explode_get_fastas_file \
        --input all-taxa-incomplete.fasta \
        --output-dir exploded-fastas \
        --by-taxon

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

Aligning UCE loci
=================

You have lots of options when aligning UCE loci.  You can align the loci and use
those alignments with no trimming, you can edge-trim the alignments following
some algorithm, and you can end+internally trim alignments following some
algorithm. It's hard to say what is best in all situations.  When taxa are
"closely" related (< 30-50 MYA, perhaps), I think that edge-trimming alignments
is reasonable.  When the taxa you are interested in span a wider range of
divergence times (> 50 MYA), you may want to think about internal trimming.

How you accomplish you edge- or internal trimming is also a decision you need to
make. In phyluce_, we implement our edge-trimming algorithm by running the
alignment program "as-is" (i.e., without the `--no-trim`) option.  We do
internal-trimming by turning off trimming using `--no-trim`, then passing the
resulting alignments (in FASTA format) to a parallel wrapper around Gblocks_.

You also have a choice of aligner - mafft_ or muscle_ (or you can externally
align UCE loci using a tool like SATé, as well).

Generally, I would use mafft_.

Edge trimming
-------------

Edge trimming your alignments is a relatively simple matter.  You can run edge
trimming, as follows:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # align the data
    phyluce_align_seqcap_align \
        --fasta all-taxa-incomplete.fasta \
        --output mafft-nexus-edge-trimmed \
        --taxa 4 \
        --aligner mafft \
        --cores 12 \
        --incomplete-matrix \
        --log-path log

.. warning:: Note that I am using 12 physical CPU cores here.  You need to use
    the number of physical cores available on *your* machine.

The output should look like this:

.. code-block:: bash

    2015-04-12 13:01:31,919 - phyluce_align_seqcap_align - INFO - ============== Starting phyluce_align_seqcap_align ==============
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Version: git a6a957a
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --aligner: mafft
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --ambiguous: False
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --cores: 12
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --fasta: /scratch/uce-tutorial/taxon-sets/all/all-taxa-incomplete.fasta
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 13:01:31,920 - phyluce_align_seqcap_align - INFO - Argument --max_divergence: 0.2
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --min_length: 100
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --no_trim: False
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --notstrict: True
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-edge-trimmed
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --output_format: nexus
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --proportion: 0.65
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --taxa: 4
    2015-04-12 13:01:31,921 - phyluce_align_seqcap_align - INFO - Argument --threshold: 0.65
    2015-04-12 13:01:31,922 - phyluce_align_seqcap_align - INFO - Argument --verbosity: INFO
    2015-04-12 13:01:31,922 - phyluce_align_seqcap_align - INFO - Argument --window: 20
    2015-04-12 13:01:31,922 - phyluce_align_seqcap_align - INFO - Building the locus dictionary
    2015-04-12 13:01:31,922 - phyluce_align_seqcap_align - INFO - Removing ALL sequences with ambiguous bases...
    2015-04-12 13:01:34,051 - phyluce_align_seqcap_align - WARNING - DROPPED locus uce-6169. Too few taxa (N < 3).
    [many more loci dropped here]
    2015-04-12 13:01:34,596 - phyluce_align_seqcap_align - INFO - Aligning with MAFFT
    2015-04-12 13:01:34,598 - phyluce_align_seqcap_align - INFO - Alignment begins. 'X' indicates dropped alignments (these are reported after alignment)
    .................[continued]
    2015-04-12 13:02:05,847 - phyluce_align_seqcap_align - INFO - Alignment ends
    2015-04-12 13:02:05,848 - phyluce_align_seqcap_align - INFO - Writing output files
    2015-04-12 13:02:06,567 - phyluce_align_seqcap_align - INFO - ============== Completed phyluce_align_seqcap_align =============

The `.` values that you see represent loci that were aligned and succesfully
trimmed. Any `X` values that you see represent loci that were removed
because trimming reduced their length to effectively nothing.

.. attention:: The number of potential alignments dropped here is abnormally
    large becase our **sample size is so small (n=4)**.

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
    │           └── mafft-nexus-edge-trimmed
    │               ├── uce-1008.nexus
    │               ├── uce-1014.nexus
    │               ├── uce-1039.nexus
    │               ...
    │               └── uce-991.nexus
    ...
    └── uce-search-results

We can output summary stats for these alignments by running the following
program:

.. code-block:: bash

    phyluce_align_get_align_summary_data \
        --alignments mafft-nexus-edge-trimmed \
        --cores 12 \
        --log-path log

.. warning:: Note that I am using 12 physical CPU cores here.  You need to use
    the number of physical cores available on *your* machine.

The output from the program should look like:

.. code-block:: bash

    2015-04-12 13:40:56,410 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
    2015-04-12 13:40:56,410 - phyluce_align_get_align_summary_data - INFO - Version: git a6a957a
    2015-04-12 13:40:56,410 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-edge-trimmed
    2015-04-12 13:40:56,410 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
    2015-04-12 13:40:56,410 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
    2015-04-12 13:40:56,411 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 13:40:56,411 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2015-04-12 13:40:56,411 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
    2015-04-12 13:40:56,411 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
    2015-04-12 13:40:56,416 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2015-04-12 13:40:56,716 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:      913
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:    559,055
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:      612.33
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:    13.92
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:       157
    2015-04-12 13:40:56,717 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:       1,266
    2015-04-12 13:40:56,719 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2015-04-12 13:40:56,719 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:        3.38
    2015-04-12 13:40:56,719 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:      0.03
    2015-04-12 13:40:56,719 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:         3
    2015-04-12 13:40:56,719 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:         4
    2015-04-12 13:40:56,720 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2015-04-12 13:40:56,720 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:     9.61
    2015-04-12 13:40:56,720 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:   0.42
    2015-04-12 13:40:56,720 - phyluce_align_get_align_summary_data - INFO - [Missing] min:      0.00
    2015-04-12 13:40:56,720 - phyluce_align_get_align_summary_data - INFO - [Missing] max:      27.98
    2015-04-12 13:40:56,732 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2015-04-12 13:40:56,732 - phyluce_align_get_align_summary_data - INFO - [All characters]    1,852,854
    2015-04-12 13:40:56,732 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]       1,600,577
    2015-04-12 13:40:56,733 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2015-04-12 13:40:56,733 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]        913 alignments
    2015-04-12 13:40:56,733 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]        913 alignments
    2015-04-12 13:40:56,733 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]        913 alignments
    2015-04-12 13:40:56,733 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]        913 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]        913 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]        913 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]        346 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]        346 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]        346 alignments
    2015-04-12 13:40:56,734 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]        346 alignments
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 72,168 times
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] '?' is present 180,109 times
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 494,411 times
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 308,796 times
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 295,171 times
    2015-04-12 13:40:56,735 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 502,199 times
    2015-04-12 13:40:56,736 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========

.. attention:: Note that there are only 2 sets of counts in the `Data matrix
    completeness` section because (1) we dropped all loci having fewer than 3
    taxa and (2) that only leaves two remaining options.

The most important data here are the number of loci we have and the number of
loci in data matrices of different completeness.  The locus length stats are
also reasonably important, but they can also be misleading because edge-trimming
does not remove internal gaps that often inflate the length of alignments.

Internal trimming
-----------------

Now, let's do the same thing, but run internal trimming on the resulting
alignments.  We will do that by turning off trimming `--no-trim` and outputting
FASTA formatted alignments with `--output-format fasta`.

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # align the data - turn off trimming and output FASTA
    phyluce_align_seqcap_align \
        --fasta all-taxa-incomplete.fasta \
        --output mafft-nexus-internal-trimmed \
        --taxa 4 \
        --aligner mafft \
        --cores 12 \
        --incomplete-matrix \
        --output-format fasta \
        --no-trim \
        --log-path log

.. attention:: The number of UCE loci dropped here is abnormally large becase
    our **sample size is so small (n=4)**.

.. warning:: Note that I am using 12 physical CPU cores here.  You need to use
    the number of physical cores available on *your* machine.

The output from the program should be the roughly the same as what we saw
before.  The current directory structure should look like (I've collapsed a
number of branches in the tree):

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
    │           └── mafft-nexus-internal-trimmed
    │               ├── uce-1008.nexus
    │               ├── uce-1014.nexus
    │               ├── uce-1039.nexus
    │               ...
    │               └── uce-991.nexus
    ...
    └── uce-search-results

Now, we are going to trim these loci using Gblocks_:

.. code-block:: bash

    # run gblocks trimming on the alignments
    phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft-nexus-internal-trimmed \
        --output mafft-nexus-internal-trimmed-gblocks \
        --cores 12 \
        --log log

The output should look like this:

.. code-block:: bash

    2015-04-12 13:51:52,450 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO -  Starting phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed
    2015-04-12 13:51:52,450 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Version: git a6a957a
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b1: 0.5
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b2: 0.85
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b3: 8
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b4: 10
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --cores: 12
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --input_format: fasta
    2015-04-12 13:51:52,451 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 13:51:52,452 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2015-04-12 13:51:52,452 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --output_format: nexus
    2015-04-12 13:51:52,452 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --verbosity: INFO
    2015-04-12 13:51:52,452 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO -  Starting phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed
    2015-04-12 13:51:52,452 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Getting aligned sequences for trimming
    2015-04-12 13:51:52,460 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Alignment trimming begins.
    .................[continued]
    2015-04-12 13:51:53,063 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Alignment trimming ends
    2015-04-12 13:51:53,063 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Writing output files
    2015-04-12 13:51:53,725 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO -  Completed phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed

The `.` values that you see represent loci that were aligned and succesfully
trimmed. Andy `X` values that you see represent loci that were aligned and
trimmed so much that there was nothing left.

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
    │           ├── mafft-nexus-internal-trimmed
    │           └── mafft-nexus-internal-trimmed-gblocks
    │               ├── uce-1008.nexus
    │               ├── uce-1014.nexus
    │               ├── uce-1039.nexus
    │               ...
    │               └── uce-991.nexus
    ...
    └── uce-search-results

We can output summary stats for these alignments by running the following
program:

.. code-block:: bash

    phyluce_align_get_align_summary_data \
        --alignments mafft-nexus-internal-trimmed-gblocks \
        --cores 12 \
        --log-path log

.. warning:: Note that I am using 12 physical CPU cores here.  You need to use
    the number of physical cores available on *your* machine.

The output from the program should look like:

.. code-block:: bash

    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Version: git a6a957a
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 13:54:48,251 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2015-04-12 13:54:48,252 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
    2015-04-12 13:54:48,252 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
    2015-04-12 13:54:48,257 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:      913
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:    499,929
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:      547.57
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:    14.10
    2015-04-12 13:54:48,523 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:       101
    2015-04-12 13:54:48,524 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:       1,180
    2015-04-12 13:54:48,525 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2015-04-12 13:54:48,525 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:        3.38
    2015-04-12 13:54:48,525 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:      0.03
    2015-04-12 13:54:48,525 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:         3
    2015-04-12 13:54:48,526 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:         4
    2015-04-12 13:54:48,526 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2015-04-12 13:54:48,526 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:     0.00
    2015-04-12 13:54:48,527 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:   0.00
    2015-04-12 13:54:48,527 - phyluce_align_get_align_summary_data - INFO - [Missing] min:      0.00
    2015-04-12 13:54:48,527 - phyluce_align_get_align_summary_data - INFO - [Missing] max:      0.00
    2015-04-12 13:54:48,537 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - [All characters]    1,717,433
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]       1,621,491
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]        913 alignments
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]        913 alignments
    2015-04-12 13:54:48,538 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]        913 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]        913 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]        913 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]        913 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]        346 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]        346 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]        346 alignments
    2015-04-12 13:54:48,539 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]        346 alignments
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 95,942 times
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 498,643 times
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 314,200 times
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 300,494 times
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 508,154 times
    2015-04-12 13:54:48,540 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========

Alignment cleaning
==================

If you look at the alignments we currently have, you will notice that each
alignment contains the locus name along with the taxon name.  This is not what
we want downstream, but it does enable us to ensure the correct data went into
each alignment.  So, we need to clean our alignments. For the remainder of
this tutorial, we will work with the Gblocks_ trimmed alignments, so we will
clean those alignments:

 .. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # align the data - turn off trimming and output FASTA
    phyluce_align_remove_locus_name_from_nexus_lines \
        --alignments mafft-nexus-internal-trimmed-gblocks \
        --output mafft-nexus-internal-trimmed-gblocks-clean \
        --cores 12 \
        --log-path log

The output should be similar to:

.. code-block:: bash

    2015-04-12 14:08:19,925 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - === Starting phyluce_align_remove_locus_name_from_nexus_lines ===
    2015-04-12 14:08:19,925 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Version: git a6a957a
    2015-04-12 14:08:19,925 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2015-04-12 14:08:19,925 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --cores: 1
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --input_format: nexus
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --log_path: None
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --output_format: nexus
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --taxa: None
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Argument --verbosity: INFO
    2015-04-12 14:08:19,926 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Getting alignment files
    Running............[continued]
    2015-04-12 14:08:22,001 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - Taxon names in alignments: alligator_mississippiensis,anolis_carolinensis,gallus_gallus,mus_musculus
    2015-04-12 14:08:22,001 - phyluce_align_remove_locus_name_from_nexus_lines - INFO - === Completed phyluce_align_remove_locus_name_from_nexus_lines ==

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
    │           ├── mafft-nexus-internal-trimmed
    │           ├── mafft-nexus-internal-trimmed-gblocks
    │           └── mafft-nexus-internal-trimmed-gblocks-clean
    │               ├── uce-1008.nexus
    │               ├── uce-1014.nexus
    │               ├── uce-1039.nexus
    │               ...
    │               └── uce-991.nexus
    ...
    └── uce-search-results

Now, if you look at the alignments, you will see that the locus names are
removed.  We're ready to generate our final data matrices.

Final data matrices
===================

For the most part, I analyze 75% and 95% complete matrices, where "completeness"
for the 75% matrix means than, in a study of 100 taxa (total), all alignments
will contain at least 75 of these 100 taxa.  Similarly, for the 95% matrix, in a
study of 100 taxa, all alignments will contain 95 of these 100 taxa.

.. attention:: Notice that this metric for completeness does not pay attention
    to which taxa are in which alignments - so the 75%, above, does **not** mean
    that a given taxon will have data in 75 of 100 alignments.

To create a 75% data matrix, run the following.  Notice that the integer
following `--taxa` is the **total** number of organisms in the study.

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # the integer following --taxa is the number of TOTAL taxa
    # and I use "75p" to denote the 75% complete matrix
    phyluce_align_get_only_loci_with_min_taxa \
        --alignments mafft-nexus-internal-trimmed-gblocks-clean \
        --taxa 4 \
        --percent 0.75 \
        --output mafft-nexus-internal-trimmed-gblocks-clean-75p \
        --cores 12 \
        --log-path log

The output should look like the following:

.. code-block:: bash

    2015-04-12 14:17:37,589 - phyluce_align_get_only_loci_with_min_taxa - INFO - ======= Starting phyluce_align_get_only_loci_with_min_taxa ======
    2015-04-12 14:17:37,589 - phyluce_align_get_only_loci_with_min_taxa - INFO - Version: git a6a957a
    2015-04-12 14:17:37,589 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean
    2015-04-12 14:17:37,589 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --cores: 12
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --input_format: nexus
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --percent: 0.75
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --taxa: 4
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --verbosity: INFO
    2015-04-12 14:17:37,590 - phyluce_align_get_only_loci_with_min_taxa - INFO - Getting alignment files
    2015-04-12 14:17:37,780 - phyluce_align_get_only_loci_with_min_taxa - INFO - Copied 913 alignments of 913 total containing ≥ 0.75 proportion of taxa (n = 3)
    2015-04-12 14:17:37,780 - phyluce_align_get_only_loci_with_min_taxa - INFO - ====== Completed phyluce_align_get_only_loci_with_min_taxa ======

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
    │           ├── mafft-nexus-internal-trimmed
    │           ├── mafft-nexus-internal-trimmed-gblocks
    │           ├── mafft-nexus-internal-trimmed-gblocks-clean
    │           └── mafft-nexus-internal-trimmed-gblocks-clean-75p
    │               ├── uce-1008.nexus
    │               ├── uce-1014.nexus
    │               ├── uce-1039.nexus
    │               ...
    │               └── uce-991.nexus
    ...
    └── uce-search-results

Preparing data for RAxML and ExaML
==================================

Now that we have our `75p` data matrix completed, we can generate input files
for subsequent phylogenetic analysis.  For the most part, I favor ExaBayes_,
RAxML_, and ExaML_ (usually in that order).  Formatting our `75p` data into
a phylip file for these programs is rather easy.  To do that, run:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # build the concatenated data matrix
    phyluce_align_format_nexus_files_for_raxml \
        --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
        --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
        --charsets \
        --log-path log

The output from this program will look like:

.. code-block:: bash

    2015-04-12 14:40:52,276 - phyluce_align_format_nexus_files_for_raxml - INFO - ====== Starting phyluce_align_format_nexus_files_for_raxml ======
    2015-04-12 14:40:52,276 - phyluce_align_format_nexus_files_for_raxml - INFO - Version: git a6a957a
    2015-04-12 14:40:52,276 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --alignments: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --charsets: True
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --log_path: /scratch/uce-tutorial/taxon-sets/all/log
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --nexus: False
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --output: /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Argument --verbosity: INFO
    2015-04-12 14:40:52,277 - phyluce_align_format_nexus_files_for_raxml - INFO - Reading input alignments in NEXUS format
    2015-04-12 14:40:53,291 - phyluce_align_format_nexus_files_for_raxml - INFO - Concatenating files
    2015-04-12 14:40:54,242 - phyluce_align_format_nexus_files_for_raxml - INFO - Writing charsets to /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml/mafft-nexus-internal-trimmed-gblocks-clean-75p.charsets
    2015-04-12 14:40:54,242 - phyluce_align_format_nexus_files_for_raxml - INFO - Writing concatenated PHYLIP alignment to /scratch/uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml/mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip
    2015-04-12 14:40:54,245 - phyluce_align_format_nexus_files_for_raxml - INFO - ====== Completed phyluce_align_format_nexus_files_for_raxml =====

.. attention:: Notice that using the `--charsets` flag will output the charsets
    as well as the concatenated PHYLIP file.  You generally want these and the
    cost for them is low.  I would always run this option - even if you later
    do not use them.

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
    │           ├── mafft-nexus-internal-trimmed
    │           ├── mafft-nexus-internal-trimmed-gblocks
    │           ├── mafft-nexus-internal-trimmed-gblocks-clean
    │           ├── mafft-nexus-internal-trimmed-gblocks-clean-75p
    │           └── mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml
    │               ├── mafft-nexus-internal-trimmed-gblocks-clean-75p.charsets
    │               └── mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip
    ...
    └── uce-search-results

Usually, I will create a separate directory to hold the alignment and
input/output files for ExaBayes:

.. code-block:: bash

    cp -R mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml mafft-nexus-internal-trimmed-gblocks-clean-75p-exbayes

RAxML
-----

The above data are ready to analysis in RAxML_.  I usually run RAxML_ in two
pieces - first running the "best" tree search, then running the bootstrap
replicates, and I almost always run these analyses on the HPC.

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml

    # get two random numbers
    for i in 1 2; do echo $RANDOM; done

    # that output the following
    # 19877
    # 7175

    # run the search for the "best" ML tree
    raxmlHPC-PTHREADS-SSE3 \
        -m GTRGAMMA \
        -N 20 \
        -p 19877 \
        -n best \
        -s mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip \
        -T 12

    # analyze boostrap data sets using the autoMRE function of RAxML
    raxmlHPC-PTHREADS-SSE3 \
        -m GTRGAMMA \
        -N autoMRE \
        -p 19877 \
        -b 7175 \
        -n bootreps \
        -s mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip \
        -T 12

.. warning:: Note that I am using 12 physical CPU cores here (`-T 12`).  You
    need to use the number of physical cores available on *your* machine.

Once those are finished running, we need to reconcile the "best" tree with the
ML bootstrap replicates:

.. code-block:: bash

    # reconcile the "best" ML tree with the bootreps
    raxmlHPC-SSE3 \
        -m GTRGAMMA \
        -f b \
        -t RAxML_bestTree.best \
        -z RAxML_bootstrap.bootreps

ExaBayes
--------

I have not included ExaBayes_ in the conda package for phyluce_ because it
generally requires MPI and should really be run on an HPC or a large local
machine.  If you would like to run it on your machine, you'll need to get it
installed on your own.

To run data in ExaBayes_, you need the PHYLIP file you just created as well as
2 other files in order to the the program - one giving partition information and
another giving configuration options.  Create these.

aln.part
^^^^^^^^

.. code-block:: bash

    DNA, p1=1-499929

config.nexus
^^^^^^^^^^^^

.. code-block:: bash

    #NEXUS

    begin run;
     numruns 4
     numgen 1e6
    end;

Run ExaBayes
^^^^^^^^^^^^

Now, we can run ExaBayes:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-exabayes

    # get a random number
    echo $RANDOM

    # this outputs
    5341

    # ensure the files we created above are in this directory, then
    mpirun -np 12 exabayes -f mafft-nexus-internal-trimmed-gblocks-clean-75p.phylip -n run1 -q aln.part -s 5341 -c config.nexus -R 4

.. warning:: Note that I am using 12 physical CPU cores here (`-T 12`).  You
    need to use the number of physical cores available on *your* machine.

Once those are finished running, we can summarize the posterior using the
`consense` and `postProcParams` program:

.. code-block:: bash

    consense -f ExaBayes_topologies.run1.phylip.* -n some_descriptive_name_here
    postProcParam -f ExaBayes_parameters.run1.phylip.* -n some_descriptive_name_here
