.. include:: ../global.rst

.. _TutorialOne:

*****************************
Tutorial I: UCE Phylogenomics
*****************************

In the following example, we are going to process raw read data from UCE enrichments performed against several divergent taxa so that you can get a feel for how a typical analysis goes.  For more general analysis notes, see the :ref:`UCE Processing` chapter.  That said, this is a good place to start.

The taxa we are working with will be:

* Mus musculus (PE100)
* Anolis carolinensis (PE100)
* Alligator mississippiensis (PE150)
* Gallus gallus (PE250)

Download the data
=================

You can download the data from figshare (http://dx.doi.org/10.6084/m9.figshare.1284521).  If you want to use the command line, you can use something like:

.. code-block:: bash

    # create a project directory
    mkdir uce-tutorial

    # change to that directory
    cd uce-tutorial

    # download the data into a file names fastq.zip
    wget -O fastq.zip https://ndownloader.figshare.com/articles/1284521/versions/2

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

    -rw-r--r--. 1 bcf users 4.4M Feb 22 14:14 Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf users 4.3M Feb 22 14:14 Alligator_mississippiensis_GGAGCTATGG_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf users 4.9M Feb 22 14:14 Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf users 4.9M Feb 22 14:15 Anolis_carolinensis_GGCGAAGGTT_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf users 7.6M Feb 22 14:15 Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf users 8.4M Feb 22 14:15 Gallus_gallus_TTCTCCTTCA_L001_R2_001.fastq.gz
    -rw-r--r--. 1 bcf users 4.9M Feb 22 14:16 Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
    -rw-r--r--. 1 bcf users 4.9M Feb 22 14:16 Mus_musculus_CTACAACGGC_L001_R2_001.fastq.gz

Alternatively, if you think of the filesystem as a tree-like structure, the directory in which we are working (`uce-tutorial`) would look like:

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


If you do not want to use the command line, you can download the data using the figshare interface or by clicking:

http://downloads.figshare.com/article/public/1284521

Count the read data
===================

Usually, we want a count of the actual number of reads in a given sequence file for a given species. We can do this several ways, but here, we'll use tools from  unix, because they are fast. The next line of code will count the lines in each R1 file (which should be equal to the reads in the R2 file) and divide that number by 4 to get the number of sequence reads.

.. code-block:: bash

    for i in *_R1_*.fastq.gz; do echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done

You should see:

.. code-block:: bash

    Alligator_mississippiensis_GGAGCTATGG_L001_R1_001.fastq.gz
    50000
    Anolis_carolinensis_GGCGAAGGTT_L001_R1_001.fastq.gz
    50000
    Gallus_gallus_TTCTCCTTCA_L001_R1_001.fastq.gz
    50000
    Mus_musculus_CTACAACGGC_L001_R1_001.fastq.gz
    50000

Notice that all the read counts are equal - that is because these 50,000 reads in each R1 and R2 file were subsampled, randomly, from a file of many more reads.

Clean the read data
===================

The data you just downloaded are actual, raw, untrimmed fastq data.  This means they contain adapter contamination and low quality bases.  We need to remove these - which you can do several ways.  We'll use another program that I wrote (illumiprocessor_) because it allows us to trim many different indexed adapters from individual-specific fastq files - something that is a pain to do by hand.  That said, you can certainly trim your reads however you would like. See the illumiprocessor_ website for instructions on installing the program.
 
To use this program, we will create a configuration file that we will use to inform the program about which adapters are in which READ1 and READ2 files.  The data we are trimming, here, are from TruSeq v3 libraries, but the indexes are 10 nucleotides long.  We will set up the trimming file with these parameters, but please see the illumiprocessor_ documentation for other options.

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

Now I run illumiprocessor_ against the data.  Note that I am using **4 physical CPU cores** to do this work.  You need to use the number of physical cores available on your machine, although there is not sense in using more cores than you have taxa (in this case).

.. code-block:: bash

    # go to the directory containing our config file and data
    cd uce-tutorial

    # run illumiprocessor

    illumiprocessor \
        --input raw-fastq/ \
        --output clean-fastq \
        --config illumiprocessor.conf \
        --cores 4

The output should look like the following:

.. code-block:: bash

    2021-02-22 14:59:26,488 - illumiprocessor - INFO - ==================== Starting illumiprocessor ===================
    2021-02-22 14:59:26,489 - illumiprocessor - INFO - Version: 2.0.9
    2021-02-22 14:59:26,489 - illumiprocessor - INFO - Argument --config: illumiprocessor.conf
    2021-02-22 14:59:26,489 - illumiprocessor - INFO - Argument --cores: 4
    2021-02-22 14:59:26,489 - illumiprocessor - INFO - Argument --input: /scratch/bfaircloth/uce-tutorial/raw-fastq
    2021-02-22 14:59:26,490 - illumiprocessor - INFO - Argument --log_path: None
    2021-02-22 14:59:26,490 - illumiprocessor - INFO - Argument --min_len: 40
    2021-02-22 14:59:26,490 - illumiprocessor - INFO - Argument --no_merge: False
    2021-02-22 14:59:26,490 - illumiprocessor - INFO - Argument --output: /scratch/bfaircloth/uce-tutorial/clean-fastq
    2021-02-22 14:59:26,490 - illumiprocessor - INFO - Argument --phred: phred33
    2021-02-22 14:59:26,491 - illumiprocessor - INFO - Argument --r1_pattern: None
    2021-02-22 14:59:26,491 - illumiprocessor - INFO - Argument --r2_pattern: None
    2021-02-22 14:59:26,491 - illumiprocessor - INFO - Argument --se: False
    2021-02-22 14:59:26,491 - illumiprocessor - INFO - Argument --trimmomatic: /home/bcf/conda/envs/phyluce/bin/trimmomatic
    2021-02-22 14:59:26,491 - illumiprocessor - INFO - Argument --verbosity: INFO
    2021-02-22 14:59:26,904 - illumiprocessor - INFO - Trimming samples with Trimmomatic
    Running....
    2021-02-22 14:59:36,754 - illumiprocessor - INFO - =================== Completed illumiprocessor ===================

Notice that the program has created a ``log`` file showing what it did, and it has also created a new directory holding the clean data that has the name ``clean-fastq`` (what you told it to name the directory). Within that new directory, there are taxon-specific folder for the cleaned reads. More specifically, your directory structure should look similar to the following (I've collapsed the list of raw-reads):

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

And, within each of those directories nested within the species-specific directory, there are additional files or links to files:

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

The `->` in the `raw-reads` directory above means there are symlinks_ to the files. I have removed the file paths and replaced them with `<PATH>` so that the figure will fit on a page.

The really important information is in the `split-adapter-quality-trimmed` directory - which now holds our reads that have had adapter-contamination and low-quality bases removed. Within this `split-adapter-quality-trimmed` directory,  the `READ1` and `READ2` files hold reads that remain in a pair (the reads are in the same consecutive order in each file).  The `READ-singleton` file holds READ1 reads **OR** READ2 reads that lost their "mate" or "paired-read" because of trimming or removal.

Quality control
---------------

You might want to get some idea of what effect the trimming has on read counts and overall read lengths. There are certainly other (better) tools out there to do this (like FastQC_), but you can get a reasonable idea of how good your reads are by running the following, which will output a CSV listing of read stats by sample:

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

    All files in dir with alligator_mississippiensis-READ1.fastq.gz,93699,8418476,89.84595353205476,0.059508742244529164,40,100,100.0
    All files in dir with anolis_carolinensis-READ-singleton.fastq.gz,92184,8401336,91.13659637247244,0.048890234925557836,40,100,100.0
    All files in dir with gallus_gallus-READ1.fastq.gz,99444,21218771,213.37406982824504,0.16122899415574637,40,251,250.0
    All files in dir with mus_musculus-READ2.fastq.gz,89841,8165734,90.89095179261139,0.052266485638855914,40,100,100.

Now, we're ready to assemble our reads.

.. _tutorial-assembly:

Assemble the data
=================

phyluce_ has several options for assembly - you can use velvet_, abyss_, or spades_.  For this tutorial, we are going to use spades_, because it seems to works best for most purposes, it is easy to install and run, and it works consistently.  The helper programs for the other assemblers use the same config file, so you can easily experiment with all of the assemblers.

To run an assembly, we need to create a another configuration file.  The assembly configuration file looks like the following, assuming we want to assemble all of our data from the organisms above:

.. code-block:: bash

    [samples]
    alligator_mississippiensis:/path/to/the/uce-tutorial/clean-fastq/alligator_mississippiensis/split-adapter-quality-trimmed/
    anolis_carolinensis:/path/to/the/uce-tutorial/clean-fastq/anolis_carolinensis/split-adapter-quality-trimmed/
    gallus_gallus:/path/to/the/uce-tutorial/clean-fastq/gallus_gallus/split-adapter-quality-trimmed/
    mus_musculus:/path/to/the/uce-tutorial/clean-fastq/mus_musculus/split-adapter-quality-trimmed/

You need to modify this file to use the path to the clean read data **on your computer** (``/path/to/the/`` is a placeholder, here). You will save this into a file named ``assembly.conf`` at the top of our ``uce-tutorial`` directory:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    └── trinity-assemblies

If you want to change the names on the left hand side of the colon in the config file, you can do so, but the paths on the right hand side need to point to our "clean" UCE raw reads.  If you have files in multiple locations, you can use any number of different paths on the right-hand side of the colon.

.. attention:: Although you can easily input new PATHs in this file, the
    **structure** of the data below the PATH you use must be the same - meaning
    that the structure and naming scheme for READ1, READ2, and READ-singleton
    must be the same.  Or, put another way, the assembly programs for phyluce
    assume the data **always** look like the following:

    .. code-block:: bash

        <some working folder name>
        └── clean-fastq
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

Now that we have that file created, copy it to our working directory, and run the ``phyluce_assembly_assemblo_spades`` program:

.. code-block:: bash

    # make sure we are at the top-level of our uce tutorial directory
    cd uce-tutorial

    # run the assembly
    phyluce_assembly_assemblo_spades \
        --conf assembly.conf \
        --output spades-assemblies \
        --cores 12

.. warning:: Note that I am using 12 physical CPU cores to do this work.  You
    need to use the number of physical cores available on **your** machine. ``phyluce_assembly_assemblo_spades`` assumes you have **at least** 8 GB of RAM 
    on your system, and it is better to have much more.  If you use more 
    CPU cores than you have or you specify more RAM than you have, the job 
    can fail.

    You can adjust the RAM dedicated to the job using the ``--memory`` option,
    which takes an integer value (in GB RAM).

.. note:: If you are wondering, Trinity_ is no longer supported in phyluce_

As the assembly proceeds, you should see output similar to the following:

.. code-block:: bash

    2021-02-26 21:12:16,615 - phyluce_assembly_assemblo_spades - INFO - =========== Starting phyluce_assembly_assemblo_spades ===========
    2021-02-26 21:12:16,615 - phyluce_assembly_assemblo_spades - INFO - Version: 1.7.0
    2021-02-26 21:12:16,615 - phyluce_assembly_assemblo_spades - INFO - Commit: None
    2021-02-26 21:12:16,615 - phyluce_assembly_assemblo_spades - INFO - Argument --config: /data/assembly.conf
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --cores: 12
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --dir: None
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --do_not_clean: False
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --log_path: None
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --memory: 8
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --output: /data/spades-assemblies
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --subfolder:
    2021-02-26 21:12:16,616 - phyluce_assembly_assemblo_spades - INFO - Argument --verbosity: INFO
    2021-02-26 21:12:16,617 - phyluce_assembly_assemblo_spades - INFO - Getting input filenames and creating output directories
    2021-02-26 21:12:16,765 - phyluce_assembly_assemblo_spades - INFO - ------------- Processing alligator_mississippiensis -------------
    2021-02-26 21:12:16,775 - phyluce_assembly_assemblo_spades - INFO - Finding fastq/fasta files
    2021-02-26 21:12:16,787 - phyluce_assembly_assemblo_spades - INFO - File type is fastq
    2021-02-26 21:12:16,787 - phyluce_assembly_assemblo_spades - INFO - Running SPAdes for PE data
    2021-02-26 21:13:35,643 - phyluce_assembly_assemblo_spades - INFO - Symlinking assembled contigs into /data/spades-assemblies/contigs
    ...[continued]...
    2021-02-26 21:19:06,618 - phyluce_assembly_assemblo_spades - INFO - Symlinking assembled contigs into /data/spades-assemblies/contigs
    2021-02-26 21:19:06,624 - phyluce_assembly_assemblo_spades - INFO - =========== Completed phyluce_assembly_assemblo_spades ==========

One the assembly is finished, have a look at the directory structure:

.. code-block:: bash

    uce-tutorial
    ├── assembly.conf
    ├── clean-fastq
    ├── illumiprocessor.conf
    ├── illumiprocessor.log
    ├── phyluce_assembly_assemblo_trinity.log
    ├── raw-fastq
    └── spades-assemblies
        ├── alligator_mississippiensis_trinity
        │   ├── contigs.fasta
        │   ├── scaffolds.fasta
        │   └── spades.log
        ├── anolis_carolinensis_trinity
        │   ├── contigs.fasta
        │   ├── scaffolds.fasta
        │   └── spades.log
        ├── contigs
        │   ├── alligator_mississippiensis.contigs.fasta -> ../alligator_mississippiensis_trinity/contigs.fasta
        │   ├── anolis_carolinensis.contigs.fasta -> ../anolis_carolinensis_trinity/contigs.fasta
        │   ├── gallus_gallus.contigs.fasta -> ../gallus_gallus_trinity/contigs.fasta
        │   └── mus_musculus.contigs.fasta -> ../mus_musculus_trinity/contigs.fasta
        ├── gallus_gallus_trinity
        │   ├── contigs.fasta
        │   ├── scaffolds.fasta
        │   └── spades.log
        └── mus_musculus_trinity
            ├── contigs.fasta
            ├── scaffolds.fasta
            └── spades.log

Your species-specific assembly files are in the `spades-assemblies` directory
nested within species-specific directories that correspond to the name you used
in the `assembly.conf` file (to the left of the colon).

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

    for i in spades-assemblies/contigs/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done

This should output something similar to the following.  I've added the header as
a comment:

.. code-block:: bash

    # samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
    alligator_mississippiensis.contigs.fasta,870,220436,253.37471264367815,6.9967833874072225,56,3831,236.0,5
    anolis_carolinensis.contigs.fasta,1136,355837,313.2367957746479,8.00195622090676,56,3182,243.0,9
    gallus_gallus.contigs.fasta,5228,2841631,543.5407421576128,2.8905564028218618,56,4117,495.5,107
    mus_musculus.contigs.fasta,1395,324498,232.61505376344087,6.895305881534584,56,1646,88.0,13

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
        --contigs spades-assemblies/contigs \
        --probes uce-5k-probes.fasta \
        --output uce-search-results

You should see output similar to the following (also stored in
`phyluce_assembly_assemblo_trinity.log`):

.. code-block:: bash

    2021-02-26 21:37:43,108 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Starting phyluce_assembly_match_contigs_to_probes =======
    2021-02-26 21:37:43,108 - phyluce_assembly_match_contigs_to_probes - INFO - Version: 1.7.0
    2021-02-26 21:37:43,109 - phyluce_assembly_match_contigs_to_probes - INFO - Commit: None
    2021-02-26 21:37:43,109 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --contigs: /data/spades-assemblies/contigs
    2021-02-26 21:37:43,109 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --csv: test-output.csv
    2021-02-26 21:37:43,109 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --dupefile: None
    2021-02-26 21:37:43,110 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --keep_duplicates: None
    2021-02-26 21:37:43,110 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --log_path: None
    2021-02-26 21:37:43,110 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_coverage: 80
    2021-02-26 21:37:43,110 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --min_identity: 80
    2021-02-26 21:37:43,110 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --output: /data/uce-search-results
    2021-02-26 21:37:43,111 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --probes: /data/uce-5k-probes.fasta
    2021-02-26 21:37:43,111 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --regex: ^(uce-\d+)(?:_p\d+.*)
    2021-02-26 21:37:43,111 - phyluce_assembly_match_contigs_to_probes - INFO - Argument --verbosity: INFO
    2021-02-26 21:37:43,225 - phyluce_assembly_match_contigs_to_probes - INFO - Creating the UCE-match database
    2021-02-26 21:37:43,266 - phyluce_assembly_match_contigs_to_probes - INFO - Processing contig data
    2021-02-26 21:37:43,266 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2021-02-26 21:37:44,680 - phyluce_assembly_match_contigs_to_probes - INFO - alligator_mississippiensis: 422 (48.51%) uniques of 870 contigs, 0 dupe probe matches, 3 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    2021-02-26 21:37:45,940 - phyluce_assembly_match_contigs_to_probes - INFO - anolis_carolinensis: 399 (35.12%) uniques of 1136 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    2021-02-26 21:37:53,686 - phyluce_assembly_match_contigs_to_probes - INFO - gallus_gallus: 3492 (66.79%) uniques of 5228 contigs, 0 dupe probe matches, 22 UCE loci removed for matching multiple contigs, 37 contigs removed for matching multiple UCE loci
    2021-02-26 21:37:54,853 - phyluce_assembly_match_contigs_to_probes - INFO - mus_musculus: 324 (23.23%) uniques of 1395 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci
    2021-02-26 21:37:54,853 - phyluce_assembly_match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2021-02-26 21:37:54,853 - phyluce_assembly_match_contigs_to_probes - INFO - The LASTZ alignments are in /data/uce-search-results
    2021-02-26 21:37:54,854 - phyluce_assembly_match_contigs_to_probes - INFO - The UCE match database is in /data/uce-search-results/probe.matches.sqlite
    2021-02-26 21:37:54,854 - phyluce_assembly_match_contigs_to_probes - INFO - ======= Completed phyluce_assembly_match_contigs_to_probes ======

The header info at the top tells us exactly what version of the code we are
running and keeps track of our options.  The important output is:

.. code-block:: bash

    alligator_mississippiensis: 422 (48.51%) uniques of 870 contigs, 0 dupe probe matches, 3 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    anolis_carolinensis: 399 (35.12%) uniques of 1136 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 2 contigs removed for matching multiple UCE loci
    gallus_gallus: 3492 (66.79%) uniques of 5228 contigs, 0 dupe probe matches, 22 UCE loci removed for matching multiple contigs, 37 contigs removed for matching multiple UCE loci
    mus_musculus: 324 (23.23%) uniques of 1395 contigs, 0 dupe probe matches, 0 UCE loci removed for matching multiple contigs, 1 contigs removed for matching multiple UCE loci

Which we can break down to the following (for `alligator_mississippiensis`):

.. code-block:: bash

    alligator_mississippiensis:
        422 (48.51%) uniques of 870 contigs
        0 dupe probe matches
        3 UCE loci removed for matching multiple contigs
        2 contigs removed for matching multiple UCE loci

These are the capture data for the `alligator_mississippiensis` sample.  We targeted 5k UCE loci in this sample and recovered roughly 422 of those loci **in this subsampled set of reads**. Before reaching that total of 422 loci, we removed 3 UCE loci and 2 contigs from the data set because they looked like duplicates (probes supposedly targeting different loci hit the same contig or two supposedly different contigs hit probes designed for a single UCE locus).

.. admonition:: Question: Why is the count of UCE loci different by sample?
    :class: admonition tip

    For these example data, we enriched some samples (alligator_mississippiensis
    and gallus_gallus) for 5k UCE loci, while we enriched others
    (anolis_carolinensis and mus_musculus) for 2.5k UCE loci.  Additionally, the
    2.5k UCE enrichments did not work very well (operator error).  Finally, because
    we have subsampled the data to make the files of reasonable size for the 
    tutorial, that process removes lots of read that **would** have assembled 
    into UCE contigs (e.g., if we use ALL the data, we recover 4000+ UCE contigs
    for alligator).

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

.. _UceExtraction:

Extracting UCE loci
===================

Now that we have located UCE loci, we need to determine which taxa we want in
our analysis, create a list of those taxa, and then generate a list of which
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

    2021-03-01 15:46:14,277 - phyluce_assembly_get_match_counts - INFO - =========== Starting phyluce_assembly_get_match_counts ==========
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Version: 1.7.0
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Commit: None
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --extend_locus_db: None
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --incomplete_matrix: True
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --keep_counts: False
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --locus_db: /data/uce-search-results/probe.matches.sqlite
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --log_path: None
    2021-03-01 15:46:14,278 - phyluce_assembly_get_match_counts - INFO - Argument --optimize: False
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --output: /data/taxon-sets/all/all-taxa-incomplete.conf
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --random: False
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --sample_size: 10
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --samples: 10
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --silent: False
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_group: all
    2021-03-01 15:46:14,279 - phyluce_assembly_get_match_counts - INFO - Argument --taxon_list_config: /data/taxon-set.conf
    2021-03-01 15:46:14,280 - phyluce_assembly_get_match_counts - INFO - Argument --verbosity: INFO
    2021-03-01 15:46:14,319 - phyluce_assembly_get_match_counts - INFO - There are 4 taxa in the taxon-group '[all]' in the config file taxon-set.conf
    2021-03-01 15:46:14,319 - phyluce_assembly_get_match_counts - INFO - Getting UCE names from database
    2021-03-01 15:46:44,488 - phyluce_assembly_get_match_counts - INFO - There are 5041 total UCE loci in the database
    2021-03-01 15:46:44,567 - phyluce_assembly_get_match_counts - INFO - Getting UCE matches by organism to generate a INCOMPLETE matrix
    2021-03-01 15:46:44,569 - phyluce_assembly_get_match_counts - INFO - There are 3653 UCE loci in an INCOMPLETE matrix
    2021-03-01 15:46:44,569 - phyluce_assembly_get_match_counts - INFO - Writing the taxa and loci in the data matrix to /data/taxon-sets/all/all-taxa-incomplete.conf
    2021-03-01 15:46:44,574 - phyluce_assembly_get_match_counts - INFO - ========== Completed phyluce_assembly_get_match_counts ==========

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
        --contigs ../../spades-assemblies/contigs \
        --locus-db ../../uce-search-results/probe.matches.sqlite \
        --match-count-output all-taxa-incomplete.conf \
        --output all-taxa-incomplete.fasta \
        --incomplete-matrix all-taxa-incomplete.incomplete \
        --log-path log

The output should look something like the following:

.. code-block:: bash

    2021-03-01 15:47:55,574 - phyluce_assembly_get_fastas_from_match_counts - INFO - ===== Starting phyluce_assembly_get_fastas_from_match_counts ====
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Version: 1.7.0
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Commit: None
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --contigs: /data/spades-assemblies/contigs
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --extend_locus_contigs: None
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --extend_locus_db: None
    2021-03-01 15:47:55,575 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --incomplete_matrix: /data/taxon-sets/all/all-taxa-incomplete.incomplete
    2021-03-01 15:47:55,576 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --locus_db: /data/uce-search-results/probe.matches.sqlite
    2021-03-01 15:47:55,576 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 15:47:55,576 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --match_count_output: /data/taxon-sets/all/all-taxa-incomplete.conf
    2021-03-01 15:47:55,576 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --output: /data/taxon-sets/all/all-taxa-incomplete.fasta
    2021-03-01 15:47:55,576 - phyluce_assembly_get_fastas_from_match_counts - INFO - Argument --verbosity: INFO
    2021-03-01 15:47:55,609 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 4 taxa in the match-count-config file named all-taxa-incomplete.conf
    2021-03-01 15:47:55,612 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 3653 UCE loci in an INCOMPLETE matrix
    2021-03-01 15:47:55,614 - phyluce_assembly_get_fastas_from_match_counts - INFO - ---------Getting UCE loci for alligator_mississippiensis---------
    2021-03-01 15:47:55,925 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 422 UCE loci for alligator_mississippiensis
    2021-03-01 15:47:55,925 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for alligator_mississippiensis
    2021-03-01 15:47:56,396 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /data/taxon-sets/all/all-taxa-incomplete.incomplete
    2021-03-01 15:47:56,399 - phyluce_assembly_get_fastas_from_match_counts - INFO - -------------Getting UCE loci for anolis_carolinensis------------
    2021-03-01 15:47:56,638 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 399 UCE loci for anolis_carolinensis
    2021-03-01 15:47:56,638 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for anolis_carolinensis
    2021-03-01 15:47:57,077 - phyluce_assembly_get_fastas_from_match_counts - INFO - Replaced <20 ambiguous bases (N) in 7 contigs for anolis_carolinensis
    2021-03-01 15:47:57,077 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /data/taxon-sets/all/all-taxa-incomplete.incomplete
    2021-03-01 15:47:57,079 - phyluce_assembly_get_fastas_from_match_counts - INFO - ----------------Getting UCE loci for gallus_gallus---------------
    2021-03-01 15:47:58,940 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 3492 UCE loci for gallus_gallus
    2021-03-01 15:47:58,940 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for gallus_gallus
    2021-03-01 15:48:01,965 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /data/taxon-sets/all/all-taxa-incomplete.incomplete
    2021-03-01 15:48:01,966 - phyluce_assembly_get_fastas_from_match_counts - INFO - ----------------Getting UCE loci for mus_musculus----------------
    2021-03-01 15:48:02,168 - phyluce_assembly_get_fastas_from_match_counts - INFO - There are 324 UCE loci for mus_musculus
    2021-03-01 15:48:02,168 - phyluce_assembly_get_fastas_from_match_counts - INFO - Parsing and renaming contigs for mus_musculus
    2021-03-01 15:48:02,508 - phyluce_assembly_get_fastas_from_match_counts - INFO - Writing missing locus information to /data/taxon-sets/all/all-taxa-incomplete.incomplete
    2021-03-01 15:48:02,519 - phyluce_assembly_get_fastas_from_match_counts - INFO - ==== Completed phyluce_assembly_get_fastas_from_match_counts ====

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

.. _ExplodingFasta:

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
        --output exploded-fastas \
        --by-taxon

    # get summary stats on the FASTAS
    for i in exploded-fastas/*.fasta;
    do
        phyluce_assembly_get_fasta_lengths --input $i --csv;
    done

    # samples,contigs,total bp,mean length,95 CI length,min length,max length,median legnth,contigs >1kb
    alligator-mississippiensis.unaligned.fasta,422,118864,281.66824644549763,9.03914154267783,206,3831,250.5,1
    anolis-carolinensis.unaligned.fasta,399,206926,518.6115288220551,11.261692813904311,206,1115,525.0,1
    gallus-gallus.unaligned.fasta,3492,2157267,617.7740549828179,3.075870593111506,307,1503,607.0,81
    mus-musculus.unaligned.fasta,324,188082,580.5,13.930004844688803,207,1058,623.5,6

Aligning UCE loci
=================

You have lots of options when aligning UCE loci.  You can align the loci and use
those alignments with no trimming, you can edge-trim the alignments following
some algorithm, and you can end+internally trim alignments following some
algorithm. It's hard to say what is best in all situations.  When taxa are
"closely" related (< 30-50 MYA, perhaps), I think that edge-trimming alignments
is reasonable.  When the taxa you are interested in span a wider range of
divergence times (> 50 MYA), you may want to think about internal trimming.

How you accomplish you edge- or internal-trimming is also a decision you need to
make. In phyluce_, we implement our edge-trimming algorithm by running the
alignment program "as-is" (i.e., without the `--no-trim`) option.  We do
internal-trimming by turning off trimming using `--no-trim`, then passing the
resulting alignments (in FASTA format) to a parallel wrapper around Gblocks_.

You also have a choice of aligner - mafft_ or muscle_ (or you can externally
align UCE loci using a tool like SATé, as well).

Generally, I would use mafft_.


.. _EdgeTrimming:

Edge trimming
-------------

Edge trimming your alignments is a relatively simple matter.  You can run edge
trimming, as follows:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # align the data
    phyluce_align_seqcap_align \
        --input all-taxa-incomplete.fasta \
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

    2021-03-01 15:51:21,854 - phyluce_align_seqcap_align - INFO - ============== Starting phyluce_align_seqcap_align ==============
    2021-03-01 15:51:21,855 - phyluce_align_seqcap_align - INFO - Version: 1.7.0
    2021-03-01 15:51:21,855 - phyluce_align_seqcap_align - INFO - Commit: None
    2021-03-01 15:51:21,856 - phyluce_align_seqcap_align - INFO - Argument --aligner: mafft
    2021-03-01 15:51:21,856 - phyluce_align_seqcap_align - INFO - Argument --ambiguous: False
    2021-03-01 15:51:21,856 - phyluce_align_seqcap_align - INFO - Argument --cores: 12
    2021-03-01 15:51:21,856 - phyluce_align_seqcap_align - INFO - Argument --input: /data/taxon-sets/all/all-taxa-incomplete.fasta
    2021-03-01 15:51:21,857 - phyluce_align_seqcap_align - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 15:51:21,857 - phyluce_align_seqcap_align - INFO - Argument --max_divergence: 0.2
    2021-03-01 15:51:21,857 - phyluce_align_seqcap_align - INFO - Argument --min_length: 100
    2021-03-01 15:51:21,857 - phyluce_align_seqcap_align - INFO - Argument --no_trim: False
    2021-03-01 15:51:21,857 - phyluce_align_seqcap_align - INFO - Argument --notstrict: True
    2021-03-01 15:51:21,858 - phyluce_align_seqcap_align - INFO - Argument --output: /data/taxon-sets/all/mafft-nexus-edge-trimmed
    2021-03-01 15:51:21,858 - phyluce_align_seqcap_align - INFO - Argument --output_format: nexus
    2021-03-01 15:51:21,858 - phyluce_align_seqcap_align - INFO - Argument --proportion: 0.65
    2021-03-01 15:51:21,858 - phyluce_align_seqcap_align - INFO - Argument --taxa: 4
    2021-03-01 15:51:21,858 - phyluce_align_seqcap_align - INFO - Argument --threshold: 0.65
    2021-03-01 15:51:21,859 - phyluce_align_seqcap_align - INFO - Argument --verbosity: INFO
    2021-03-01 15:51:21,859 - phyluce_align_seqcap_align - INFO - Argument --window: 20
    2021-03-01 15:51:21,859 - phyluce_align_seqcap_align - INFO - Building the locus dictionary
    2021-03-01 15:51:21,859 - phyluce_align_seqcap_align - INFO - Removing ALL sequences with ambiguous bases...
    2021-03-01 15:51:22,328 - phyluce_align_seqcap_align - WARNING - DROPPED locus uce-4698. Too few taxa (N < 3).
    [many more loci dropped here]
    2021-03-01 15:51:22,750 - phyluce_align_seqcap_align - INFO - Aligning with MAFFT
    2021-03-01 15:51:22,751 - phyluce_align_seqcap_align - INFO - Alignment begins. 'X' indicates dropped alignments (these are reported after alignment)
    .................[continued]
    2021-03-01 15:51:28,544 - phyluce_align_seqcap_align - INFO - Alignment ends
    2021-03-01 15:51:28,545 - phyluce_align_seqcap_align - INFO - Writing output files
    2021-03-01 15:51:28,788 - phyluce_align_seqcap_align - INFO - ============== Completed phyluce_align_seqcap_align =============

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

    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - Version: 1.7.0
    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - Commit: None
    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-edge-trimmed
    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
    2021-03-01 15:58:43,241 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
    2021-03-01 15:58:43,242 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 15:58:43,242 - phyluce_align_get_align_summary_data - INFO - Argument --output: None
    2021-03-01 15:58:43,242 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2021-03-01 15:58:43,242 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
    2021-03-01 15:58:43,242 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
    2021-03-01 15:58:43,251 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2021-03-01 15:58:43,596 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2021-03-01 15:58:43,597 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:	190
    2021-03-01 15:58:43,597 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:	86,032
    2021-03-01 15:58:43,597 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:	452.80
    2021-03-01 15:58:43,597 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:	21.51
    2021-03-01 15:58:43,597 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:	126
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:	907
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - ------------------- Informative Sites summary -------------------
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Sites] loci:	190
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Sites] total:	58
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Sites] mean:	0.31
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Sites] 95% CI:	0.18
    2021-03-01 15:58:43,598 - phyluce_align_get_align_summary_data - INFO - [Sites] min:	0
    2021-03-01 15:58:43,599 - phyluce_align_get_align_summary_data - INFO - [Sites] max:	12
    2021-03-01 15:58:43,600 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2021-03-01 15:58:43,600 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:	3.15
    2021-03-01 15:58:43,600 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:	0.05
    2021-03-01 15:58:43,600 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:		3
    2021-03-01 15:58:43,600 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:		4
    2021-03-01 15:58:43,601 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2021-03-01 15:58:43,601 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:	    9.36
    2021-03-01 15:58:43,601 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:	0.84
    2021-03-01 15:58:43,601 - phyluce_align_get_align_summary_data - INFO - [Missing] min:	    0.00
    2021-03-01 15:58:43,601 - phyluce_align_get_align_summary_data - INFO - [Missing] max:	    32.67
    2021-03-01 15:58:43,603 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2021-03-01 15:58:43,603 - phyluce_align_get_align_summary_data - INFO - [All characters]	268,288
    2021-03-01 15:58:43,603 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]		235,385
    2021-03-01 15:58:43,603 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]		190 alignments
    2021-03-01 15:58:43,604 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]		29 alignments
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]		29 alignments
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]		29 alignments
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]		29 alignments
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 7,819 times
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Characters] '?' is present 25,084 times
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 70,371 times
    2021-03-01 15:58:43,605 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 48,431 times
    2021-03-01 15:58:43,606 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 41,941 times
    2021-03-01 15:58:43,606 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 74,642 times
    2021-03-01 15:58:43,606 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========

.. attention:: Note that there are only 2 sets of counts in the ``Data matrix
    completeness`` section because (1) we dropped all loci having fewer than 3
    taxa and (2) that only leaves two remaining options.

The most important data here are the number of loci we have and the number of
loci in data matrices of different completeness. The locus length stats are
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
        --input all-taxa-incomplete.fasta \
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

    2021-03-01 16:01:38,320 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO -  Starting phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed
    2021-03-01 16:01:38,321 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Version: 1.7.0
    2021-03-01 16:01:38,321 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Commit: None
    2021-03-01 16:01:38,321 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-internal-trimmed
    2021-03-01 16:01:38,321 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b1: 0.5
    2021-03-01 16:01:38,321 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b2: 0.85
    2021-03-01 16:01:38,322 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b3: 8
    2021-03-01 16:01:38,322 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --b4: 10
    2021-03-01 16:01:38,322 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --cores: 12
    2021-03-01 16:01:38,322 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --input_format: fasta
    2021-03-01 16:01:38,323 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 16:01:38,323 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --output: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2021-03-01 16:01:38,323 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --output_format: nexus
    2021-03-01 16:01:38,323 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Argument --verbosity: INFO
    2021-03-01 16:01:38,323 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Getting aligned sequences for trimming
    2021-03-01 16:01:38,338 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Alignment trimming begins.
    .................[continued]
    2021-03-01 16:01:38,729 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Alignment trimming ends
    2021-03-01 16:01:38,730 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO - Writing output files
    2021-03-01 16:01:38,901 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - WARNING - Unable to write uce-816 - alignment too short
    2021-03-01 16:01:39,099 - phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed - INFO -  Completed phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed

The `.` values that you see represent loci that were aligned and succesfully
trimmed. Any `X` values that you see represent loci that were aligned and
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

    2021-03-01 16:03:34,322 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
    2021-03-01 16:03:34,322 - phyluce_align_get_align_summary_data - INFO - Version: 1.7.0
    2021-03-01 16:03:34,322 - phyluce_align_get_align_summary_data - INFO - Commit: None
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --output: None
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
    2021-03-01 16:03:34,323 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
    2021-03-01 16:03:34,334 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2021-03-01 16:03:34,641 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2021-03-01 16:03:34,642 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:	189
    2021-03-01 16:03:34,642 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:	70,380
    2021-03-01 16:03:34,643 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:	372.38
    2021-03-01 16:03:34,643 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:	21.61
    2021-03-01 16:03:34,643 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:	140
    2021-03-01 16:03:34,643 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:	797
    2021-03-01 16:03:34,644 - phyluce_align_get_align_summary_data - INFO - ------------------- Informative Sites summary -------------------
    2021-03-01 16:03:34,645 - phyluce_align_get_align_summary_data - INFO - [Sites] loci:	189
    2021-03-01 16:03:34,645 - phyluce_align_get_align_summary_data - INFO - [Sites] total:	69
    2021-03-01 16:03:34,645 - phyluce_align_get_align_summary_data - INFO - [Sites] mean:	0.37
    2021-03-01 16:03:34,645 - phyluce_align_get_align_summary_data - INFO - [Sites] 95% CI:	0.19
    2021-03-01 16:03:34,646 - phyluce_align_get_align_summary_data - INFO - [Sites] min:	0
    2021-03-01 16:03:34,646 - phyluce_align_get_align_summary_data - INFO - [Sites] max:	12
    2021-03-01 16:03:34,648 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2021-03-01 16:03:34,649 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:		3.15
    2021-03-01 16:03:34,649 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:	0.05
    2021-03-01 16:03:34,649 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:		3
    2021-03-01 16:03:34,650 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:		4
    2021-03-01 16:03:34,650 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2021-03-01 16:03:34,651 - phyluce_align_get_align_summary_data - INFO - [Missing] mean:	0.00
    2021-03-01 16:03:34,651 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:	0.00
    2021-03-01 16:03:34,651 - phyluce_align_get_align_summary_data - INFO - [Missing] min:	0.00
    2021-03-01 16:03:34,651 - phyluce_align_get_align_summary_data - INFO - [Missing] max:	0.00
    2021-03-01 16:03:34,655 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2021-03-01 16:03:34,655 - phyluce_align_get_align_summary_data - INFO - [All characters]	222,814
    2021-03-01 16:03:34,655 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]		215,091
    2021-03-01 16:03:34,656 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2021-03-01 16:03:34,656 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]		189 alignments
    2021-03-01 16:03:34,656 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]		189 alignments
    2021-03-01 16:03:34,657 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]		189 alignments
    2021-03-01 16:03:34,657 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]		189 alignments
    2021-03-01 16:03:34,657 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]		189 alignments
    2021-03-01 16:03:34,657 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]		189 alignments
    2021-03-01 16:03:34,658 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]		29 alignments
    2021-03-01 16:03:34,658 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]		29 alignments
    2021-03-01 16:03:34,658 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]		29 alignments
    2021-03-01 16:03:34,658 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]		29 alignments
    2021-03-01 16:03:34,659 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2021-03-01 16:03:34,659 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 7,723 times
    2021-03-01 16:03:34,659 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 64,167 times
    2021-03-01 16:03:34,659 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 44,488 times
    2021-03-01 16:03:34,660 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 37,901 times
    2021-03-01 16:03:34,660 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 68,535 times
    2021-03-01 16:03:34,660 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========

Alignment cleaning
==================

If you look in one of the files for the alignments we currently have, you will notice that each alignment contains a name that is a combination of the taxon name + the locus name for that taxon.  This is not what we want downstream, but it does enable us to ensure the correct data went into each alignment.  So, we need to clean our alignments. For the remainder of this tutorial, we will work with the Gblocks_ trimmed alignments, so we will clean those alignments:

 .. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # align the data - turn off trimming and output FASTA
    phyluce_align_remove_locus_name_from_files \
        --alignments mafft-nexus-internal-trimmed-gblocks \
        --output mafft-nexus-internal-trimmed-gblocks-clean \
        --cores 12 \
        --log-path log

The output should be similar to:

.. code-block:: bash

    2021-03-01 16:05:12,605 - phyluce_align_remove_locus_name_from_files - INFO - === Starting phyluce_align_remove_locus_name_from_files ===
    2021-03-01 16:05:12,605 - phyluce_align_remove_locus_name_from_files - INFO - Version: 1.7.0
    2021-03-01 16:05:12,605 - phyluce_align_remove_locus_name_from_files - INFO - Commit: None
    2021-03-01 16:05:12,606 - phyluce_align_remove_locus_name_from_files - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks
    2021-03-01 16:05:12,606 - phyluce_align_remove_locus_name_from_files - INFO - Argument --cores: 12
    2021-03-01 16:05:12,606 - phyluce_align_remove_locus_name_from_files - INFO - Argument --input_format: nexus
    2021-03-01 16:05:12,607 - phyluce_align_remove_locus_name_from_files - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 16:05:12,607 - phyluce_align_remove_locus_name_from_files - INFO - Argument --output: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean
    2021-03-01 16:05:12,607 - phyluce_align_remove_locus_name_from_files - INFO - Argument --output_format: nexus
    2021-03-01 16:05:12,608 - phyluce_align_remove_locus_name_from_files - INFO - Argument --taxa: None
    2021-03-01 16:05:12,608 - phyluce_align_remove_locus_name_from_files - INFO - Argument --verbosity: INFO
    2021-03-01 16:05:12,608 - phyluce_align_remove_locus_name_from_files - INFO - Getting alignment files
    Running...............[continued]
    2021-03-01 16:05:12,848 - phyluce_align_remove_locus_name_from_files - INFO - Taxon names in alignments: gallus_gallus,mus_musculus,anolis_carolinensis,alligator_mississippiensis
    2021-03-01 16:05:12,849 - phyluce_align_remove_locus_name_from_files - INFO - === Completed phyluce_align_remove_locus_name_from_files ==

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

Now, if you look in one of the individual alignment files, you will see that the locus names are removed.  We're ready to generate our final data matrices.

Final data matrices
===================

For the most part, I analyze 75% and 95% complete matrices, where "completeness"
for the 75% matrix means that, in a study of 100 taxa (total), all alignments
will contain at least 75 of these 100 taxa.  Similarly, for the 95% matrix, in a
study of 100 taxa, all alignments will contain 95 of these 100 taxa.

.. attention:: Notice that this metric for completeness does not pay attention
    to which taxa are in which alignments - so the 75%, above, does **not** mean
    that a given taxon will have data in all 75 of 100 alignments.

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

    2021-03-01 16:06:56,294 - phyluce_align_get_only_loci_with_min_taxa - INFO - ======= Starting phyluce_align_get_only_loci_with_min_taxa ======
    2021-03-01 16:06:56,294 - phyluce_align_get_only_loci_with_min_taxa - INFO - Version: 1.7.0
    2021-03-01 16:06:56,294 - phyluce_align_get_only_loci_with_min_taxa - INFO - Commit: None
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --cores: 12
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --input_format: nexus
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --output: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --percent: 0.75
    2021-03-01 16:06:56,295 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --taxa: 4
    2021-03-01 16:06:56,296 - phyluce_align_get_only_loci_with_min_taxa - INFO - Argument --verbosity: INFO
    2021-03-01 16:06:56,296 - phyluce_align_get_only_loci_with_min_taxa - INFO - Getting alignment files
    2021-03-01 16:06:56,534 - phyluce_align_get_only_loci_with_min_taxa - INFO - Copied 189 alignments of 189 total containing ≥ 0.75 proportion of taxa (n = 3)
    2021-03-01 16:06:56,534 - phyluce_align_get_only_loci_with_min_taxa - INFO - ====== Completed phyluce_align_get_only_loci_with_min_taxa ======

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

Preparing data for downstream analysis
======================================

Now that we have our `75p` data matrix completed, we can generate input files for subsequent phylogenetic analysis.  For the most part, I use RAxML or IQTree, both of which will take a phylip-formatted file as input. Formatting our `75p` data into a phylip file for these programs is rather easy.  To do that, run:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # build the concatenated data matrix
    phyluce_align_concatenate_alignments \
        --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
        --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
        --phylip \
        --log-path log

The output from this program will look like:

.. code-block:: bash

    2021-03-01 16:09:01,391 - phyluce_align_concatenate_alignments - INFO - ========= Starting phyluce_align_concatenate_alignments =========
    2021-03-01 16:09:01,391 - phyluce_align_concatenate_alignments - INFO - Version: 1.7.0
    2021-03-01 16:09:01,391 - phyluce_align_concatenate_alignments - INFO - Commit: None
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --alignments: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --input_format: nexus
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --log_path: /data/taxon-sets/all/log
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --nexus: False
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --output: /data/taxon-sets/all/mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --phylip: True
    2021-03-01 16:09:01,392 - phyluce_align_concatenate_alignments - INFO - Argument --verbosity: INFO
    2021-03-01 16:09:01,393 - phyluce_align_concatenate_alignments - INFO - Reading input alignments in NEXUS format
    2021-03-01 16:09:01,393 - phyluce_align_concatenate_alignments - INFO - Getting alignment files
    2021-03-01 16:09:01,733 - phyluce_align_concatenate_alignments - INFO - Concatenating files
    2021-03-01 16:09:01,796 - phyluce_align_concatenate_alignments - INFO - Writing concatenated alignment to PHYLIP format (with charsets)
    2021-03-01 16:09:01,803 - phyluce_align_concatenate_alignments - INFO - ========= Completed phyluce_align_concatenate_alignments ========

.. attention:: Charsets are now output by default for all data sets.
    You generally want these and the cost for them is low.

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

If you need to output concatenated data in ``NEXUS`` format, you can simply run the same program as above, but change ``--phylip`` to ``--nexus``, like:

.. code-block:: bash

    # make sure we are in the correct directory
    cd uce-tutorial/taxon-sets/all

    # build the concatenated data matrix
    phyluce_align_concatenate_alignments \
        --alignments mafft-nexus-internal-trimmed-gblocks-clean-75p \
        --output mafft-nexus-internal-trimmed-gblocks-clean-75p-raxml \
        --nexus \
        --log-path log

Downstream Analysis
-------------------

The above data are ready to analyze in a program like RAxML_ or IQTree_.  See the documentation for those programs to learn how to use them.  

If you are performing locus-based inference (e.g. so-called gene-tree, species-tree methods, you can analyze the individual locus alignments present in the folder you created just prior to concatenation (e.g. ``mafft-nexus-internal-trimmed-gblocks-clean-75p``).  If you need to convert those files into another format (by default, they are in ``NEXUS`` format), you can use ``phyluce_align_convert_one_align_to_another``.  For performing analyses to infer a locus-tree from individual alignments, you might look into a program like pargenes_.  You can also reasonably script IQTree_ to do this on, for example, an HPC system (we do it using GNU Parallel to send each alignment to one compute core, where we run IQTree_ on that alignment).

Next Steps
==========

After completing the tutorial, you should have a reasonably good idea of how to use phyluce_ in a day-to-day situation.  If you want to know more about specifics, you can read through the :ref:`Daily Use` sections, which provide additional detail.  Also be sure to poke around the other programs that come with phyluce_ - short list of which you can find in the :ref:`List of Programs`.
