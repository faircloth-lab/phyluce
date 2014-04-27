.. include:: global.rst

********************
Program organization
********************

phyluce_ is organized hierarchically.  There are *sub-commands* of the phyluce
command that cover different types of analyses you made need to run, and within
each of these *sub-commands* there are logically groups *sub-sub-commands*.
An example should make this clearer.

When you run the main program, you see output similar to::

    $ phyluce
    usage: phyluce [-h] [-V] command ...

    phyluce is a software package for processing UCEand other phylogenomic data
    for systematics andpopulation genetics.

    positional arguments:
      command
        assemble     Assemble cleaned/trimmed sequencing reads.
        fetch        Fetch (AKA get) UCE contigs from assembled data.
        align        Alignment routines for UCE (and other) FASTA data.
        convert      Convert files to different formats.
        misc         Miscellaneous utilities for files and directories.
        help         Get help info on a phyluce command.

    optional arguments:
      -h, --help     show this help message and exit
      -V, --version  show program's version number and exit

The commands that are listed are *sub-commands* of the main phyluce_ command.
In this example there are six of them covering 5 categories of analysis

* assembly
* UCE discovery
* alignment
* conversion
* miscellaneous

Sub-commands have sub-commands
==============================

When you run one of these subcommands you will see that each has its own set of
*sub-sub-commands*.  So, if we run `phyluce assembly`, we will see two assembly
options - assembly using velvet and assembly using abyss::

    $ phyluce assemble
    usage: phyluce assemble [-h] command ...

    Assemble cleaned/trimmed sequencing reads.

    positional arguments:
      command
        velvet    Assemble reads using velvet.
        abyss     Assemble reads using abYss.

    optional arguments:
      -h, --help  show this help message and exit

Sub-sub-commands have their own options
=======================================

Now, for each *sub-sub-command* there are additional options to run that
*sub-sub-command* against some data.  So, if you run::

    $ phyluce assemble velvet --help
    usage: phyluce assemble velvet [-h] (--config CONFIG | --dir DIR) --output
                                   OUTPUT [--kmer KMER] [--cores CORES]
                                   [--subfolder SUBFOLDER]
                                   [--verbosity {INFO,WARN,CRITICAL}]
                                   [--log-path LOG_PATH] [--clean]

    Assemble reads using velvet.

    optional arguments:
      -h, --help            show this help message and exit
      --config CONFIG       A configuration file containing reads to assemble
      --dir DIR             A directory of reads to assemble
      --output OUTPUT       The directory in which to store the assembly data
      --kmer KMER           The kmer value to use
      --cores CORES         The number of compute cores/threads to run with velvet
      --subfolder SUBFOLDER
                            A subdirectory, below the level of the group,
                            containing the reads
      --verbosity {INFO,WARN,CRITICAL}
                            The logging level to use
      --log-path LOG_PATH   The path to a directory to hold logs.
      --clean               Cleanup all intermediate velvet files

This will show you the options you need to use to run an assembly using velvet.


Number of sub-commands
=======================

Some sub-commands of phyluce_ have more sub-sub-commands than others.  For
example, if we run `phyluce align` we will see all of the sub-commands of the
*alignment* sub-command::

    $ phyluce align
    usage: phyluce align [-h] command ...

    Alignment routines for UCE (and other) FASTA data.

    positional arguments:
      command
        mafft     Align UCE FASTA data using MAFFT.
        muscle    Align UCE FASTA data using MUSCLE.
        clean     Clean a set of newly-aligned UCE data.
        check     Check alignments for problems.
        remove    Remove problematic alignments.
        repair    Repair problematic alignments (convert to N-bases).
        adjust    Adjust a set of alignment files by adding missing taxa.
        sift      Sift alignments and output those with more than --percent-taxa
                  or more than --min-taxa.
        pull      Pull taxon/taxa from a set of alignments. Output new alignments
                  without pulled taxon/taxa.
        pluck     Pluck taxon from a set of alignments. Output plucked taxon as
                  FASTA sequence file.
        explode   Explode a set of alignment files by-locus or by-taxon into the
                  consituent set of FASTA sequence files.
        stats     Compute alignment stats.
        sites     Compute the number of informative sites.

    optional arguments:
      -h, --help  show this help message and exit

Description of sub-commands
===========================

You can find a description of the sub-commands in the :ref:`Subcommands`
section.
