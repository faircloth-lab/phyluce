.. include:: global.rst

Software Installation Overview
==============================

phyluce_ (phy-loo-chee) is a software package for analyzing data collected from
ultraconserved elements in organismal genomes.  The package includes a number of
tools spanning assembly of raw read data to contigs (using velvet_) to
identification of UCE from assembled contigs to alignment generation, trimming,
and data preparation for analysis.

phyluce_ also has many parts, some necessary for certain operations, and others
useful only for very specific operations.  phyluce_ is also under constant
development, and the code-base changes rather rapidly, as we need new features,
generate new data formats, and fix bugs.

The code within this package assumes that you have a familiarity with the
command line (AKA Terminal in OSX) and **all** of the tools within the package
are driven from the command-line.  This allows us to focus on generating new
tools which are largely useable across different operating systems, without
getting bogged down with the intricacies of generating platform-specific
graphical user interfaces.

Below, we've outlined the necessary steps for installing the software and
pre-requisites.  We have focused on installation steps for OSX and Unix/Linux
because those are the operating systems that we use.

Installing the necessary pre-requisites is the most challenging part of getting
the various software components working.  On OSX, we've tried to simplify these
steps by providing packages through the homebrew_ package manager. On
Unix/Linux, many of these packages are available through your 
distribution-specific package manager or you can install from source.

Needed Python packages
**********************

* python 2.7 (installed by default on OSX 10.7+)
* numpy (>= 1.5)
* biopython (>= 1.54)
* seqtools

* dendropy (optional - used in some programs)
* ete2 (optional - used in some programs)
* scipy (optional - used in some programs)
* bx-python (optional - used in some programs)

Needed external applications
****************************

All external applications should be placed in your `$PATH`, so that they are
available to the code within the PHYLUCE_ package.

* lastz
* mafft
* velvet and velvetoptimiser (for contig assembly)
* bioperl (required by velvetoptimiser)

* splitaake (optional - used for demultiplexing Illumina reads)
* illumiprocessor (optional - used for trimming and renaming Illumina reads)
* AMOS (optional - used for coverage calculations of assembled UCE loci)
* muscle (optional - used as an alignment option for UCE contigs)
* dialign (optional - used as an alignment option for UCE contigs)
* saté (optional - used for alternate alignment of UCE contigs)
* raxml (optional - used for tree inference)
* mrbayes (optional - used for tree inference)
* cloudforest (optional - used for tree inference)

Useful tools
************

* "Kent Source" packages
* A good text editor (VIM, textmate2, etc.)


Software Installation (OSX)
===========================

