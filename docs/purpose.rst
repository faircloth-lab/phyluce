.. include:: global.rst

Purpose
=======

`Phylogenomics`_ offers the possibility of helping to resolve the `Tree of
Life`_.  To do this, we often need either very cheap sources of organismal
genome data or decent methods of subsetting larger genomes (e.g., vertebrates; 1
Gbp) such that we can collect data from across the genome in an efficient and
economical fashion, find the regions of each genome that are shared among
organisms, and attempt to infer the evolutionary history of the organisms in
which we're interested using the data we collect.

Genome reduction techniques offer one way to collect these types of data from
both small- and large-genome organisms.  These "reduction" techniques include
various flavors of `amplicon sequencing`_, `RAD-seq`_ (**R**\ estriction site
**A**\ ssociated **D**\ NA markers), `RNA-seq`_ (transcriptome sequencing), and
sequence capture methods.

phyluce_ is a software package for working with data generated from sequence
capture of UCE (**u**\ ltra-\ **c**\ onserved **e**\ lement) loci, as first
published in [BCF2012]_.  Specifically, phyluce_ is a suite of programs to:

* assemble raw sequence reads from Illumina platforms into contigs
* determine which contigs represent UCE loci
* filter potentially paralagous UCE loci
* generate different sets of UCE loci across taxa of interest

Additionally, phyluce_ is capable of the following tasks, which are generally
suited to any number of phylogenomic analyses:

* produce large-scale alignments of these loci in parallel
* manipulate alignment data prior to further analysis
* convert alignment data between formats
* compute statistics on alignments and other data

phyluce_ is written to process data/individuals/samples/species in parallel,
where possible, to speed execution.  Parallelism is achieved through the use
of the Python_ `multiprocessing`_ module, and most computations are suited to
workstation-class machines or servers (i.e., rather than clusters).  Where
cluster-based analyses are needed, phyluce_ will produce the necessary outputs
for input to the cluster/program that you are running so that you can setup
the run according to your cluster design, job scheduling system, etc.  Clusters
are simply too heterogenous to do a good job at this part of the analytical
workflow.

Short-term goals (`cli branch`_, v2.0.0+)
-----------------------------------------

We are currently working on a new release (this documentation) to:

* ease the burden of installing dependencies using conda_
* simplify the `CLI`_ (command-line interface) of phyluce_
* standardize **HOW** you run various analyses
* improve test coverage of the code by `unittests`_
* improve and standardize the documentation

Longer-term goals (v2.0.x and beyond)
-------------------------------------

We are also working towards adding:

* SNP-calling pipelines (sensu [BTS2013]_)
* sequence capture bait design
* identification of UCE loci

Much of this code is already written and in use by several of the
:ref:`Contributors`. As we test and improve these functions, we will add them
to the code in the future.

Who wrote this?
---------------

This documentation was written primarily by Brant Faircloth
(http://faircloth-lab.org). Brant is also responsible for the development of
most of the phyluce_ code.  Faults within the code are his.

You can find additional authors and contributors in the :ref:`Authors` and
:ref:`Contributors` sections.
