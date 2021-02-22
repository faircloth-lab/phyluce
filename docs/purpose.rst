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

Longer-term goals (v2.0.0+ and beyond)
---------------------------------------

We are also working towards adding:

* simplify the `CLI`_ (command-line interface) of phyluce_
* add additioanl ``workflows`` for multi-step analyses

Who wrote this?
---------------

This documentation was written primarily by Brant Faircloth
(http://faircloth-lab.org). Brant is also responsible for the development of
most of the phyluce_ code.  Bugs within the code are usually his.

You can find additional authors and contributors in the :ref:`Attributions`
section.

How do I report bugs?
----------------------

To report a bug, please post an issue to
https://github.com/faircloth-lab/phyluce/issues.  Please also ensure that you
are using one of the "supported" operating systems and a supported installation
method.  Please see the :ref:`Installation` section for more details.
