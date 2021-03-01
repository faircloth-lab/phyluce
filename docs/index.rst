.. phyluce documentation master file, created by
   sphinx-quickstart on Sat Jan 22 09:04:05 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. include:: global.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

phyluce:  software for UCE (and general) phylogenomics
======================================================

Release v\ |version|

:Author: Brant C. Faircloth
:Date: |date|
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

phyluce_ (phy-**loo**-chee) is a software package that was initially developed
for analyzing data collected from ultraconserved elements in organismal genomes
(see :ref:`References` and http://ultraconserved.org for additional
information).

The package includes a number of tools spanning:

* the assembly of raw read data to contigs
* the separation of UCE loci from assembled contigs
* parallel alignment generation, alignment trimming, and alignment data summary
  methods in preparation for analysis
* SNP calling and contig correction using raw-read data

As it stands, the phyluce_ package is useful for analyzing both data collected
from UCE loci and also data collection from other types of loci for phylogenomic
studies at the species, population, and individual levels.

Contributions
--------------

phyluce_ is open-source (see :ref:`License`) and we welcome contributions from
anyone who is interested.  Please make a pull request on github_.

Issues
------
The issue tracker for phyluce_ is `available on github <https://github.com/faircloth-lab/phyluce/issues>`_.
If you have an issue, please ensure that you are experiencing this issue on a 
supported OS (see :ref:`Installation`) using the conda_ installation of
phyluce_. When submitting issues, please include a test case demonstrating 
the issue and indicate which operating system and phyluce version
you are using.


Guide
=====

.. toctree::
   :caption: Contents:
   :maxdepth: 3

   purpose
   installation
   tutorials/index
   daily-use/index


Project info
============
.. toctree::
   :maxdepth: 1

   citing
   license
   attributions
   funding
   ack
