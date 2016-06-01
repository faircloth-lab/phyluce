.. phyluce documentation master file, created by
   sphinx-quickstart on Sat Jan 22 09:04:05 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. include:: global.rst
.. |date| date:: %d %B %Y %H:%M %Z (%z)

phyluce:  software for UCE (and general) phylogenomics
======================================================

Release v\ |version|. (:ref:`Changelog`)

:Author: Brant C. Faircloth
:Date: |date|
:Copyright: This documentation is available under a Creative Commons (`CC-BY`_) license.

phyluce_ (phy-**loo**-chee) is a software package that was initially developed
for analyzing data collected from ultraconserved elements in organismal genomes
(see :ref:`References` and http://ultraconserved.org for additional
information).

The package now includes a number of tools spanning:

* the assembly of raw read data to contigs
* the separation of UCE loci from assembled contigs
* parallel alignment generation, alignment trimming, and alignment data summary
  methods in preparation for analysis
* alignment and SNP calling using UCE or other types of raw-read data.

As it stands, the phyluce_ package is useful for analyzing both data collected
from UCE loci and also data collection from other types of loci for phylogenomic
studies at the species, population, and individual levels.

Contributions
--------------

phyluce_ is open-source (see :ref:`License`) and we welcome contributions from
anyone who is interested.  Please make a pull request on github_.  The issue
tracker for phyluce_ is also `available on github <https://github.com/faircloth-
lab/phyluce/issues>`_.

Issues
------

If you have an issue, please ensure that you are experiencing this issue on a
supported OS (see :ref:`Installation`) using the conda_ installation of
phyluce_.  If possible, please submit a test case demonstrating the issue and
indicate which platform, git checkout, and phyluce version you are using.


Guide
=====

.. toctree::
   :maxdepth: 2

   purpose
   installation
   quality-control
   assembly
   uce-processing
   tutorial-one
   tutorial-three
   tutorial-four

Project info
============
.. toctree::
   :maxdepth: 1

   citing
   license
   changelog
   attributions
   funding
   ack

Supporting documents
====================

.. toctree::
   :maxdepth: 2

   list-of-programs
