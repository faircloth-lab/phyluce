.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/faircloth-lab/phyluce
   :target: https://gitter.im/faircloth-lab/phyluce?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge


phyluce: software for UCE (and general) phylogenomics
-----------------------------------------------------

phyluce_ (phy-**loo**-chee) is a software package that was initially developed
for analyzing data collected from ultraconserved elements in organismal genomes.

The package now includes a number of tools spanning:

- the assembly of raw read data to contigs
- the identification of UCE contigs
- parallel alignment generation, alignment trimming, and alignment data summary
  methods in preparation for analysis
- alignment and SNP calling using UCE or other types of raw-read data.

As it stands, the phyluce_ package is useful for analyzing both data collected
from UCE loci and also data collection from other types of loci for phylogenomic
studies at the species, population, and individual levels.

Please see the `Documentation <http://faircloth-lab.github.com/phyluce/>`_ for
additional information on using phyluce_

.. _Installation:

Installing
----------

We **strongly** suggest installaing phyluce_ using anaconda_/miniconda_ and the
conda_ package manager for python_.  We build and test phyluce_ dependencies and
function on the following platforms.

- Apple OSX 10.9.x
- CentOS 6.x
- Ubuntu 14.04 LTS

Quick version
^^^^^^^^^^^^^

#. Install JAVA for your platform
#. Install anaconda_ for your platform (Instructions_)
#. Check that anaconda_/miniconda_ are installed and in your ``$PATH``::

    python -V
    Python 2.7.6 :: Anaconda 1.8.0 (x86_64)

#. Add the `faircloth-lab binstar`_ repository to your ``.condarc``::

    conda config --add channels http://conda.binstar.org/faircloth-lab

#. Install phyluce_::

    conda install phyluce

For more in-depth installation instructions, see the `Installation documents`_.

License
-------

3-clause BSD. See `License.txt`_ for more information.

Contributions
--------------

phyluce_ is open-source (see License_), and we welcome contributions from anyone
who is interested in contributing.  To contribute code, please make a pull
request on github_.  The issue tracker for phyluce_ is also `available on github
<https://github.com/faircloth-lab/phyluce/issues>`_.

Issues
------

If you have an issue, please ensure that you are experiencing this issue on a
supported OS (see :ref:`Installation`) using the conda_ installation of
phyluce_.  If possible, please submit a test case demonstrating the issue and
indicate which platform, git checkout, and phyluce version you are using.

Citing
------

If you use the `phyluce`_ code in any form, please cite the following manuscript
(which is "in review"):

.. [BCF2015] Faircloth BC. 2015. PHYLUCE is a software package for the analysis
   of conserved genomic loci.  Bioinformatics.
   doi:`10.1093/bioinformatics/btv646 <10.1093/bioinformatics/btv646>`_.

If you are processing UCE data that you have collected by targeted
enrichment using our probes/protocols, please cite the following
manuscripts, which describes the first use of the general approach:

.. [BCF2012] BC Faircloth, McCormack JE, Crawford NG, Harvey MG, Brumfield RT,
   Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers
   spanning multiple evolutionary timescales. Systematic Biology 61: 717–726.
   doi:`10.1093/sysbio/SYS004 <http://doi.org/10.1093/sysbio/SYS004>`_.

.. _phyluce: https://github.com/faircloth-lab/phyluce
.. _conda: http://docs.continuum.io/conda/
.. _anaconda: http://docs.continuum.io/anaconda/install.html
.. _miniconda: http://repo.continuum.io/miniconda/
.. _License: https://github.com/faircloth-lab/phyluce/blob/master/LICENSE.txt
.. _License.txt: https://github.com/faircloth-lab/phyluce/blob/master/LICENSE.txt
.. _Instructions: http://docs.continuum.io/anaconda/install.html
.. _Installation documents: http://phyluce.readthedocs.org/en/installation.html
.. _python: http://www.python.org
.. _faircloth-lab binstar: http://binstar.org/faircloth-lab/

