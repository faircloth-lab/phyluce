.. include:: global.rst

.. _Installation:

************
Installation
************

phyluce_ uses a number of tools that allow it to assemble data, search for UCE
loci, align results reads, manipulate alignments, prepare alignments for
analysis, etc.  To accomplish these goals, phyluce_ uses wrappers around a
number of programs that do each of these tasks (sometimes phyluce can use
several different programs that accomplish the same task in different ways).
As a result, the
`dependency chain <http://en.wikipedia.org/wiki/Dependency_hell>`_ , or the
programs that phyluce_ requires to run, is reasonably complex.

In previous versions (< 1.4.x), we required users to install a number of
different software packages, get those into the user's `$PATH`, cross their
fingers, and hope that everything ran smoothly (it usually did not).

In the current versions (> 1.4.x), we have removed a number of dependencies,
**and we very strongly suggest** that users install phyluce_ using either the
anaconda_ or miniconda_ Python distributions, along with bioconda_.

.. attention:: We do not support installing phyluce through means other than the
    conda_ installer.  This means that we do not test phyluce_ against any
    binaries, other than those we build and distribute through conda_.
    You will eventually be able to configure phyluce_ to use binaries
    of different provenance, although this will not be officically supported,
    other than providing a mechanism to do so.

.. note:: We build and test the binaries available through conda_ using
    64-bit operating systems that include:

    - Apple OSX 10.9.x
    - CentOS 6.x, 7.x
    - Ubuntu 14.04 LTS

Why conda?
==========

It may seem odd to impose a particular disitribution on users, and we largely
agree.  However, conda_ makes it very easy for us to distribute both Python_ and
non-Python packages (e.g. velvet, ABySS, etc.), setup identical environments
across very heterogenous platforms (linux, osx), make sure all the `$PATHs` are
correct, and have things run largely as expected. Using conda_ has several other
benefits, including environment separation similar to virtualenv_. In short,
using conda_ gets us as close to a "one-click" install that we will probably
ever get.


Install Process
===============

.. attention:: We do not support phyluce_ on Windows.

.. note:: We build and test the binaries available through conda_ using
   64-bit operating systems that include:

   - Apple OSX 10.9.x
   - CentOS 6.x, 7.x
   - Ubuntu 14.04 LTS

The installation process is a 3-step process.  You need to:

#. Install conda_ (either anaconda_ or miniconda_)
#. Edit your ``~/.condarc`` to add the necessary bioconda repositories
#. Install phyluce_

Installing phyluce_ will install all of the required binaries, libraries, and
Python_ dependencies.


Install Anaconda or miniconda
-----------------------------

First, you need to install anaconda_ or miniconda_ **with Python 2.7**.  Whether
you choose miniconda_ or anaconda_ is up to you, your needs, how much disk
space you have, and if you are on a fast/slow connection.

.. attention:: You can easily install anaconda_ or miniconda_ in your ``$HOME``,
    although you should be aware that this setup can sometimes cause problems in
    cluster-computing situations.

.. tip:: Do I want anaconda_ or miniconda_?
    :class: admonition tip

    The major difference between the two python distributions is that anaconda_
    comes with many, many packages pre-installed, while miniconda_ comes with
    almost zero packages pre-installed.  As such, the beginning anaconda_
    distribution is roughly 200-500 MB in size while the beginning miniconda_
    distribution is 15-30 MB in size.

    **We suggest that you install miniconda.**

.. tip:: What version of miniconda_ or anaconda_ do I need?
    :class: admonition tip

    Right now, phyluce_ **only runs with Python 2.7**.  This means that you need
    to install a version of miniconda_ or anaconda_ that uses Python 2.7.  The
    easiest way to do this is to choose carefully when you download a
    particular distribution for your OS (be sure to choose the Python 2.7
    version).


miniconda
^^^^^^^^^

Follow the instructions here for your platform:
https://conda.io/docs/user-guide/install/index.html

.. note:: Once you have installed either Miniconda or Anaconda, we will refer
    to the install as `conda` throughout the remainder of this documentation.

anaconda
^^^^^^^^

Follow the instructions here for your platform:
http://docs.continuum.io/anaconda/install.html


.. note:: Once you have installed either Miniconda or Anaconda, we will refer
    to the install as `conda` throughout the remainder of this documentation.


Checking your `$PATH`
^^^^^^^^^^^^^^^^^^^^^

Regardless of whether you install miniconda_ or anaconda_, you need to check
that you've installed the package correctly.  To ensure that the correct
location for anaconda_ or miniconda_ are added to your $PATH (this occurs
automatically on the $BASH shell), run the following::

    $ python -V

The output should look similar to (`x` will be replaced by a version)::

    Python 2.7.x :: Anaconda x.x.x (x86_64)

Notice that the output shows we're using the `Anaconda x.x.x` version of
Python_. If you do not see the expected output (or something similar), then you
likely need to edit your $PATH variable to add anaconda_ or miniconda_.

The easiest way to edit your path, if needed is to open ``~/.bashrc`` with a
text editor (if you are using ZSH, this will be ``~/.zshrc``) and add, as the
last line::

    export PATH=$HOME/path/to/conda/bin:$PATH

where ``$HOME/path/to/conda/bin`` is the location of anaconda/miniconda on your
system (usually ``$HOME/anaconda/bin`` or ``$HOME/miniconda/bin``).

.. warning:: If you have previously set your ``$PYTHONPATH`` elsewhere in your
   configuration, it may cause problems with your anaconda_ or miniconda_
   installation of phyluce_.  The solution is to remove the offending library
   (-ies) from your ``$PYTHONPATH``.


Add the necessary bioconda repositories to conda
------------------------------------------------

You need to add the location of the bioconda_ repositories to your conda_
installation.  To do that, you can follow the instructions `at the bioconda
site <https://bioconda.github.io/#set-up-channels>`_ or you can simply edit
your ``~/.condarc`` file to look like:

.. code-block:: bash

    channels:
      - defaults
      - conda-forge
      - bioconda

Once you do this, you have access to all of the packages installed at
bioconda_ and conda-forge_.  The order of this file is important - conda_ will
first search in it's default repositories for package, then it will check
conda-forge, finally it will check bioconda.

How to install phyluce
----------------------

You now have two options for installing phyluce_.  You can install phyluce_ in
what is known as a `conda environment <https://conda.io/docs/user-guide/tasks
/manage-environments.html>`_, which lets you keep code for different
applications separated into different environments.  **We suggest this
route**.

You can also install all of the phyluce_ code and dependencies in
your default conda_ environment.

Install phyluce in it's own conda environment
---------------------------------------------

We can install everything that we need for phyluce_ in it's own environment by running:

.. code-block:: bash

    conda create --name phyluce phyluce

This will create an environment named ``phyluce``, then download and install
everything you need to run phyluce_ into this `phyluce` conda environment. To
use this phyluce environment, you **must** run:

.. code-block:: bash

    conda activate phyluce

To stop using this phyluce environment, you **must** run:

.. code-block:: bash

    conda deactivate

Install phyluce in the default conda environment
------------------------------------------------

We can simply install everything that we need in our default conda_
environment, as well.  In some ways, this is easier, but it could be viewed as
a less-ideal option in terms of repeatability and separability of functions.
To install phyluce_ in the default environment, after making sure that you
have miniconda_ or anaconda_ in your $PATH, and after adding the bioconda
repositories, run:

.. code-block:: bash

    conda install phyluce

If you need GATK
================

GATK changed its licensing policies, which means that there are some extra
steps you need to take to install GATK alongside phyluce.  First, follow `this
link to download GATK 3.5
<https://software.broadinstitute.org/gatk/download/auth?package=GATK-
archive&version=3.5-0-g36282e4>`_.  Once that is downloaded, you need to
unzip/expand the archive, and ``GenomeAnalysisTK.jar`` will be inside.  To
install this, run the follwing:

.. code-block:: bash

    # if phyluce is installed in its own environment (if not, skip this)
    source activate phyluce

    # install GATK
    gatk-register /path/to/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar

That should take care of everything you need, and you should be able to run GATK on the command-line.


What conda installs
===================

When you install phyluce, it specifies a number of dependencies that it needs
to run.  conda_ is great because it will pull specific **versions** of the
required programs from the bioconda_ repository and install those on your machine,
setup the paths, etc.

Below is a list of what phyluce_ currently (1.6.2) installs:

3rd-party dependencies and packages installed
---------------------------------------------

.. code-block:: text

    - python
    - abyss 1.5.2
    - bcftools
    - bedtools
    - biopython
    - bwa
    - bx-python
    - dendropy 3.12.3
    - gblocks
    - lastz
    - mafft
    - muscle
    - pandas
    - picard
    - pysam
    - pyvcf
    - raxml
    - samtools
    - seqtk
    - trimal
    - trinity # [not osx]
    - velvet
    - illumiprocessor
    - spades
    - itero


Added benefits
==============

An added benefit of using conda_ and installing packages in this way is that you
can also run all of the 3rd-party binaries without worrying about setting the
correct $PATH, etc.

For example, phyluce_ requires MUSCLE for installation, and MUSCLE was installed
by conda_ as a dependency of phyluce_. Because ``$HOME/anaconda/bin`` (which we
will now call `$CONDA`) is part of our path now, and because phyluce_ installed
MUSCLE, we can also just run MUSCLE on the command-line, with::

    $ muscle

    MUSCLE v3.8.31 by Robert C. Edgar

    http://www.drive5.com/muscle
    This software is donated to the public domain.
    Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.


    Basic usage

        muscle -in <inputfile> -out <outputfile>

    Common options (for a complete list please see the User Guide):

        -in <inputfile>    Input file in FASTA format (default stdin)
        -out <outputfile>  Output alignment in FASTA format (default stdout)
        -diags             Find diagonals (faster for similar sequences)
        -maxiters <n>      Maximum number of iterations (integer, default 16)
        -maxhours <h>      Maximum time to iterate in hours (default no limit)
        -html              Write output in HTML format (default FASTA)
        -msf               Write output in GCG MSF format (default FASTA)
        -clw               Write output in CLUSTALW format (default FASTA)
        -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
        -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
        -quiet             Do not write progress messages to stderr
        -version           Display version information and exit

    Without refinement (very fast, avg accuracy similar to T-Coffee): -maxiters 2
    Fastest possible (amino acids): -maxiters 1 -diags -sv -distance1 kbit20_3
    Fastest possible (nucleotides): -maxiters 1 -diags

This is true for other binaries you install from our repository (e.g. velveth,
velvetg, abyss-pe, mafft) or any other conda_ repository - those binaries are
all stored in ``$CONDA/bin``.


$PATH configuration
===================

As of v1.5, phyluce_ uses a configuration file to keep track of paths to relvant
binaries, as well as some configuration information.  This file is located at
`$CONDA/config/phyluce.conf`.  Although you can edit this file directly, you can
also create a user-specific configuration file at `~/.phyluce.conf` (**note the
preceding dot**), which will override the default values with different paths.

So, if you need to use a slightly different binary or you want to experiment
with new binaries (e.g. for assembly), then you can change the paths in this
file rather than deal with hard-coded paths.

.. attention:: You do NOT **need** to to anything with this file - $PATHs should
    automatically resolve.

.. warning:: Changing the `$PATHs` in the config file can break things pretty
    substantially, so please use with caution (and edit the copy at
    `~/.phyluce.conf`) rather than the default copy.

The format of the config file as of v1.6 looks like the following:

.. code-block:: bash

    [binaries]
    abyss:$CONDA/bin/ABYSS
    abyss-pe:$CONDA/bin/abyss-pe
    bcftools:$CONDA/bin/bcftools
    bedtools:$CONDA/bin/bedtools
    bwa:$CONDA/bin/bwa
    gatk:$CONDA/bin/gatk
    gblocks:$CONDA/bin/gblocks
    lastz:$CONDA/bin/lastz
    mafft:$CONDA/bin/mafft
    muscle:$CONDA/bin/muscle
    picard:$CONDA/bin/picard
    raxmlHPC-SSE3:$CONDA/bin/raxmlHPC-SSE3
    raxmlHPC-PTHREADS-SSE3:$CONDA/bin/raxmlHPC-PTHREADS-SSE3
    samtools:$CONDA/bin/samtools
    seqtk:$CONDA/bin/seqtk
    spades:$CONDA/bin/spades.py
    trimal:$CONDA/bin/trimal
    trinity:$CONDA/bin/Trinity
    vcfutils:$CONDA/bin/vcfutils.pl
    velvetg:$CONDA/bin/velvetg
    velveth:$CONDA/bin/velveth

    #----------------
    #    Advanced
    #----------------

    [headers]
    trinity:comp\d+_c\d+_seq\d+|c\d+_g\d+_i\d+|TR\d+\|c\d+_g\d+_i\d+
    velvet:node_\d+
    abyss:node_\d+
    idba:contig-\d+_\d+
    spades:NODE_\d+_length_\d+_cov_\d+.\d+

    [trinity]
    max_memory:8G
    kmer_coverage:2

    [spades]
    max_memory:2
    cov_cutoff:5


Other useful tools
==================

You will need to be familiar with the command-line/terminal, and it helps to
have a decent text editor for your platform.  Here are some suggestions:

- `vscode <https://code.visualstudio.com/>`_
- `Sublime Text <https://www.sublimetext.com/>`_
- `atom <https://atom.io/>`_
