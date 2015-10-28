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
anaconda_ or miniconda_ Python distributions.

.. attention:: We do not support installing phyluce through means other than the
    conda_ installer.  This means that we do not test phyluce_ against any
    binaries, other than those we build and distribute through conda_.
    You will eventually be able to configure phyluce_ to use binaries
    of different provenance, although this will not be officically supported,
    other than providing a mechanism to do so.

.. note:: We build and test the binaries available through conda_ using
    64-bit operating systems that include:

    - Apple OSX 10.9.x
    - CentOS 6.x
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
   - CentOS 6.x
   - Ubuntu 14.04 LTS

The installation process is now a 3-step process.  You need to:

#. Install JAVA
#. Install conda_ (either anaconda_ or miniconda_)
#. Install phyluce_

Installing phyluce_ will install all of the required binaries, libraries, and
Python_ dependencies.

Install JAVA
------------

Although we're using conda_, you need to install a JAVA distribution for your
platform.  You should install **JAVA7**, as that will ensure a number of future
tools in phyluce_ will work on your platform.  Installing JAVA is a little
tricky across different platforms (and somewhat beyond the scope of this
document), but we describe, below, how we usually do it.

Apple OS X
^^^^^^^^^^

To install JAVA 1.7, download and install the JAVA 1.7 package from Oracle here:
http://www.java.com/en/download/manual.jsp

CentOS 6.5.x linux
^^^^^^^^^^^^^^^^^^^

You can install the JRE with the following `yum` command::

    su -c "yum update"
    su -c "yum install java-1.7.0-openjdk"

Ubuntu 14.04 linux
^^^^^^^^^^^^^^^^^^^

You can install the JRE with the following `apt-get` command::

    sudo apt-get update
    sudo apt-get install openjdk-7-jre


Install Anaconda or miniconda
-----------------------------

After you installed `JAVA`, you need to install anaconda_ or miniconda_.  Which
one you choose is up to you, your needs, how much disk space you have, and if
you are on a fast/slow connection.

.. attention:: You can easily install anaconda_ or miniconda_ in your $HOME,
    although you should be aware that this setup can cause problems in some
    cluster-computing situations.

.. tip:: Do I want anaconda_ or miniconda_?
    :class: admonition tip

    The major difference between the two python distributions is that anaconda_
    comes with many, many packages pre-installed, while miniconda_ comes with
    almost zero packages pre-installed.  As such, the beginning anaconda_
    distribution is roughly 200-500 MB in size while the beginning miniconda_
    distribution is 15-30 MB in size.


anaconda
^^^^^^^^

Follow the instructions here for your platform:
http://docs.continuum.io/anaconda/install.html

miniconda
^^^^^^^^^

Find the correct `miniconda-x.x.x` file for your platform from
http://repo.continuum.io/miniconda/ and download that file.  Be sure you **do
not** get one of the packages that has a name starting with `miniconda3-`. When
that has completed, run one of the following::

    bash Miniconda-x.x.x-Linux-x86_64.sh  [linux]
    bash Miniconda-x.x.x-MacOSX-x86_64.sh [osx]

.. note:: Once you have installed Miniconda, we will refer to it as **anaconda**
   throughout the remainder of this documentation.


Checking your `$PATH`
^^^^^^^^^^^^^^^^^^^^^

Regardless of whether you install anaconda_ or miniconda_, you need to check
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

Add the faircloth-lab repository to conda
-----------------------------------------

You need to add the location of the packages we need to your conda
distributions.  to do taht, you have to add the
`faircloth-lab conda repository <http://conda.binstar.org/faircloth-lab>`_
to conda.  You can do that with the following command, which automatically edits
your ``~/.condarc`` file)::

    conda config --add channels https://conda.binstar.org/faircloth-lab


Install phyluce
---------------

Now, you install phyluce_ and all of it's dependencies by running::

    conda install phyluce

This step will add the phyluce_ library to your $PYTHONPATH and also add a
number of scripts from the phyluce_ code to your `$HOME/anaconda/bin` directory.


What conda installs
===================

When you install phyluce, it specifies a number of dependencies that it needs
to run.  conda_ is great because it will pull specific **versions** of the
required programs from the
`faircloth-lab conda repository <http://conda.binstar.org/faircloth-lab>`_ and
install those on your machine, setup the paths, etc.

Below is a list of what phyluce_ currently (v1.4.x) requires for installation.

3rd-party dependencies
----------------------

- abyss 1.3.7 (max kmer=96; built with boost and google-sparsehash)
- bowtie 1.1.1
- bedtools 2.18.1
- bwa 0.7.7
- bx-python 0.7.1
- dendropy 3.12.0
- gatk-lite 2.3.0
- gblocks 0.91b
- illumiprocessor 2.0.7
- lastz 1.02.00
- mafft 7.130
- muscle 3.8.31
- picard 1.106
- pysam 0.7.7
- pyvcf 0.6.4
- raxml 8.0.19
- samtools 0.1.19
- trinity 2.0.6
- velvet 1.2.10 (max kmer=96)

Python packages
---------------

- Python 2.7 (sets conda default Python to 2.7)
- numpy 1.7
- BioPython 1.63
- dendropy 3.12.0
- illumiprocessor 2.0.7
- phyluce 1.4.x


Added benefits
==============

An added benefit of using conda_ and installing packages in this way is that you
can also run all of the 3rd-party binaries without worrying about setting the
correct $PATH, etc.

For example, phyluce_ required MUSCLE for installation, and MUSCLE was installed
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

We have setup conda to install other files in a standard location as well.  So
JAR files are stored in ``$CONDA/jar``; libraries that you install from our repo
are stored in ``$CONDA/lib``, etc. The locations and versions are standardized
within our conda_ distribution so that we always know where things are
installed, hopefully avoiding lots of the problems with `dependency hell
<http://en.wikipedia.org/wiki/Dependency_hell>`_ and making our lives easier.

The structure of the conda repository that we use looks like the following::

    $CONDA/ (<= can be main CONDA or CONDA env)
        bin/
            Binaries we build
        config/
            Holds a default configuration file for phyluce
        docs/
            Pertinent documentation, if any.
        include/
            Header files (Qt, Boost)
        jar/
            Java Archive (JAR) Files
        lib/
            Libraries
        libexec/
            Helpers called by other programs (e.g. mafft)
        share/
            Additional support files
        test/
            Test data required by some programs.


$PATH configuration
===================

As of v1.5, phyluce_ uses a configuration file to keep track of paths to relvant
binaries, as well as some configuration information.  This file is located at
`$CONDA/config/phyluce.conf`.  Although you can edit this file directly, you can
also create a user-specific configuration file at `~/.phyluce.conf` (**note the
preceding dot**), which will override the default values for different paths.
So, if you need to use a slightly different binary or you want to experiment
with new binaries (e.g. for assembly), then you can change the paths in this
file rather than deal with hard-coded paths.

.. attention:: You do NOT **need** to to anything with this file - $PATHs should
    automatically resolve.

.. warning:: Changing the `$PATHs` in the config file can break things pretty
    substantially, so please use with caution (and edit the copy at
    `~/.phyluce.conf`) rather than the default copy.

The format of the config file as of v1.5 looks like the following:

.. code-block:: bash

    [abyss]
    abyss:$CONDA/bin/ABYSS
    abyss-pe:$CONDA/bin/abyss-pe

    [bowtie]
    bowtie:$CONDA/bin/bowtie

    [bwa]
    bwa:$CONDA/bin/bwa

    [gblocks]
    gblocks:$CONDA/bin/gblocks

    [java]
    executable:java
    mem:-Xmx8g
    jar:$CONDA/jar
    gatk:GenomeAnalysisTKLite.jar

    [lastz]
    lastz:$CONDA/bin/lastz

    [mafft]
    mafft:$CONDA/bin/mafft

    [muscle]
    muscle:$CONDA/bin/muscle

    [raxml]
    raxmlHPC-SSE3:$CONDA/bin/raxmlHPC-SSE3
    raxmlHPC-PTHREADS-SSE3:$CONDA/bin/raxmlHPC-PTHREADS-SSE3

    [samtools]
    samtools:$CONDA/bin/samtools

    [trinity]
    trinity:$CONDA/bin/Trinity
    max_memory:8G
    kmer_coverage:2

    [velvet]
    velvetg:$CONDA/bin/velvetg
    velveth:$CONDA/bin/velveth

    #----------------
    #    Advanced
    #----------------

    [headers]
    trinity:comp\d+_c\d+_seq\d+|c\d+_g\d+_i\d+|TR\d+\|c\d+_g\d+_i\d+
    velvet:node_\d+


Other useful tools
==================

You will need to be familiar with the command-line/terminal, and it helps to
have a decent text editor for your platform:

- gedit [linux]
- Sublime Text [linux, osx]
- TextMate [osx]
