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
`dependency chain <http://en.wikipedia.org/wiki/Dependency_hell>`_  (the
programs that phyluce_ requires to run) is reasonably complex.

In the current versions (> 1.7.x), **and we very strongly suggest** that 
users install phyluce_ using the miniconda_ Python distribution.

.. attention:: We do not support installing phyluce through means other than the
    conda_ installer.  This means that we do not test phyluce_ against any
    binaries, other than those we build and distribute through conda_.
    Although you can configure phyluce_ to use binaries
    of different provenance, this is not officially supported.

.. note:: We build and test the binaries available through conda_ using
    64-bit operating systems that include:

    - MacOS 10.15
    - Ubuntu 20.04 LTS

phyluce_ is also available for use as a docker_ image.  Underneath the hood
the docker_ image runs Ubuntu 20.04 LTS and installs phyluce_ and related
packages using conda_.


Install Process
===============

.. attention:: We do not support phyluce_ on Windows, although you technically
    should be able to install phyluce_ on Windows using the Windows Subsystem
    for Linux (WSL) and installing Ubuntu 20.04 LTS from the Windows Store. 
    You should also be able to use the docker_ image.

Using Conda
-----------

.. note:: We build and test the binaries available through conda_ using
    64-bit operating systems that include the following.  We will officially
    support MacOS 10.16 when the github_ build system offers this platform
    for automated tests.

    - MacOS 10.15
    - Ubuntu 20.04 LTS

The installation process is a 2-step process.  You need to:

#. Install miniconda_
#. Install phyluce_

Installing phyluce_ will install all of the required binaries, libraries, and
Python_ dependencies that you need to run the program.


Install miniconda
^^^^^^^^^^^^^^^^^

First, you need to install miniconda_. Follow the instructions for your platform that
are available from `conda.io <https://conda.io/docs/user-guide/install/index.html>`_. 
After you have run the install process **be sure** that you:

#. close and re-open your terminal window
#. run ``conda list`` which should produce output

Install phyluce
^^^^^^^^^^^^^^^

Current practice with conda_ is to keep all environments separate and **not**
to use the base environment as a "default" environment.  So, we will be 
installing phyluce_ into an environment named ``phyluce-vX.X`` where the 
``X.x`` represents the version you choose.

#. Go to the `phyluce github release page <https://github.com/faircloth-lab/phyluce/releases>`_
#. Download the appropriate ``*.yml`` file for the phyluce_ version 
   you want and the **operating system you are using**
#. Install that into an environment corresponding to the phyluce 
   version, e.g. ``phyluce-1.7`` following the instructions on the phyluce_ release page

This will create an environment named ``phyluce-X.x``, then download and install
everything you need to run phyluce_ into this ``phyluce-X.x`` conda_ environment. 

To use your new phyluce_ environment, you **must** run (replace ``X.x`` with the
correct version):

.. code-block:: bash

    conda activate phyluce-X.x

To stop using this phyluce environment, you **must** run:

.. code-block:: bash

    conda deactivate

What conda installs
^^^^^^^^^^^^^^^^^^^

When you install phyluce, it specifies a number of dependencies that it needs
to run.  If you would like to know **everything** that conda_ has 
installed, you can open up the ``*.yml`` you downloaded (it is simply a text file)
and take a look at the contents.

From within the conda_ environment, you can also run

.. code-block:: bash

   conda activate phyluce-X.x
   conda list

Added benefits
^^^^^^^^^^^^^^

An added benefit of using conda_ is that you can also run all of the 3rd-party
binaries without worrying about setting the correct $PATH, etc.

For example, phyluce_ requires MUSCLE for installation, and MUSCLE was installed
by conda_ as a dependency of phyluce_. Because conda puts all of these binaries
in our ``$PATH`` when the environment is activateed, we can also just run MUSCLE
on the command-line, with, e.g.,:

.. code-block:: bash

    $ muscle -version

    MUSCLE v3.8.1551 by Robert C. Edgar

Using Docker
------------

We also provide phyluce_ as a docker_ image, which means you can run the phyluce_
installation anywhere that you can run docker_. The docker_ image is built on
Ubuntu 20.04 LTS using conda_.  To pull the docker image:

#. Go to the `phyluce github release page <https://github.com/faircloth-lab/phyluce/releases>`_
#. Find the phyluce_ release you want (usually the most recent)
#. Run the ``docker pull`` command listed

Although using docker_ is beyond the scope of this guide, you can run phyluce_
within the docker container with, e.g.:

.. code-block:: bash

    docker run phyluce <phyluce_program_name>

When you run this, all commands are run in the default directory ``/work``.  So,
you will very likely want to mount a local directory to the docker container like
so:

.. code-block:: bash

    docker run -v $(pwd):/data phyluce <phyluce_program_name>

Which will store the output of the commands you are running the direcotry 
from which you ran the command.


phyluce configuration
=====================

As of v1.5, phyluce_ uses a configuration file to keep track of paths to relevant
binaries, as well as some configuration information.  This file is located at
``$CONDA_PREFIX/phyluce/config``.  Although you can edit this file directly, you
can also create a user-specific configuration file at `~/.phyluce.conf` (**note the
preceding dot**), which will override the default values with different paths.

So, if you need to use a slightly different binary or you want to experiment
with new binaries (e.g. for assembly), then you can change the paths in this
file rather than deal with hard-coded paths.

.. attention:: This **WILL NOT** work for the docker_ image by default.  You 
    also do NOT **need** to to anything with this file - $PATHs should 
    automatically resolve.

.. warning:: Changing the `$PATHs` in the config file can break things pretty
    substantially, so please use with caution. If you are making changes, 
    edit the copy at ``~/.phyluce.conf``) rather than the default copy.

The format of the config file as of v1.7 looks like the following:

.. code-block:: bash

    [binaries]
    abyss:$CONDA/bin/ABYSS
    abyss-pe:$CONDA/bin/abyss-pe
    bcftools:$CONDA/bin/bcftools
    bedtools:$CONDA/bin/bedtools
    bwa:$CONDA/bin/bwa
    gblocks:$CONDA/bin/Gblocks
    lastz:$CONDA/bin/lastz
    mafft:$CONDA/bin/mafft
    muscle:$CONDA/bin/muscle
    pilon:$CONDA/bin/pilon
    raxml-ng:$CONDA/bin/raxml-ng
    samtools:$CONDA/bin/samtools
    seqtk:$CONDA/bin/seqtk
    spades:$CONDA/bin/spades.py
    trimal:$CONDA/bin/trimal
    velvetg:$CONDA/bin/velvetg
    velveth:$CONDA/bin/velveth
    snakemake:$CONDA/bin/Snakemake

    [workflows]
    mapping:$WORKFLOWS/mapping/Snakefile
    correction:$WORKFLOWS/contig-correction/Snakefile
    phasing:$WORKFLOWS/phasing/Snakefile

    #----------------
    #    Advanced
    #----------------

    [headers]
    trinity:comp\d+_c\d+_seq\d+|c\d+_g\d+_i\d+|TR\d+\|c\d+_g\d+_i\d+|TRINITY_DN\d+_c\d+_g\d+_i\d+
    velvet:node_\d+
    abyss:node_\d+
    idba:contig-\d+_\d+
    spades:NODE_\d+_length_\d+_cov_\d+.\d+

    [spades]
    max_memory:4
    cov_cutoff:5


Other useful tools
==================

You will need to be familiar with the command-line/terminal, and it helps to
have a decent text editor for your platform.  Here are some suggestions that
are free:

- `vscode <https://code.visualstudio.com/>`_
- `atom <https://atom.io/>`_
