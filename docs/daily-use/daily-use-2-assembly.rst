.. include:: ../global.rst

.. _Assembly:

*********
Assembly
*********

Setup
=====

Once your reads are clean, you're ready to assemble. At the moment, you can use
velvet_, ABySS_, and spades_ for assembly.

Most of the assembly process is automated using code within phyluce_,
specifically the following 3 scripts:

- `phyluce_assembly_assemblo_abyss`
- `phyluce_assembly_assemblo_spades`
- `phyluce_assembly_assemblo_velvet`

The code of each of the above programs **always** expects your input directories
to have the following structure (from the :ref:`Quality Control` section)::

    uce-clean/
        genus_species1/
            adapters.fasta
            raw-reads/
                genus_species1-READ1.fastq.gz (symlink)
                genus_species1-READ2.fastq.gz (symlink)
            split-adapter-quality-trimmed/
                genus_species1-READ1.fastq.gz
                genus_species1-READ2.fastq.gz
                genus_species1-READ-singleton.fastq.gz
            stats/
                genus_species1-adapter-contam.txt

And, each of these assembly helper programs take the same configuration files as
input.  You should format the configuration file for input according to the
following scheme::

    [samples]
    name_you_want_assembly_to_have:/path/to/uce-clean/genus_species1

In practice, this means you need to create a configuration file that looks
like::

    [samples]
    anas_platyrhynchos1:/path/to/uce-clean/anas_platyrhynchos1
    anas_carolinensis1:/path/to/uce-clean/anas_carolinensis1
    dendrocygna_bicolor1:/path/to/uce-clean/dendrocygna_bicolor1

The `assembly name` on the left side of the colon can be whatever you want.
The `path name` on the right hand side of the colon must be a valid path to a
directory containing read data in a format similar to that described above.

.. attention:: Assembly names **MUST** be unique.

.. admonition:: Question: How do I name my samples/assemblies?
    :class: admonition tip

    Naming samples is a contentious issue and is also a hard thing to deal with
    using computer code.  You should **never** have a problem if you name your
    samples as follows, where the genus and specific epithet are separated by
    an underscore, and multiple individuals of a given species are indicated
    using a trailing integer value::

        anas_platyrhynchos1
        anas_carolinensis1
        dendrocygna_bicolor1

    You should also not have problems if you use a naming scheme that suffixes
    the species binomial(s) with an accession number that is **simply**
    formatted (e.g. no slashes, dashes, etc.)::

        anas_platyrhynchos_KGH2267
        anas_carolinensis_KGH2269
        dendrocygna_bicolor_DWF4597
    
    The above is the recommended working format. When you search for UCE contigs
    phyluce_ should screen your taxon name to ensure they do not contain 
    restricted characters.  This includes ``.+:"'-?!*@%^&#=/\`` or names that **begin**
    with a number.  It's probably best to get that all squared away now.


Running the assembly
====================

Once your configuration file is created (best to use a decent text editor that
will not cause you grief), you are ready to start assembling your read data into
contigs that we will search for UCEs.  The code to do this for the three helper
scripts is below.

General process
---------------

The general process that the helper scripts use is:

#. Create the output directory (AKA $ASSEMBLY, below)
#. Create a ``contigs`` folder within the output directory
#. For each taxon create ``$ASSEMBLY/genus-species`` directory, based on config
   file entries
#. Find the correct fastq files for a given sample
#. Input those fastq files to whichever assembly program
#. Assemble reads
#. Strip contigs of potentially problematic bases (ABySS-only)
#. Normalize contig names
#. Link all assembly files with normalized names in $ASSEMBLY/genus-species/ into
   $ASSEMBLY/contigs/genus-species.contigs.fasta, so that all assemblies are linked
   in the same output directory.

velvet
------

.. code-block:: bash

    # make a directory for log files
    mkdir log
    # run the assembly
    phyluce_assembly_assemblo_velvet \
        --config config_file_you_created.conf \
        --output /path/where/you/want/assemblies \
        --kmer 35 \
        --subfolder split-adapter-quality-trimmed \
        --cores 12 \
        --clean \
        --log-path log

Results
^^^^^^^

The directory structure created for velvet_-based assemblies looks like::

    path-to-output-directory/
        contigs/
            genus-species1 -> ../genus-species1/out_k31/contigs.fa
        genus-species1/
            contigs.fasta -> out_k31/contigs.fa
            out_k31
            velvetg-k31.err.log
            velvetg-k31.out.log
            velveth-k31.err.log
            velveth-k31.out.log

ABySS
-----

.. code-block:: bash

    # make a directory for log files
    mkdir log
    # run the assembly
    phyluce_assembly_assemblo_abyss \
        --config config_file_you_created.conf \
        --output /path/where/you/want/assemblies \
        --kmer 35 \
        --subfolder split-adapter-quality-trimmed \
        --cores 12 \
        --clean \
        --log-path log

.. attention:: Following assembly, `phyluce_assembly_assemblo_abyss` modifies
    the assemblies by replacing degenerate base codes with standard nucleotide
    encodings.  We do this because lastz_, which we use to match contigs to
    targeted UCE loci, is not compatible with degenerate IUPAC codes.

    The phyluce_assembly_assemblo_abyss code makes these substitutions for every site having a
    degenerate code by selecting the appropriate nucleotide encoding randomly.
    The code also renames the ABySS assemblies using the velvet_ naming
    convention.  The modified contigs are them symlinked into
    ``$ASSEMBLY/contigs``.  Unmodified contigs are available in ``$ASSEMBLY/genus-
    species/out_k*-contigs.fa``


Results
^^^^^^^

The directory structure created for ABySS_-based assemblies looks like::

    path-to-output-directory/
        contigs/
            genus-species1 -> ../genus-species1/out_k31-contigs-velvet.fa
        genus-species1/
            abyss-k31.err.log
            contigs.fasta -> out_k31-contigs-velvet.fa
            out_k31-contigs.fa
            out_k31-scaffolds.fa
            out_k31-unitigs.fa
            abyss-k31.out.log
            coverage.hist
            out_k31-contigs-velvet.fa
            out_k31-stats


Spades
------

.. code-block:: bash

    # make a directory for log files
    mkdir log
    # run the assembly
    phyluce_assembly_assemblo_spades \
        --config config_file_you_created.conf \
        --output /path/where/you/want/assemblies \
        --subfolder split-adapter-quality-trimmed \
        --clean \
        --cores 12 \
        --log-path log

.. admonition:: Question: What is the –clean option?
    :class: admonition tip
    The ``--clean`` option removes all extraneous files produced by spades_.


Results
^^^^^^^

The directory structure created for spades_-based assemblies looks like::

    path-to-output-directory/
        contigs/
            genus-species1 -> ../genus-species1/scaffolds.fasta
        genus-species1/
            contigs.fasta -> Trinity.fasta
            Trinity.fasta
            trinity.log


Common questions
----------------

.. admonition:: Question: Which assembly program do I pick?
    :class: admonition tip

    Generally, I would suggest that you use spades_.  It produces
    reasonable contig assemblies that are longer than the assemblies built by
    velvet_, ABySS_, or Trinity_ (now removed from phyluce_).  It arguable
    produces assemblies that are more accurate than assemblies from
    these other programs.


.. admonition:: Question: For ABySS and velvet, what --kmer value do I use?
    :class: admonition tip

    Also a hard question.  Part of the reason that it is hard is due to the fact
    that we are trying to assemble data of heterogenous read depth (i.e., our
    reads are spread across (mostly) UCE loci, but the depth of coverage of each
    locus is varaible due to capture efficiency).  Longer kmer values can give
    you longer (but fewer) contigs, while shorter kmer values produce fewer,
    more abundant contigs.  In most cases, your assemblies will be decent with a
    kmer value around 55-65.
