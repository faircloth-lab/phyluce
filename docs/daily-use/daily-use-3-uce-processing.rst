.. include:: ../global.rst

.. _UCE Processing:

********************************
UCE Processing for Phylogenomics
********************************

The process described below is meant for users who are analyzing UCE data in
phylogenetic contexts - meaning that you are interested in addressing questions
at or deeper than the species-level.


Identifying UCE loci
====================

Once we have assembled our fastq data (see :ref:`Assembly`), we need to process
those contigs to (a) determine which represent enrichend UCE loci and (b)
remove any potential paralogs from the data set.  Before we can do that, we
need to to a little preparatory work by downloading a FASTA file representing
the bait/probe set that we used.

Get the probe set FASTA
-----------------------

To identify which of the contigs we've assembled are UCE loci (and which UCE
loci they might be), we are going to match our assembled contigs to the probes
we used to enrich UCE loci.  Before we do that, however, we need to download
a copy of probe set we used for matching purposes.

.. attention:: We archive official probe sets at
    https://github.com/faircloth-lab/uce-probe-sets/, but you need to be
    careful about which one you grab - probe sets can be of different sizes
    (e.g. 2,500 or 5,500 loci) and for different groups of taxa (e.g., amniotes,
    fish)

Download the probe set
^^^^^^^^^^^^^^^^^^^^^^^

To download a given probe set for phyluce_, you need to figure out which probe
set you need.  Then, you can use a command like ``wget`` on the command-line (or
navigate with your browser to the URL and save the file):

.. code-block:: bash

    # to get the 2.5k, amniote probe set
    wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-2.5k-probe-set/uce-2.5k-probes.fasta

    # to get the 5k, amniote probe set
    wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta


.. _contigs-matching:


Match contigs to probes
-----------------------

Once we've downloaded the probe set we used to enrich UCE loci, we need to find
which of our assembled contigs are the UCE loci that we enriched.  During this
process, the code will also remove any contigs that appear to be duplicates as a
result of assembly/other problems **or** a biological event(s).

The way that this process works is that phyluce_ aligns (using lastz_) the
contigs you assembled to the probes you input on a taxon-by-taxon (or otu-by-otu)
 basis.  Then, the code parses the alignment file to determine which contigs
matched which probes, whether any probes from a single locus matched multiple
contigs or whether a single contig matched probes designed from muliple UCE
loci.  Either of these latter two events suggests that the locus in question is
problematic.

.. hint:: **ADVANCED**: The default regular expression assumes probes in your
    file are named according to ``uce-NNN_pN``, where ``uce-`` is just a text
    string, ``NNN`` is an integer value denoting each unique locus, ``_p`` is a
    text string denoting a "probe" targeting locus ``NNN``, and the trailing
    ``N`` is an integer value denoting each unique probe targeting the same
    locus.

    If you are using a custom probe file, then you will either need to ensure
    that your naming scheme conforms to this approach **OR** you will need to
    input a different regular expression to convert the probe names to locus
    names using the ``--regex`` flag.  **It is up to you to determine what 
    is the appropriate regular expression**.

To identify which of your assembled contigs are UCE contigs, run:

.. code-block:: bash

    # make a directory for log files
    mkdir log
    # match contigs to probes
    phyluce_assembly_match_contigs_to_probes \
        --contigs /path/to/assembly/contigs/ \
        --probes uce-5k-probes.fasta \
        --output /path/to/uce/output \
        --log-path log

When you run this code, you should see output similar to::

    2014-04-24 14:38:15,979 - match_contigs_to_probes - INFO - ================ Starting match_contigs_to_probes ===============
    2014-04-24 14:38:15,979 - match_contigs_to_probes - INFO - Version: git 7aec8f1
    2014-04-24 14:38:15,979 - match_contigs_to_probes - INFO - Argument --contigs: /path/to/assembly/contigs/
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --keep_duplicates: None
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --log_path: None
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --min_coverage: 80
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --min_identity: 80
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --output: /path/to/uce/output
    2014-04-24 14:38:15,980 - match_contigs_to_probes - INFO - Argument --probes: uce-5k-probes.fasta
    2014-04-24 14:38:15,981 - match_contigs_to_probes - INFO - Argument --regex: ^(uce-\d+)(?:_p\d+.*)
    2014-04-24 14:38:15,981 - match_contigs_to_probes - INFO - Argument --verbosity: INFO
    2014-04-24 14:38:16,138 - match_contigs_to_probes - INFO - Checking probe/bait sequences for duplicates
    2014-04-24 14:38:19,022 - match_contigs_to_probes - INFO - Creating the UCE-match database
    2014-04-24 14:38:19,134 - match_contigs_to_probes - INFO - Processing contig data
    2014-04-24 14:38:19,134 - match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2014-04-24 14:38:25,713 - match_contigs_to_probes - INFO - genus_species1: 1031 (70.14%) uniques of 1470 contigs, 0 dupe probe matches, 48 UCE probes matching multiple contigs, 117 contigs matching multiple UCE probes
    2014-04-24 14:38:32,846 - match_contigs_to_probes - INFO - genus_species2: 420 (68.52%) uniques of 613 contigs, 0 dupe probe matches, 30 UCE probes matching multiple contigs, 19 contigs matching multiple UCE probes
    2014-04-24 14:38:39,184 - match_contigs_to_probes - INFO - genus_species3: 1071 (63.15%) uniques of 1696 contigs, 0 dupe probe matches, 69 UCE probes matching multiple contigs, 101 contigs matching multiple UCE probes
    2014-04-24 14:49:59,654 - match_contigs_to_probes - INFO - -----------------------------------------------------------------
    2014-04-24 14:49:59,654 - match_contigs_to_probes - INFO - The LASTZ alignments are in /path/to/uce/output/
    2014-04-24 14:49:59,654 - match_contigs_to_probes - INFO - The UCE match database is in /path/to/uce/output/probe.matches.sqlite
    2014-04-24 14:49:59,655 - match_contigs_to_probes - INFO - =============== Completed match_contigs_to_probes ===============

.. note:: The ``*.log`` files for each operation are always printed to the
    screen AND also written out to the ``$CWD`` (current working directory).
    You can keep these files more orderly by specifying a ``$LOG`` on the
    command line using the ``--log-path`` option.

Results
^^^^^^^

The resulting files will be in the::

    /path/to/output

directory. If you look in this directory, you'll see that it contains species-
specific `lastz_` files as well as an sqlite_ database::

    $ ls /path/to/output

    genus_species1.contigs.lastz
    genus_species2.contigs.lastz
    genus_species3.contigs.lastz
    probe.matches.sqlite

The ``*.lastz`` files within the ``/path/to/output`` directory are basically for
reference and individual review (they are text files that you can open using a
text editor to view).  The really important data from the lastz_ files are
summarized in the::

    probe.matches.sqlite

database file.  It's probably a good idea to have some knowledge of how this
database is structured, since it's basically what makes the next few steps work.
So, let's go over the structure and contents of this database.

.. _Database:

The probe.matches.sqlite database
.................................

``probe.matches.sqlite`` is a `relational database`_ that summarizes all
**valid** matches of contigs to UCE loci across the set of taxa that you fed it.
The database is created by and for a program named sqlite_, which is a very
handy, portable SQL database. For more info on SQL and SQLITE, see this
`sqlite-tutorial`_. I'll briefly cover the database contents and use below.

First, take a look at the contents of the database by running:

.. code-block:: bash

    sqlite3 probe.matches.sqlite

You'll now see something like::

    SQLite version 3.7.3
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite>

It's often easier to change some defaults for better viewing, so at the prompt,
paste in the following::

    sqlite> .mode columns
    sqlite> .headers on
    sqlite> .nullvalue .

.. tip:: For more info on sqlite_ "dot" commands, you can type
    ``.help``.

Now that that's done, let's see which tables the database contains by running
the ``.tables`` command::

    sqlite> .tables
    match_map  matches

This tells us there's two tables in the database, named ``match_map`` and
``matches``.


The ``matches`` table
.....................

Let's take a look at the contents of the ``matches`` table.  Once you've started
the sqlite interface, run:

.. code-block:: sql

    sqlite> SELECT * FROM matches LIMIT 10;

This query select all rows (``SELECT *``) from the ``matches`` table (``FROM
matches``) and limits the number of returned rows to 10 (``LIMIT 10``). This
will output data that look something like::

    uce         genus_species1  genus_species2  genus_species3
    ----------  --------------  --------------  --------------
    uce-500     1               .               .
    uce-501     1               .               .
    uce-502     1               .               .
    uce-503     1               1               1
    uce-504     1               .               .
    uce-505     1               .               .
    uce-506     .               .               .
    uce-507     1               .               .
    uce-508     1               1               .
    uce-509     1               1               1

Basically, what this indicates is that you enriched 9 of 10 targeted UCE loci
from ``genus_species1``, 3 of 10 UCE loci in the list from ``genus_species2``,
and 2 of 10 UCE loci from ``genus_species3``. The locus name is given in the
``uce column``.  Remember that we've limited the results to 10 rows for the sake
of making the results easy to view.

If we wanted to see only those loci that enriched in all species, we could run:

.. code-block:: sql

    sqlite> SELECT * FROM matches WHERE genus_species1 = 1
       ...> AND genus_species2 = 1 AND genus_species3 = 1;

Assuming we only had those 10 UCE loci listed above in the database, if we ran
this query, we would see something like::

    uce         genus_species1  genus_species2  genus_species3
    ----------  --------------  --------------  --------------
    uce-503     1               1               1
    uce-509     1               1               1

Basically, the ``matches`` table and this query are what we run to generate
**complete** (only loci enriched in all taxa) and **incomplete** (all loci
enriched from all taxa) matrices very easily and quickly (see
:ref:`locus-counts`).

The ``match_map`` table
.......................

The ``match_map`` table shows us which species-specific, contigs match which UCE
loci. Because each assembly program assigns an arbitrary designator to each
assembled contig, we need to map these arbitrary designators (which also differ
for each taxon/OTU) to the UCE locus to which it corresponds. Because assembled
contigs are also not in any particular orientation relative to each other across
taxa/OTUs (i.e., they may be 5' - 3' or 3' - 5'), the database also records the
orientation of all contigs relative to orientation of each probe in the probes
file.

Let's take a quick look at the ``match_map`` table:

.. code-block:: sql

    SELECT * FROM match_map LIMIT 10;

This query is similar to the one that we ran against ``matches`` and returns the
first 10 rows of the ``match_map`` table::

    uce         genus_species1  genus_species2  genus_species3
    ----------  --------------  --------------  --------------
    uce-500     node_233(+)     .               .
    uce-501     node_830(+)     .               .
    uce-502     node_144(-)     .               .
    uce-503     node_1676(+)    node_243(+)     node_322(+)
    uce-504     node_83(+)      .               .
    uce-505     node_1165(-)    .               .
    uce-506     .               .               .
    uce-507     node_967(+)     .               .
    uce-508     node_671(+)     node_211(-)     .
    uce-509     node_544(-)     node_297(+)     node_37(+)

As stated above, these results show which assembled contigs "hit" particular UCE
loci. So, if we were to open the
``$ASSEMBLY/contigs/genus_species1.contigs.fasta`` symlink the contig named
``node_1676`` corresponds to UCE locus ``uce-503``.  Because contigs are named
arbitrarily during assembly, this same UCE locus is also found in
genus_species2, but it is named ``node-243``.

Each entry in the rows also provides the orientation for particular contigs
``(-)`` or ``(+)``. This orientation is relative to the orientation of the UCE
probes/locus in the source genome (e.g., chicken for tetrapod probes).

We use this table to generate a FASTA file of UCE loci for alignment (see :ref
:`fasta-file`), after we've identified the loci we want in a particular data set
(see :ref:`locus-counts`). The code for this step also uses the associated
orientation data to ensure that all the sequence data have the same orientation
prior to alignment (some aligners will force alignment of all reads using the
given orientation rather than also trying the reverse complement and picking the
better alignment of the two).

Now that we know the taxa for which we've enriched UCE loci and which contigs
we've assembled match which UCE loci, we're ready to generate some data
matrices.

The data matrix generation process consists of two distinct parts:

#. Getting locus counts and generating a taxon set
#. Extracting FASTA data from our ``$ASSEMBLY/contigs`` based on the taxon set

.. _locus-counts:

Creating a data matrix configuration file
==========================================

After we identify the UCE loci we enriched, but before we extract fasta data
from our ``$ASSEMBLY/contigs`` corresponding to those loci, we need to create a
data matrix configuration file that denotes (1) which taxa we want to include in
a given analysis and (2) which loci will be included with this taxon set.

The taxa included in the data matrix configuration file are determined by the
user - you input a list of taxa you want to the analysis.  The UCE loci included
in the data matrix configuration are then determined by the software which
compares the requested taxa to UCE match results in ``probe.matches.sqlite`` and
two flags that you pass either one requesting **complete data matrix** or one
requesting an **incomplete data matrix**.

    complete matrix
        A phylogenetic matrix (typically sequence data) in which there are no
        missing data at any locus for any taxon/OTU.

    incomplete matrix
        A phylogenetic matrix (typocally sequence data) in which data may be
        missing from a given taxon or a given loci (or both).

During the creation of the data matrix configuration file you can also include
additional data from pre-existing UCE match databases and contigs (see :ref
:`outgroup-data`).

We'll start very simply.

Complete taxon set
------------------

First, let's generate a data matrix configuration file from only the current UCE
enrichments that will be **complete** - meaning that we will not include loci where
certain taxa have no data (either the locus was not enriched for that taxon or
removed during the filtering process for duplicate loci).

To do this, you need to create a starting taxon-configuration file (a text-based
file) denoting the taxa we want in the data set.  The taxon-configuration file
should look exactly like this (substitute in your taxon names)::

    [dataset1]
    genus_species1
    genus_species2
    genus_species3

Let's assume you save this file as ``datasets.conf``.  Now, to create the data
matrix configuration file from this taxon-configuration file, run:

.. code-block:: bash

    # create the output directory for this taxon set
    mkdir /path/to/uce/taxon-set1/

    # create the data matrix configuration file
    phyluce_assembly_get_match_counts \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --taxon-list-config datasets.conf \
        --taxon-group 'dataset1' \
        --output /path/to/uce/taxon-set1/dataset1.conf

This will basically run a query against the database, and pull out those loci
for those taxa in the `datasets.conf` file having UCE contigs.


Results
^^^^^^^

The output printed to the screen and ``$LOG`` file should look something like::

    2014-04-24 17:25:08,145 - get_match_counts - INFO - =================== Starting get_match_counts ===================
    2014-04-24 17:25:08,145 - get_match_counts - INFO - Version: git 7aec8f1
    2014-04-24 17:25:08,145 - get_match_counts - INFO - Argument --extend_locus_db: None
    2014-04-24 17:25:08,145 - get_match_counts - INFO - Argument --incomplete_matrix: False
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --keep_counts: False
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --locus_db: /path/to/uce/output/probes.matches.sqlite
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --log_path: /path/to/uce
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --optimize: False
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --output: /path/to/uce/taxon-set1/dataset1.conf
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --random: False
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --sample_size: 10
    2014-04-24 17:25:08,146 - get_match_counts - INFO - Argument --samples: 10
    2014-04-24 17:25:08,147 - get_match_counts - INFO - Argument --silent: False
    2014-04-24 17:25:08,147 - get_match_counts - INFO - Argument --taxon_group: dataset1
    2014-04-24 17:25:08,147 - get_match_counts - INFO - Argument --taxon_list_config: datasets.conf
    2014-04-24 17:25:08,147 - get_match_counts - INFO - Argument --verbosity: INFO
    2014-04-24 17:25:08,150 - get_match_counts - INFO - There are 3 taxa in the taxon-group '[dataset1]' in the config file dataset1.conf
    2014-04-24 17:25:08,151 - get_match_counts - INFO - Getting UCE names from database
    2014-04-24 17:25:08,407 - get_match_counts - INFO - There are 1314 total UCE loci in the database
    2014-04-24 17:25:11,046 - get_match_counts - INFO - Getting UCE matches by organism to generate a COMPLETE matrix
    2014-04-24 17:25:11,051 - get_match_counts - INFO - There are 306 shared UCE loci in a COMPLETE matrix
    2014-04-24 17:25:11,051 - get_match_counts - INFO -     Failed to detect 428 UCE loci in genus_species1
    2014-04-24 17:25:11,051 - get_match_counts - INFO -     Failed to detect 380 UCE loci in genus_species2
    2014-04-24 17:52:54,850 - get_match_counts - INFO - Writing the taxa and loci in the data matrix to /path/to/uce/taxon-set1/dataset1.conf
    2014-04-24 17:52:54,862 - get_match_counts - INFO - =================== Completed get_match_counts ==================

This basically says that although we've detected a total of 1,314 UCE loci in
the 3 taxa in which we are interested, when we boil those down to a complete
matrix, the complete matrix is only going to contain 306 UCE loci (of the
1,314). We had to drop 428 loci because we did not detect them in genus_species1
and we had to drop another 380 loci because we did not detect them in
genus_species2.

The output written to the ``/path/to/uce/taxon-set1/dataset1.conf`` will look
something like::

    [Organisms]
    genus_species1
    genus_species2
    genus_species3
    [Loci]
    uce-1005
    uce-1018
    uce-1025
    uce-1028
    uce-1042
    uce-1055
    uce-1060
    uce-107
    uce-1073
    uce-1074
    uce-1076
    uce-108
    ...

Taxon set membership and locus number
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now, you might think that increasing the locus count is simply a matter of
removing `genus_species1` from the list of taxa.  This is not strictly true,
however, given the vagaries of hits and misses among taxa.
`phyluce_assembly_get_match_counts`
has several other options to help you determine which taxa may be causing
problems, but picking the best combination of taxa to give you the highest
number of loci is a reasonably hard optimization problem.

.. _incomplete-matrix:

Incomplete data matrix
----------------------

You may not always want a complete data matrix. Or generating a complete matrix
drops too many loci for your tastes (it often does).  In that case, you can
easily generate an incomplete dataset using the following:

.. code-block:: bash

    # create the data matrix configuration file
    phyluce_assembly_get_match_counts \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --taxon-list-config datasets.conf \
        --taxon-group 'dataset1' \
        --output /path/to/uce/taxon-set1/dataset1.conf \
        --incomplete-matrix

.. attention:: Note the addition of the ``--incomplete-matrix`` flag.

This will generate a dataset that includes any loci enriched across the taxa
in the `datasets.conf` file.

.. note:: You do not determine the "completeness" of the final data matrix
    that you want to create during this stage - that happens later, after
    alignment (see :ref:`finalize-matrix`).  As a result, we are alinging data
    from any and all UCE loci having >= 3 taxa, which allows us to flexibly
    select the level of incompleteness later, without having to re-run our
    alignments.


Creating additional data matrix configuration files for other analyses
----------------------------------------------------------------------

If you want to generate/evaluate many data matrix configuration files containing
different taxa, you can simply create new lists within the `datasets.conf` file
like so::

    [dataset1]
    genus_species1
    genus_species2
    genus_species3

    [dataset2]
    genus_species2
    genus_species3
    genus_species4
    genus_species5
    genus_species6

And then you can run ``phyluce_assembly_get_match_counts`` against this new
section to output the data matrix configuration files:

.. code-block:: bash

    # create the data matrix configuration file
    phyluce_assembly_get_match_counts \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --taxon-list-config datasets.conf \
        --taxon-group 'dataset2' \
        --output /path/to/uce/taxon-set2/dataset2.conf

In this way, you can get some idea of how different taxon-set memberships
affect the resulting data matrix configuration files *prior to* extracting the
relevant FASTA data from ``$ASSEMBLY/contigs`` - which is a reasonably slow
process.

Incorporating outgroup/other data
---------------------------------

You may want to include outgroup data from another source into your datasets.
This can be from the pre-processed outgroup data files, but it doesn't need to
be these outgroup data. These additional data can also be contigs previously
assembled from a different set of taxa.

.. hint:: **ADVANCED**: If you want to include outgroup data from
    genome-enabled taxa, we have already created several repositories
    containing these data.  We maintaing these data under version control at:
    https://github.com/faircloth-lab/uce-probe-sets.  To download these data
    and use them in your analyses, you can clone the data using git::

        git clone https://github.com/faircloth-lab/uce-probe-sets

    Then update your ``--taxon-list-config`` file and provide the proper paths
    to the cloned data, as detailed below.

The first step of this process is to setup your ``--taxon-list-config`` slightly
differently - by indicating taxa from external data sources using asterisks::

    [dataset3]
    genus_species1
    genus_species2
    genus_species3
    genus_species4*
    genus_species5*

Here, ``genus_species4*`` and ``genus_species5*`` come from an external data
source.

Then, you need to pass ``phyluce_assembly_get_match_counts`` the location of
both **your** ``--locus-db`` and the ``--extend-locus-db``.  For example:

.. code-block:: bash

    # create the data matrix configuration file
    phyluce_assembly_get_match_counts \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --taxon-list-config datasets.conf \
        --taxon-group 'dataset3' \
        --extend-locus-db /path/to/some/other/probe.matches.sqlite \
        --output /path/to/uce/taxon-set3/dataset3.conf

To keep all this extension from getting too crazy, I've limited the ability to
include external data to a single set.  If you have lots of data from many
different enrichments, you'll need to generate a `contigs` folder containing all
these various assemblies (or symlinks to them), then align the probes to these
data (see :ref:`contigs-matching`).  Once you do that, you can extend your
current data set with all of these other data.

.. _extracting-fasta:

Extracting FASTA data using the data matrix configuration file
==============================================================

Once we have created the data matrix configuration file containing data for our
taxa of interest and those loci of interest, we need to extract the appropriate
FASTA sequences from each assembly representing the taxon/OTU of interest (e.g.
in ``$ASSEMBLY/contigs``).  This is a reasonably straightforward process that
differs only slightly based on whether you are extracting a complete matrix of
data, an incomplete matrix of data, and/or whether you are incorporating any
external data sources.


Complete data matrix
--------------------

To generate FASTA file containing the sequence data from a complete data matrix
configuration file, run:

.. code-block:: bash

    phyluce_assembly_get_fastas_from_match_counts \
        --contigs /path/to/assembly/contigs/ \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --match-count-output /path/to/uce/taxon-set1/dataset1.conf \
        --output /path/to/uce/taxon-set1/dataset1.fasta


Incomplete data matrix
----------------------

Similarly, to generate a FASTA file containing the sequence data from a complete
data matrix configuration file, run:

.. code-block:: bash

    phyluce_assembly_get_fastas_from_match_counts \
        --contigs /path/to/assembly/contigs/ \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --match-count-output /path/to/uce/taxon-set/dataset1.conf \
        --incomplete-matrix /path/to/uce/taxon-set1/dataset1.incomplete \
        --output /path/to/uce/taxon-set1/dataset1.fasta

.. attention:: Note the addition of the ``--incomplete-matrix`` option.  This
    creates an output file that contains the names of the **missing** loci by
    taxon/OTU. You can name this file anything you like.  I tend to use
    ``.incomplete`` as the extension so that it is clear what this file
    contains.

Incorporating outgroup/other data
---------------------------------

When we're incorporating external data, we need to pass the name of the external
database as well as the name of the external ``contigs``.  To generate a FASTA
file containing the sequence data from a complete data matrix configuration
that includes exeternal data sources, run:

.. code-block:: bash

    phyluce_assembly_get_fastas_from_match_counts \
        --contigs /path/to/assembly/contigs/ \
        --locus-db /path/to/uce/output/probe.matches.sqlite \
        --match-count-output /path/to/uce/taxon-set1/dataset1.conf \
        --incomplete-matrix /path/to/uce/taxon-set1/dataset1.incomplete \
        --extend-locus-db /path/to/some/other/probe.matches.sqlite \
        --extend-locus-contigs /path/to/some/other/contigs \
        --output /path/to/uce/taxon-set3/dataset3.fasta


Aligning and trimming FASTA data
================================

With all of that out of the way, things get much easier to deal with.  Now, we
need to align our data across loci, and once we're done with that, the remaining
operations we can run on the data are format-conversions, QC steps, matrix
trimming for completeness, and any number of other fun things.

Aligning the amount of data generated by enrichment approaches is reasonably
computationally intensive - so the alignment step goes fastest if you have a
multicore machine.  You also have several alignment options available, although
I would suggest sticking with MAFFT.

.. attention:: The alignment process, as implemented by phyluce_, includes
    trimming steps that trim ragged edges and remove alignments that become to
    short following trimming.

    To turn trimming off and trim alignments using another approach, pass the
    ``--no-trim`` option.  There are also several more options related to
    trimming that you can tweak.  To view these, run
    ``phyluce_align_seqcap_align --help``.

    Also see the :ref:`TutorialOne` for more info on trimming options. 

Complete data matrix
--------------------

Alignment
^^^^^^^^^^

To align the loci, by taxon, in the FASTA file you just created, run:

.. code-block:: bash

    phyluce_align_seqcap_align \
        --fasta /path/to/uce/taxon-set1/dataset1.fasta \
        --output /path/to/uce/taxon-set1/mafft-nexus/ \
        --taxa 3 \
        --aligner mafft \
        --cores 8

.. attention:: If you pass more ``--cores`` than your machine has, you will
    receive an error.

.. note:: Here, we are accepting the default, output alignment format ("nexus").
    To change that format to something else, pass the ``--output-format`` option
    with a choice of {fasta,nexus,phylip,clustal,emboss,stockholm}.

Alignment stats
^^^^^^^^^^^^^^^

Once you have alignments, it's nice to get a general sense of their length and
composition.  You can quickly (with a multicore machine) summarize thousands of
alignments by running:

.. code-block:: bash

    phyluce_align_get_align_summary_data \
        --alignments /path/to/uce/taxon-set1/mafft-nexus/ \
        --cores 12

This will produce output that looks similar to::

    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - ================ Starting get_align_summary_data ================
    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - Version: git 7aec8f1
    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - Argument --alignments: /path/to/uce/taxon-set1/mafft-nexus/
    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - Argument --cores: 12
    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - Argument --input_format: nexus
    2014-04-24 17:31:15,724 - get_align_summary_data - INFO - Argument --log_path: /path/to/uce/taxon-set1/log
    2014-04-24 17:31:15,725 - get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2014-04-24 17:31:15,725 - get_align_summary_data - INFO - Argument --verbosity: INFO
    2014-04-24 17:31:15,725 - get_align_summary_data - INFO - Getting alignment files
    2014-04-24 17:31:15,729 - get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2014-04-24 17:31:16,653 - get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] loci:    306
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] length:  223,929
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] mean:    731.79
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] 95% CI:  17.01
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] min:     275
    2014-04-24 17:31:16,654 - get_align_summary_data - INFO - [Alignments] max:     1,109
    2014-04-24 17:31:16,655 - get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2014-04-24 17:31:16,655 - get_align_summary_data - INFO - [Taxa] mean:          27.00
    2014-04-24 17:31:16,655 - get_align_summary_data - INFO - [Taxa] 95% CI:        0.00
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - [Taxa] min:           27
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - [Taxa] max:           27
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - [Missing] mean:       7.61
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - [Missing] 95% CI:     0.24
    2014-04-24 17:31:16,656 - get_align_summary_data - INFO - [Missing] min:        1.13
    2014-04-24 17:31:16,657 - get_align_summary_data - INFO - [Missing] max:        15.79
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - [All characters]      6,046,083
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - [Nucleotides]         4,924,129
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - [Matrix 50%]          306 alignments
    2014-04-24 17:31:16,661 - get_align_summary_data - INFO - [Matrix 55%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 60%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 65%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 70%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 75%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 80%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 85%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 90%]          306 alignments
    2014-04-24 17:31:16,662 - get_align_summary_data - INFO - [Matrix 95%]          306 alignments
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] '-' is present 651,009 times
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] '?' is present 470,945 times
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] 'A' is present 1,386,821 times
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] 'C' is present 1,089,729 times
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] 'G' is present 1,094,159 times
    2014-04-24 17:31:16,663 - get_align_summary_data - INFO - [Characters] 'T' is present 1,353,420 times
    2014-04-24 17:31:16,664 - get_align_summary_data - INFO - ================ Completed get_align_summary_data ===============

Locus name removal
^^^^^^^^^^^^^^^^^^

For historical reasons, and also for users to ensure that the sequence data
aligned together are from the same loci, each sequence line in the alignment
file output by ``seqcap_align`` contains the ``genus_species1`` designator,
but the ``genus_species1`` designator is also prepended with the locus name
(e.g. ``uce-1005_genus_species1``).  We need to remove these if we plan to
concatenate the loci (:ref:`raxml-concat`).  More generally, it is a good idea
to remove locus names from sequence lines before running any analyses. To do
this, run:

 .. code-block:: bash

    phyluce_align_remove_locus_name_from_files \
        --alignments /path/to/uce/taxon-set1/mafft-nexus/ \
        --output /path/to/uce/taxon-set1/mafft-nexus-clean/ \
        --taxa 3


Incomplete data matrix
----------------------

Alignment
^^^^^^^^^

The only difference for an alignment of incomplete data is that we also pass the
``--incomplete-matrix`` flag, which tells the code to expect that some loci will
not contain data across all taxa:

.. code-block:: bash

    phyluce_align_seqcap_align \
        --fasta /path/to/uce/taxon-set2/dataset2.fasta \
        --output /path/to/uce/taxon-set2/mafft-nexus/ \
        --taxa 34 \
        --aligner mafft \
        --cores 12 \
        --incomplete-matrix

Alignment stats
^^^^^^^^^^^^^^^

Once you have alignments, it's nice to get a general sense of their length and
composition.  You can quickly (with a multicore machine) summarize thousands of
alignments by running:

.. code-block:: bash

    phyluce_align_get_align_summary_data \
        --alignments /path/to/uce/taxon-set1/mafft-nexus/ \
        --cores 12

This will produce output that looks similar to::

    2014-04-24 20:11:18,208 - get_align_summary_data - INFO - ================ Starting get_align_summary_data ================
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Version: git 7aec8f1
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --alignments: /path/to/uce/taxon-set1/mafft-nexus/
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --cores: 12
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --input_format: nexus
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --log_path: /path/to/uce/taxon-set1/log
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2014-04-24 20:11:18,209 - get_align_summary_data - INFO - Argument --verbosity: INFO
    2014-04-24 20:11:18,210 - get_align_summary_data - INFO - Getting alignment files
    2014-04-24 20:11:18,253 - get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2014-04-24 20:11:20,573 - get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] loci:    1,104
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] length:  752,617
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] mean:    681.72
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] 95% CI:  13.03
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] min:     169
    2014-04-24 20:11:20,574 - get_align_summary_data - INFO - [Alignments] max:     4,520
    2014-04-24 20:11:20,576 - get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2014-04-24 20:11:20,576 - get_align_summary_data - INFO - [Taxa] mean:          24.29
    2014-04-24 20:11:20,576 - get_align_summary_data - INFO - [Taxa] 95% CI:        0.26
    2014-04-24 20:11:20,576 - get_align_summary_data - INFO - [Taxa] min:           3
    2014-04-24 20:11:20,576 - get_align_summary_data - INFO - [Taxa] max:           27
    2014-04-24 20:11:20,577 - get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2014-04-24 20:11:20,577 - get_align_summary_data - INFO - [Missing] mean:       7.97
    2014-04-24 20:11:20,577 - get_align_summary_data - INFO - [Missing] 95% CI:     0.16
    2014-04-24 20:11:20,578 - get_align_summary_data - INFO - [Missing] min:        0.44
    2014-04-24 20:11:20,578 - get_align_summary_data - INFO - [Missing] max:        19.71
    2014-04-24 20:11:20,592 - get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2014-04-24 20:11:20,592 - get_align_summary_data - INFO - [All characters]      18,541,550
    2014-04-24 20:11:20,592 - get_align_summary_data - INFO - [Nucleotides]         14,713,956
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 50%]          1048 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 55%]          1044 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 60%]          1035 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 65%]          1027 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 70%]          1024 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 75%]          1010 alignments
    2014-04-24 20:11:20,594 - get_align_summary_data - INFO - [Matrix 80%]          998 alignments
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Matrix 85%]          994 alignments
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Matrix 90%]          906 alignments
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Matrix 95%]          794 alignments
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Characters] '-' is present 2,301,454 times
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Characters] '?' is present 1,526,140 times
    2014-04-24 20:11:20,595 - get_align_summary_data - INFO - [Characters] 'A' is present 4,092,085 times
    2014-04-24 20:11:20,596 - get_align_summary_data - INFO - [Characters] 'C' is present 3,267,550 times
    2014-04-24 20:11:20,596 - get_align_summary_data - INFO - [Characters] 'G' is present 3,286,742 times
    2014-04-24 20:11:20,596 - get_align_summary_data - INFO - [Characters] 'T' is present 4,067,579 times
    2014-04-24 20:11:20,596 - get_align_summary_data - INFO - ================ Completed get_align_summary_data ===============

.. note::  The alignment summary stats give you some idea of data matrix
    composition at varying levels of completeness in the ``Data matrix completeness
    summary`` section.

Locus name removal
^^^^^^^^^^^^^^^^^^

For historical reasons, and also for users to ensure that the sequence data
aligned together are from the same loci, each sequence line in the alignment
file output by ``seqcap_align`` contains the ``genus_species1`` designator,
but the ``genus_species1`` designator is also prepended with the locus name
(e.g. ``uce-1005_genus_species1``).  We need to remove these if we plan to
concatenate the loci (:ref:`raxml-concat`).  More generally, it is a good idea
to remove locus names from sequence lines before running any analyses. To do
this, run:

 .. code-block:: bash

    phyluce_align_remove_locus_name_from_files \
        --alignments /path/to/uce/taxon-set1/mafft-nexus/ \
        --output /path/to/uce/taxon-set1/mafft-nexus-clean/ \
        --cores 12

.. _finalize-matrix:

Finalize matrix completeness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After checking the resulting alignment summary stats and checking your
alignments for quality, you will generally want to cull the data set to reach
your desired level of completeness. That is easily done by running the
following, while inputting the set of alignments just generated using:

.. code-block:: bash

    # the integer following --taxa is the number of TOTAL taxa
    phyluce_align_get_only_loci_with_min_taxa \
        --alignments /path/to/uce/taxon-set1/mafft-nexus-clean/ \
        --taxa 34 \
        --percent 0.75 \
        --output /path/to/uce/taxon-set1/mafft-nexus-min-25-taxa/ \
        --cores 12

.. attention:: This program computes the ``floor(taxa * percent)`` and uses the
    resulting number to determine the ``min(taxa)`` allowed in an alignment of
    ``--percent`` completeness.

This will produce output that looks similar to::

    2014-04-24 20:12:33,386 - get_only_loci_with_min_taxa - INFO - ============== Starting get_only_loci_with_min_taxa =============
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Version: git 7aec8f1
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Argument --alignments: /path/to/uce/taxon-set1/mafft-nexus
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Argument --cores: 12
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Argument --input_format: nexus
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Argument --log_path: None
    2014-04-24 20:12:33,387 - get_only_loci_with_min_taxa - INFO - Argument --output: /path/to/uce/taxon-set1/mafft-nexus-min-25-taxa
    2014-04-24 20:12:33,388 - get_only_loci_with_min_taxa - INFO - Argument --percent: 0.75
    2014-04-24 20:12:33,388 - get_only_loci_with_min_taxa - INFO - Argument --taxa: 27
    2014-04-24 20:12:33,388 - get_only_loci_with_min_taxa - INFO - Argument --verbosity: INFO
    2014-04-24 20:12:33,388 - get_only_loci_with_min_taxa - INFO - Getting alignment files
    2014-04-24 20:12:35,293 - get_only_loci_with_min_taxa - INFO - Copied 1010 alignments of 1104 total containing >= 0.75 proportion of taxa (n = 20)
    2014-04-24 20:12:35,294 - get_only_loci_with_min_taxa - INFO - ============= Completed get_only_loci_with_min_taxa =============

.. _missing data:

Add missing data designators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Sometimes**, depending on how you will handle the alignments, you will 
need to add missing data designators for taxa missing from
each alignment of a given locus. This will basically allow you to generate
concatenated data sets and it may reduce error messages from other programs
about files having unequal numbers of taxa. To do this, run:

.. code-block:: bash

    phyluce_align_add_missing_data_designators \
        --alignments /path/to/uce/taxon-set1/mafft-nexus-min-25-taxa \
        --output /path/to/uce/taxon-set1/mafft-nexus-min-25-taxa \
        --match-count-output /path/to/uce/taxon-set/dataset1.conf \
        --incomplete-matrix /path/to/uce/taxon-set1/dataset1.incomplete \
        --log-path log \
        --cores 12

.. note:: Here, we're inputting the ``--match-count-output`` and the
    ``--incomplete-matrix`` we created earlier in the :ref:`incomplete-matrix` and
    :ref:`extracting-fasta` sections.


Operations on alignments
========================

Many workflows for phylogenetics simply involve converting one alignment format
to another or changing something about the contents of a given alignment. We
use many of these manipulations in the next section (see :ref:`data-analysis`),
as well.

Converting one alignment format to another
------------------------------------------

To convert one alignment type (e.g., nexus) to another (e.g., fasta), we have a
relative simple bit of code to achieve that process. You can greatly speed this
processing step up on a multicore machine with the ``--cores`` option:

.. code-block:: bash

    phyluce_align_convert_one_align_to_another \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --output /path/to/uce/taxon-set1/mafft-fasta \
        --input-format nexus \
        --output-format fasta \
        --cores 8 \
        --log-path log

You can convert from/to:

#. fasta
#. nexus
#. phylip
#. phylip-relaxed (probably the one you want)
#. clustal
#. emboss
#. stockholm


Shortening taxon names
----------------------

You can shorten taxon names (e.g. for use with strict phylip) by modifying the
above command slightly to add ``--shorten-names``:

.. code-block:: bash

    phyluce_align_convert_one_align_to_another \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --output /path/to/uce/taxon-set1/mafft-fasta-shortnames \
        --input-format nexus \
        --output-format fasta \
        --cores 8 \
        --shorten-names \
        --log-path log


Excluding loci or taxa
----------------------

You may want to exclude loci less than a certain length or having fewer than
a particular number of taxa, or only containing certain taxa.  You can
accomplish that using:

.. code-block:: bash

    phyluce_align_filter_alignments \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --output /path/to/a/new/directory \
        --input-format nexus \
        --containing-data-for genus_species1 genus_species2 \
        --min-length 100 \
        --min-taxa 5 \
        --log-path log

This will filter alignments that do not contain the taxa requested, those
alignments shorter than 100 bp, and those alignments having fewer than 5 taxa
(taxa with only missing data are not counted).


Extracting taxon data from alignments
-------------------------------------

Sometimes you may have alignments from which you want to extract data from a
given taxon, format the alignment string as fasta, and do something with the
fasta results:

.. code-block:: bash

    phyluce_align_extract_taxon_fasta_from_alignments \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --taxon genus_species1 \
        --output /path/to/output/file.fasta


.. _data-analysis:

Preparing concatenated alignment data for analysis
==================================================

Formatting data for analysis generally involves slight differences from the
steps described above.  There are several application-specific programs in
phyluce_.

.. _raxml-concat:

RAxML
-----

For RAxML, you need a concatenated phylip file.  This is pretty easily created
if you have an input directory of nexus alignments.  To create a concatenated
phylip file from many input alignments, run:

.. code-block:: bash

    phyluce_align_concanatenate_alignments \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --output /path/to/uce/taxon-set1/mafft-raxml \
        --phylip

This will output a concatenated file named ``mafft-raxml.phylip`` in
``/path/to/uce/taxon-set1/mafft-raxml``.  It will also include a charset file, 
``mafft-raxml.charsets``.

.. _mrbayes:

MrBayes (Nexus format)
----------------------

You can create a Nexus-formatted file for programs like MrBayes and Paup_ with:

.. code-block:: bash

    phyluce_align_concanatenate_alignments \
        --alignments /path/to/uce/taxon-set1/mafft-nexus \
        --output /path/to/uce/taxon-set1/mafft-raxml \
        --nexus

The charsets will be included in the Nexus file.