.. include:: global.rst

UCE Processing for Phylogenomics
================================

The workflow described below is meant to outline the process needed for
analyzing UCE in phylogenetic contexts - meaning that you are interested in
addressing questions at the species level or deeper.

Probe sets
**********

Depending on the taxa you are targeting and the probe set that you are using,
you will need to adjust the probe set name to reference the fasta file of those
probes you used. Below, I have used the `uce-5k-probes.fasta` probe set.
However, the probe set that use could be one of:

#. uce-2k-probes.fasta (tetrapods/birds/mammals)
#. uce-5k-probes.fasta (tetrapods/birds/mammals)

.. _outgroup-data:

Outgroup data and probe set downloads
*************************************

We have started to standardize and provide prepared sets of UCE probes and
outgroup data. The outgroup data are sliced from available genome sequences,
and the probe sets and outgroup data are version controlled for a particular
set of taxa (e.g., tetrapods, fish). If you would like to use these in your own
analyses, to extend your current data set or to provide an outgroup, you can
download them from uce-probe-sets_

.. _contigs-matching:

Indentifying contigs matching UCE loci
**************************************

After assembly, we have generated contigs from raw reads.  These contigs reside
in the `contigs` resulting from assembly.  During the next part of the process,
we need to determine which of the assembled contigs match UCE loci and which do
not.  We also need to remove any contigs that appear to be duplicates as a
result of assembly/other problems **or** a duplication event(s).

The first thing to do is to make sure that our probe set does not contain any
duplicates.  So, you probably want to align the file of probe sequences to
itself (if you're using a probe-set from github, the this file should be
included):

.. code-block:: bash

    python phyluce/bin/share/easy_lastz.py \
        --target uce-5k-probes.fasta \
        --query uce-5k-probes.fasta \
        --identity 85 \
        --output uce-5k-probes.fasta.toself.lastz
        
Now, what we need to do is to align our probes to our contigs.  First, you want
to make a directory to hold our output:
        
.. code-block:: bash

    mkdir /path/to/output/lastz
    
Then, we want to align the contigs we assembled to the UCE loci represented
in the uce-5k-probes.fasta file.  Since we're using "new-style" (standardized)
probe names, we want to include the flags::

    --regex "_p[1-9]+$" --repl ""
    
Which we use to strip the probe numbers off of particular loci in the
`uce-5k-probes.fasta` file (stripping off the probe numbers allows us to
merge all probes down to a single locus).  Note, too, that we're passing the
`uce-5k-probes.fasta.toself.lastz` to the code so that we can also exclude
any UCE loci whose probes happen to overlap themselves:

.. code-block:: bash


    python phyluce/bin/assembly/match_contigs_to_probes.py \
        /path/to/velvet/assembly/contigs/ \
        /path/to/uce-5k-probes.fasta \
        /path/to/output/lastz \
        --regex "_p[1-9]+$" --repl "" \
        --dupefile uce-5k-probes.fasta.toself.lastz
        
When you run this code, you will see output similar to::

    genus_species1: 1031 (70.14%) uniques of 1470 contigs, 0 dupe probe matches, 48 UCE probes matching multiple contigs, 117 contigs matching multiple UCE probes
    genus_species2: 420 (68.52%) uniques of 613 contigs, 0 dupe probe matches, 30 UCE probes matching multiple contigs, 19 contigs matching multiple UCE probes
    genus_species3: 1071 (63.15%) uniques of 1696 contigs, 0 dupe probe matches, 69 UCE probes matching multiple contigs, 101 contigs matching multiple UCE probes

Now, what this program does is to use `lastz_` to align all probes to the
contigs. It basically ignores those contigs that don't match probes (no
matches) and screens the results to ensure that, of the matches, only one
contig matches probes from one UCE locus and that only probes from one UCE
locus match one contig. **Everything outside of these parameters is dropped**.

The resulting files will be in the::

    /path/to/output/lastz
    
directory. You'll see that this directory contains species-specific `lastz_`
files as well as an sqlite database::

    /path/to/output/lastz
        genus_species1.contigs.lastz
        genus_species2.contigs.lastz
        genus_species3.contigs.lastz
        probe.matches.sqlite
        
The `*.lastz` files are basically for reference and individual review.  The
really important data are actually summarized in the::

    probe.matches.sqlite
    
database.  It's probably a good idea to have some knowledge of how this database
is structured, since it's basically what makes the next few steps work.  So, I'll
spend some time describing the structure and contents.

The probe.matches.sqlite database
*********************************

`probe.matches.sqlite` is a relational database that summarizes all **valid**
matches of contigs to UCE loci across the set of taxa that you fed it. The
database is created by and for sqlite_, which is a very handy, portable SQL
database. For more info on SQL and SQLITE, see this `sqlite-tutorial`_. I'll briefly cover the
database contents and use below.

First, to take a look at the contents of the database run:

.. code-block:: bash

    sqlite3 probe.matches.sqlite
    
You'll now see something like::

    SQLite version 3.7.3
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite>
    
It's often easier to change some defaults for better viewing, so at the prompt, 
past in the following (for more info on sqlite_ "dot" commands, you can type
`.help`)::

    sqlite> .mode columns
    sqlite> .headers on
    sqlite> .nullvalue .
    
Now that that's done, let's see what tables the database contains::

    sqlite> .tables
    match_map  matches
    
This tells us there's two tables in the database, named `match_map` and
`matches`.  We'll look at `matches`, first.  To get some data out of `matches`,
run (the use of uppercase is convention for SQL, but not required):

The `matches` table
-------------------

Let's take a look at the contents of the `matches` table.  Once you've started
the sqlite interface, run:

.. code-block:: sql

    sqlite> SELECT * FROM matches LIMIT 10;
    
This query select all rows (`SELECT *`) from the `matches` table (`FROM
matches`) and limits the number of returned rows to 10 (`LIMIT 10`). This will
output data that look something like::

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
from `genus_species1`, 3 of 10 UCE loci in the list from `genus_species2`, and
2 of 10 UCE loci from `genus_species3`. The locus name is given in the `uce
column`.  Remember that we've limited the results to 10 rows for the sake of
making the results easy to view.

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

Basically, the `matches` table and this query are what we run to generate
**complete** (only loci enriched in all taxa) and **incomplete** (all loci
enriched from all taxa) datasets (see :ref:`locus-counts`).

The `match_map` table
---------------------

The `match_map` table shows us which species-specific, velvet-assembled contigs
match which UCE loci. Because velvet assigns an arbitrary designator to each
assembled contig, we need to map these arbitrary designators (which differ for
each taxon) to the UCE locus to which it corresponds. Because velvet contigs
are not in any particular orientation (i.e., they may be 5' - 3' or 3' - 5'),
we also need to determine the orientation of all contigs relative to the source
probe file.

Let's take a quick look:

.. code-block:: sql

    SELECT * FROM match_map LIMIT 10;

This query is similar to the one that we ran against `matches` and returns the
first 10 rows of the `match_map` table::

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
    
As stated above, these results show the "hits" of velvet-assembled contigs to
particular UCE loci. So, if we were to open the `genus_species1.contigs.fasta`
symlink (which connects to the assembly) in the `contigs` folder, the contig
named `node_233` corresponds to UCE locus `uce-500`.

Additionally, each entry in the rows also provides the orientation for
particular contigs `(-)` or `(+)`. This orientation is relative to the
orientation of the UCE probes/locus in the source genome (e.g., chicken for
tetrapod probes).

We use this table to generate a FASTA file of UCE loci for alignment (see
:ref:`locus-counts`), after we've identified the loci we want in a particular
data set. The code for this step also uses the associated orientation data to
ensure that all the sequence data have the same orientation prior to alignment
(some aligners will force alignment of all reads using the given orientation
rather than also trying the reverse complement and picking the better alignment
of the two).

.. _locus-counts:

Determining locus counts and generating a taxon-set
***************************************************

Now that we know the taxa for which we've enriched UCE loci and which 
contigs we've assembled match which UCE loci, we're ready to generate some data
sets.  The data set generation process is pretty flexible - you can select which
taxa you would like to group together for an analysis, you can generate complete
and incomplete data matrices, and you can also include additional data from the
provided outgroup files and data (see :ref:`outgroup-data`) or previous runs.
We'll start simple.

Complete matrix data set
------------------------

First, we'll generate a data set from only the current UCE enrichments,
and it will be complete - meaning that we will not include loci where certain
taxa have no data (either the locus was not enriched for that taxon or removed
during the filtering process for duplicate loci).

The first step of generating a data set is to identify those loci present in the
taxa with which we're working.  First, you need to create a configuration (text)
file denoting the taxa we want in the data set.  It should look like this::

    [dataset1]
    genus_species1
    genus_species2
    genus_species3
    
Let's assume you name this file `datasets.conf`.  Now, you want to run the
following against this file, along with several other files we've created
previously::

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/your/datasets.conf \
        'dataset1' \
        --output /path/to/some/output-file/dataset1.conf
        
This will basically run a query against the database, and pull out those loci
for those taxa in the `datasets.conf` file having UCE contigs.  The output will
look something like::

    Shared UCEs: 500

    genus_species1:108
    genus_species2:93
    genus_species3:71
    
This means that 500 loci are shared amongst the 3 taxa in `datasets.conf`.  We
might have had more, but `genus_species1` caused us to drop 108 loci,
`genus_species2` caused us to drop 93 loci, and `genus_species3` caused us to
drop 71 loci.

Now, you might think that increasing the locus count is simply a matter of
removing `genus_species1` from the list of taxa.  This is not strictly true,
however, given the vagaries of hits and misses among taxa. `get_match_counts.py`
has several other options to help you determine which taxa may be causing
problems, but picking the best combination of taxa to give you the highest
number of loci is a somewhat hard optimization problem.

If you want to generate/evalaute additional data sets with different taxa, you
can simply append that list to the `datasets.conf` file like so::

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

and then run `get_match_counts.py` against this new section::

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/your/datasets.conf \
        'dataset2' \
        --output /path/to/some/output-file/dataset2.conf

Incomplete data matrix
----------------------

You may not always want a complete data matrix or generating a complete matrix
drops too many loci for your tastes.  That's cool.  You can generate an
incomplete dataset like so:

.. code-block:: bash

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/your/datasets.conf \
        'dataset1' \
        --output /path/to/some/output-file/dataset1-incomplete.conf
        --incomplete-matrix

This will generate a dataset that includes any loci enriched across the taxa
in the `datasets.conf` file.  This will also include a file named
`dataset1-incomplete.notstrict` that contains those loci enriched for
each taxon.  We'll need that in a minute (see :ref:`extracting-fasta`)

Incorporating outgroup/other data
---------------------------------

You may want to include outgroup data from another source into your datasets.
This can be from the pre-processed outgroup data files (see
:ref:`outgroup-data`), but it doesn't need to be these outgroup data. These
additional data can also be contigs previously assembled from a different set
of taxa.

The first step of this process is to setup your `datasets.conf` slightly
differently - by indicating these external data with asterisks::

    [dataset3]
    genus_species1
    genus_species2
    genus_species3
    genus_species4*
    genus_species5*
    
Then, you need to pass `get_match_counts.py` the location of the
`probe.matches.sqlite` database previously generated as described in
:ref:`contigs-matching` or downloaded as part of :ref:`outgroup-data`::

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/your/datasets.conf \
        'dataset3' \
        --extend /path/to/some/other/probe.matches.sqlite \
        --output /path/to/some/output-file/dataset3-with-external.conf
        
To keep all this extension from getting too terribly crazy, I've limited the
ability to include external data to essentially a single set.  If you have lots
of data from many different enrichments, you'll need to generate a `contigs`
folder containing all these various assemblies (or symlinks to them), then
align the probes to these data (see :ref:`contigs-matching`).  Once you do that,
you can extend your current data set with all of these other data.

.. _extracting-fasta:

Extracting relevant FASTA data
******************************

After selecting the set of loci in which you're interested, you need to
generate a FASTA file containing the reads from each species-specific contig
that corresponds to a locus in the set.  This is reasonable easy.

Complete data matrix
--------------------

To generate a FASTA file, we're passing several previously used paths and 
the name of the output file from `get_match_counts.py` on the third line below
(`/path/to/some/output-file/dataset1.conf`)::

    python phyluce/bin/assembly/get_fastas_from_match_counts.py \
        /path/to/velvet/assembly/contigs/ \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/some/output-file/dataset1.conf \
        --output /path/to/some/output.fasta

Incomplete data matrix
----------------------

To generate a FASTA file, we're passing several previously used paths plus the
name of the output file from `get_match_counts.py` on the fourth line **AND**
the name of the `*.notstrict` file on the fifth line::

    python phyluce/bin/assembly/get_fastas_from_match_counts.py \
        /path/to/velvet/assembly/contigs/ \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/some/output-file/dataset1-incomplete.conf \
        --incomplete-matrix /path/to/some/output-file/dataset1-incomplete.notstrict
        --output /path/to/some/output.fasta

Incorporating outgroup/other data
---------------------------------

Because we're incorporating external data, we need to pass the name of the
external database, as before, as well as the name of the external `contigs`
directory::

    python phyluce/bin/assembly/get_fastas_from_match_counts.py \
        /path/to/velvet/assembly/contigs/ \
        /path/to/output/lastz/probe.matches.sqlite \
        /path/to/some/output-file/dataset3-with-external.conf \
        --extend-db /path/to/some/other/probe.matches.sqlite \
        --extend-dir /path/to/some/other/contigs/ \
        --output /path/to/some/output.fasta

Aligning and trimming FASTA data
********************************

With all of that out of the way, things get much easier to deal with.  We
basically need to align our data across loci, and we're largely ready to go.
The remaining operations we can run on the data are format-conversions, QC steps
or any number of other fun things.

Aligning this much data is reasonably computationally intensive - so this
alignment step goes fastest if you have a multicore machine.  You also have
several alignment options available, although I would suggest sticking with
MAFFT.

First, make a folder for the alignment output::

    mkdir /path/to/alignment/output
    
Complete data matrix
--------------------

The second line is the fasta created above (see :ref:`extracting-fasta`), the
second line is the path to the output, the third line gives the number of taxa
in the alignment, `--aligner mafft` determines the alignment program, and
`--cores 8` denoted the number of cores to use for this step::

    python phyluce/bin/align/seqcap_align_2.py \
        /path/to/some/output.fasta \
        /path/to/alignment/output \
        3 \
        --aligner mafft \
        --cores 8

Incomplete data matrix
----------------------
 
The only difference for an alignment of incomplete data is that we also pass
the `--notstrict` flag, which tells the code to expect that some loci will not
have data for all taxa::

    python phyluce/bin/align/seqcap_align_2.py \
        /path/to/some/output.fasta \
        /path/to/alignment/output \
        3 \
        --aligner mafft \
        --incomplete-matrix \
        --cores 8
        
After checking the resulting alignment QC (see :ref:`alignment-QC`), you will
generally need to add in missing data designators for taxa missing from the
alignment of a given locus. This will basically allow you to generate
concatenated data sets and it may reduce error messages from other programs
about files having unequal numbers of taxa. To do this, you need to run::

    python phyluce/bin/align/add_missing_data_designators.py \
        /path/to/alignment/output  \
        /path/to/alignment/output-with-missing-data/ \
        /path/to/some/output-file/dataset3-with-external.conf \
        /path/to/some/output-file/dataset3-with-external.notstrict
        
Alignment trimming
------------------

The alignment code "trims" alignments by default.  Basically, this means that
it removes ragged 5' and 3' edges from a given alignment.  However, you may not
want to run the trimming and just deal with the raw alignments output by
mafft/muscle/dialign. No problem, you run `seqcap_align_2.py` just as above, but
you add the `--notrim` option::

    python phyluce/bin/align/seqcap_align_2.py \
        /path/to/some/output.fasta \
        /path/to/alignment/output \
        3 \
        --aligner mafft \
        --cores 8 \
        --notrim
        
Saté alignment
--------------

There is also an option to run Saté alignments instead of the default code.  For
the moment, this code lives in `mpi_sate.py` and you can run it locally with
something like::

    python phyluce/bin/align/mpi_sate.py \
        /path/to/some/output.fasta \
        /path/to/alignment/output \
        3 \
        /path/to/sate \
        /path/to/sate.cfg \
        --parallelism multiprocessing \
        --cores 8
        
This code will also run on MPI enabled machines, but that is generally 
beyond the scope of this HOWTO.

Alignment trimming only
-----------------------

If you have untrimmed (ragged) alignments that you would like to trim with the
phyluce_ trimming procedures, you can also run that::

    python phyluce/bin/align/get_trimmed_alignments_from_untrimmed.py \
        /path/to/alignment/input \
        /path/to/output/for/trimmed/data/ \
        --input-format nexus
        --output-format nexus \
        --multiprocessing

.. _alignment-qc:

Alignment quality control
*************************

There are many ways to QC alignments.  The best way is to do it visually, but
that gets somewhat hard when you have thousands of loci.  There are several
programs in the `phyluce`_ package that help you QC alignments.  You probably
always want to run::

    python ~/git/phyluce/bin/align/get_align_summary_data.py \
        /path/to/alignment/output-renamed \
        --input-format nexus
        
This will output a number of stats that look somewhat like (these examle data
are from an incomplete matrix)::

    uce-1071.nex is < 100 bp long
    uce-720.nex is < 100 bp long

    Lengths
    -----
    Total length(aln)        256066              
    Average length(aln)      318.490049751       
    95 CI length(aln)        10.6004273805       
    Minimum length(aln)      64                  
    Maximum length(aln)      933                 

    Taxa
    -----
    Average(taxa)            11.526119403        
    95 CI(taxa)              0.37345195065       
    min(taxa)                3                   
    max(taxa)                21                  
    Count(taxa:# alns)       {3: 77, 4: 44, 5: 38, 6: 35, 7: 36, 8: 33, 9: 47, 10: 35, 11: 31, 12: 34, 13: 50, 14: 62, 15: 54, 16: 58, 17: 38, 18: 45, 19: 41, 20: 35, 21: 11}

    Base composition
    -----
    Bases                    {'A': 657715, 'C': 533302, '-': 380870, 'T': 667908, 'G': 523342}
    Sum(all)                 2763137             
    Sum(nucleotide only)     2382267             
    Missing data from trim (%)5.41
    
Sometimes, loci will contain bases that are not in the standard set of IUPAC
base code (e.g. "X" or "N").  To identify these loci, you can run::

    python phyluce/bin/align/screen_alignments_for_problems.py \
        /path/to/alignment/output-renamed \
        --input-format nexus

Alignment name cleaning
***********************

So that you can visually check the resulting alignments to make sure the correct
reads for each taxon are included in a given alignment, the `seqcap_align_2.py`
program writes output files that contain the locus name as part of the taxon
name in the output nexus files.

This is likely to change in the near future.  However, in the meantime, you
probably want to remove this designation from the resulting alignment files.
You can easily do this with::

    python phyluce/bin/align/remove_locus_name_from_nexus_lines.py \
        /path/to/alignment/output \
        /path/to/alignment/output-renamed \
        3

The second line gives the path to the output created during alignment, the third
line gives the path to store the cleaned alignments, and the third line gives
the number of taxa in each alignment.

Alignment manipulation
**********************

Many workflows for phylogenetics simply involve converting one alignment format
to another or changing something about the contents of a given alignment. We
use many of these manipulations in the next section (see :ref:`data-analysis`),
as well.

Converting one alignment format to another
------------------------------------------

To convert one alignment type (e.g., nexus) to another (e.g., fasta), we have a
relative simple bit of code to achieve that process. You can also speed this
processing step on a multicore machine with the `--cores` option::

    python phyluce/bin/align/convert_one_align_to_another.py \
        /path/to/input/alignments \
        /path/to/output/alignments \
        --input-format nexus \
        --output-format fasta \
        --cores 8
        
You can convert from/to:

#. fasta
#. nexus
#. phylip
#. clustal
#. emboss
#. stockholm

Shortening taxon names
----------------------

You can shorten taxon names (e.g. for use with strict phylip) by modifying the
above command slightly to add `--shorten-names`::

    python phyluce/bin/align/convert_one_align_to_another.py \
        /path/to/input/alignments \
        /path/to/output/alignments \
        --input-format nexus \
        --output-format fasta \
        --cores 8 \
        --shorten-names

Excluding loci or taxa
----------------------

You may want to exclude loci less than a certain length or having fewer than
a particular number of taxa, or only containing certain taxa.  You can
accomplish that using::

    python phyluce/bin/align/filter_alignments.py \
        /path/to/alignment/output-renamed \
        --input-format nexus \
        --containing-data-for genus_species1 genus_species2 \
        --min-length 100 \
        --min-taxa 5 \
        --output /path/to/a/new/directory
        
This will filter alignments that do not contain the taxa requested, those
alignments shorter than 100 bp, and those alignments having fewer than 5 taxa
(taxa with only missing data are not counted).

Extracting taxon data from alignments
-------------------------------------

Sometimes you may have alignments from which you want to extract data from a
given taxon, format the alignment string as fasta, and do something with the
fasta results::

    python phyluce/bin/align/extract_taxon_data_from_alignments.py \
        /path/to/alignment/ \
        genus_species1 \
        /path/to/output/file.fasta \
        --input-format nexus
        
.. _data-analysis:

Preparing alignment data for analysis
*************************************

Formatting data for analysis generally involves slight differences from the
steps described above.  There are several application-specific programs in
phyluce_.

RAxML
-----

For RAxML, you need a concatenated phylip file.  This is pretty easily created
if you have an input directory of nexus alignments.  First, make an output
directory::

    mkdir raxml
    
Then run::

    python phyluce/bin/align/format_nexus_files_for_raxml.py \
        /path/to/alignment/output-renamed \
        raxml/output-file-name.phylip

.. _strict-phylip:

PHYLIP/CloudForest
------------------

PHYLIP, PhyML, and other programs like CloudForest_ require input files to be in
strict phylip format for analysis.  Converting alignment files to this format
was discussed above, and is simple a matter of (use `--cores` if you have 
a multicore machine as that will greatly speed processing)::

    python phyluce/bin/align/convert_one_align_to_another.py \
        /path/to/input/alignments \
        /path/to/output/alignments \
        --input-format nexus \
        --output-format phylip \
        --shorten-names

MrBayes
--------

MrBayes is a little more challenging to run.  This is largely due to the fact
that we usually estimate the substitution models for all loci, then we partition
loci by substitution model, concatenate the data, and format an appropriate
file to be input to MrBayes.

The tricky part of this process is estimating the locus-specific substitution
models.  Generally speaking, I do this with CloudForest_ now, then I strip the
best-fitting substitution model from the CloudForest_ output, and input that
file to the program that creates a nexus file for MrBayes.

First, estimate the substitution models using cloudforest (this will also give
you genetrees for all loci, as a bonus).  You will need your alignments in 
strict phylip format::

    python cloudforest/cloudforest_mpi.py \
        /path/to/strict/phylip/alignments/ \
        /path/to/store/cloudforest/output/ \
        genetrees \
        $HOME/git/cloudforest/cloudforest/binaries/PhyML3linux64 \
        --parallelism multiprocessing \
        --cores 8
        
In the above, `genetrees` is a keyword that tells CloudForest_ that you mean to
estimate genetrees (instead of bootstraps).  Depending on the size of your
dataset (and computer), this may take some time.  Once this is done::

    python phyluce/bin/genetrees/split_models_from_genetrees.py \
        /path/to/cloudforest/output/genetrees.tre \
        /path/to/output_models.txt
        
Now, you're ready to go with formatting for MrBayes - note that we're inputting
the path of the models file created above (output_models.txt) on line 3::

    python phyluce/bin/align/format_nexus_files_for_mrbayes.py \
        /path/to/input/nexus/ \
        /path/to/output_models.txt \
        /path/to/output/mrbayes.nexus \
        --interleave \
        --unlink

This should create a partitioned data file for you. The partitioning will be by
model, not by locus. Should you want to fully partition by locus (which may
overparamterize), then you can run::

    python phyluce/bin/align/format_nexus_files_for_mrbayes.py \
        /path/to/input/nexus/ \
        /path/to/output_models.txt \
        /path/to/output/mrbayes.nexus \
        --interleave \
        --unlink \
        --fully-partition

CloudForest (genetree/species tree)
-----------------------------------

CloudForest_ is a program written by Nick Crawford and myself that helps you
estimate genetrees and perform bootstrap replicates for very large datasets.
Data input to CloudForest should be in strict phylip format (see
:ref:`strict-phylip`).  First, as above, run genetree analysis on your data (
if you ran this above, you don't need to run it again).  This will estimate
the genetrees for each locus in your dataset, using it's best fitting
substitution model)::

    python cloudforest/cloudforest_mpi.py \
        /path/to/strict/phylip/alignments/ \
        /path/to/store/cloudforest/output/ \
        genetrees \
        $HOME/git/cloudforest/cloudforest/binaries/PhyML3linux64 \
        --parallelism multiprocessing \
        --cores 8
        
The, to generate bootstrap replicates, you can run::

    python cloudforest/cloudforest_mpi.py \
        /path/to/strict/phylip/alignments/ \
        /path/to/store/cloudforest/output/ \
        bootstraps \
        $HOME/git/cloudforest/cloudforest/binaries/PhyML3linux64 \
        --parallelism multiprocessing \
        --cores 8 \
        --bootreps 1000 \
        --genetrees /path/to/store/cloudforest/output/genetrees.tre

**NOTE** that depending on your system, you may need to choose another value
for the path to PhyML::
    
    $HOME/git/cloudforest/cloudforest/binaries/PhyML3linux64
    
RaXML (genetree/species tree)
-----------------------------

We can also use RaXML to genrate gene trees to turn into a species tree. To keep
the taxa names similar to what I run through CloudForest_, I usually input
strict phylip formatted files to these runs (see :ref:`strict-phylip`).  Once
that's done, you can generate genetrees with::

    python phyluce/bin/genetrees/run_raxml_genetrees.py \
        /path/to/strict/phylip/alignments/ \
        /path/to/store/raxml/output/ \
        --outgroup genus_species1 \
        --cores 12 \
        --threads 1

Number of `--cores` is the number of simultaneous trees to estimate, while
`--threads` is the number of threads to use for each tree.  Although somewhat
counterintuitive, I've found that 1 `--thread` per locus and many locis being
processed at once is the fastest route to go.

Once that's finished, you can genrate bootstrap replicates for those same loci::

    python phyluce/bin/genetrees/run_raxml_bootstraps.py \
        /path/to/strict/phylip/alignments/ \
        /path/to/store/raxml/output/ \
        --bootreps 100 \
        --outgroup genus_species1 \
        --cores 12 \
        --threads 1