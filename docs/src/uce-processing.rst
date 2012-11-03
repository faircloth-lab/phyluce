UCE Processing
==============

The workflow described below is meant to outline the process needed for
analyzing UCE in phylogenetic contexts - meaning that you are interested in
addressing questions at the species level or deeper.

Probe sets
**********


Outgroup data and probe set downloads
*************************************

We have started to standardize and provide prepared sets of UCE probes and 
outgroup data.  The outgroup data are sliced from available genome sequences, 
and the probe sets and outgroup data are version controlled for a
particular set of taxa (tetrapods, fish).  If you would like to use these in
your analysis, to extend your current data set or to provide an outgroup,
you can download them from `uce-probe-sets`_

.. _uce-probe-sets: https://github.com/faircloth-lab/uce-probe-sets


Indentifying contigs matching UCE loci
**************************************

.. code-block:: bash

    python phyluce/bin/share/easy_lastz.py \
        --target uce-5k-probes.fasta \
        --query uce-5k-probes.fasta \
        --identity 85 \
        --output uce-5k-probes.fasta.toself.lastz
        
.. code-block:: bash

    python phyluce/bin/assembly/match_contigs_to_probes.py \
        /contigs/from/velvet/assembly/directory/ \
        /path/to/uce-5k-probes.fasta \
        /path/to/probe-id/lastz \
        --regex "_p[1-9]+$" --repl "" \
        --dupefile uce-5k-probes.fasta.toself.lastz
        
You will see output similar to::

    genus_species1: 1031 (70.14%) uniques of 1470 contigs, 0 dupe probe matches, 48 UCE probes matching multiple contigs, 117 contigs matching multiple UCE probes
    genus_species2: 420 (68.52%) uniques of 613 contigs, 0 dupe probe matches, 30 UCE probes matching multiple contigs, 19 contigs matching multiple UCE probes
    genus_species3: 1071 (63.15%) uniques of 1696 contigs, 0 dupe probe matches, 69 UCE probes matching multiple contigs, 101 contigs matching multiple UCE probes


Determining locus counts and generating a taxon-set
***************************************************

Complete data matrix
--------------------

.. code-block:: bash

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/probe-id/lastzprobe.matches.sqlite \
        /path/to/taxa/in/group/name.conf \
        'section header' \
        --output /path/to/output-file/name.conf

Incomplete data matrix
----------------------

.. code-block:: bash

    python phyluce/bin/assembly/get_match_counts.py \
        /path/to/probe-id/lastzprobe.matches.sqlite \
        /path/to/taxa/in/group/name.conf \
        'section header'
        --output /path/to/output-file/name.conf
        --incomplete-matrix


Extracting relevant FASTA data
******************************

Aligning and trimming FASTA data
********************************

Alignment quality control
*************************


Preparing FASTA data for analysis
*********************************









