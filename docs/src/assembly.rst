#####################
Read assembly
#####################

Once your reads are clean, you're ready to assemble.  A the moment, we use 
`velvet`_ for assembly, although I'm considering a move to `ABySS`_.  Most of
this process is automated, but you need to get a good idea of the kmer ranges
within which you will be working, prior to starting the assembly process
with the automated code.

Assuming your directory structure looks like this::

    uce-clean/
        genus_species/
            interleaved-adapter-quality-trimmed/
                genus_species-READ1and2-interleaved.fastq.gz
                genus_species-READ-singleton.fastq.gz
                
then (after removing the line continuations (`\`) below), run the following,
which assumes you have an 8-core computer (`-t 8`), a sizable amount of RAM
(>= 24 GB), and you want to try kmer values between 63 and 75:

.. code-block:: bash

    cd uce-clean/genus_species/
    mkdir assembly
    cd assembly
    
    VelvetOptimiser -s 63 -e 75 -t 8 -a -c ncon \
    -f '-fastq.gz \
    -shortPaired ../interleaved-adapter-quality-trimmed/genus_species-READ1and2-interleaved.fastq.gz \
    -short ../interleaved-adapter-quality-trimmed/genus_species-READ-singleton.fastq.gz'
    
VelvetOptimiser will assemble the read data in the interleaved-adapter-quality-trimmed
directory and when it is finished, will output the optimal kmer size for the
assembly to stdout.

You want to ensure that the optimal kmer size is **inside** of the range that
you passed to the program.  If it is, conduct several additional assemblies
with other data to ensure that this range works well across all of your samples.

Once you have a good idea of the range size for the data, you may wish to delete
the test `assembly` directories.  You can also keep them, but you will need to
manually symlink over the contigs folder (see below).  Then, you can
run the automated script to assemble your data:

.. code-block:: bash

    python ~/git/phyluce/bin/assembly/assemblo.py uce-clean 63 75 7
    
This will run VelvetOptimiser across your data in the `uce-clean` directory,
using a range of kmer values from 63-75 and 9 compute cores.  You may need to
adjust the number of compute cores for your system and/or amount of RAM.

`assemblo.py` will output the results of each assembly to `stdout`, giving the
file name, optimal kmer size, resulting n50, and an indication of whether the
optimal kmer size was reached during assembly.

If the optimal kmer size is not reached during assembly, you may wish to
re-assembly problematic contigs by hand with different parameters (a smaller or
larger range).

If you want to exclude certain directories nested within `uce-clean`, for
instance, you can run:

.. code-block:: bash

    python ~/git/phyluce/bin/assembly/assemblo.py uce-clean 63 75 7 \
    --exclude genus_species1 genus_species2
    
Within the `uce-clean` directory, `assemblo.py` will create a directory named
`contigs` that contain symlinks to the resulting assembly data for each taxon.

If you assembled any data by hand that are in `assembly` directories, you will
also want to symlink those results in the `contigs` directory.


Determining Success or Failure of an Assembly/Enrichment
********************************************************

This is generally a little bit tricky without having some practice.  Generally
speaking, if you cannot reasonably easily find an optimal kmer value for a
given assembly (the "optimal" kmer is always at the top or bottom of the range),
then the assembly (and data) are likely not as good as they should be.  The
causes of this are too many to describe here, but include low-coverage, poor
enrichment, bad libraries, contamination, etc.

You should generally see that the total number of contigs assembled is within
1-3X the number of contigs that you targeted with your enrichment probe set.
Ideally, you also want the assembled contigs to have an n50 > 300-400 bp.  You
typically **do not** want to see only very large (>10 Kbp) contigs assembled or
only very small (<300 bp) contig assembled.

Again, these results depend on a variety of factors including starting DNA
quality, library size, the efficiency of the enrichment, the number of
post-enrichment PCR cycles you used, the amount of sequence data collected for a
given library, etc., etc., etc.

.. _velvet: http://www.ebi.ac.uk/~zerbino/velvet/
.. _ABySS: http://www.bcgsc.ca/platform/bioinfo/software/abyss