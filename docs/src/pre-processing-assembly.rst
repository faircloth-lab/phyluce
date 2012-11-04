.. include:: global.rst

#####################
Read Assembly
#####################

Once your reads are clean, you're ready to assemble. At the moment, we use
velvet_ and VelvetOptimiser for assembly, although I'm considering adding
similar routines for ABySS_.

Most of this process is automated using code within phyluce_. But, before you
runt the automated routines, you need to get a good idea of the kmer ranges
within which you will assembling, prior to starting the assembly process with
the automated code.

In the following, I assume your directory structure looks like this (which
results from `illumiprocessor`_)::

    uce-clean/
        genus_species1/
            interleaved-adapter-quality-trimmed/
                genus_species-READ1and2-interleaved.fastq.gz
                genus_species-READ-singleton.fastq.gz
                
Finding a reasonable kmer range
*******************************

For the automated procedure, we need to have a reasonably good idea of the kmer
range that we want to plug into the automated process. After several rounds of
this process, you generally get a good idea of what works. However, I still
empirically check this range against a couple of test assemblies before running
many assemblies with the automation code.
                
To determine this range, you want to assemble a test batch of reads using a
range of kmer values. I usually do this for 4-5 taxa in the data set I just
cleaned using kmer values in the range of 63-75. Depending on the so-called
"optimal", kmer value reported by VelvetOptimiser, I'll adjust the upper and
lower limits of this range, hopefully settling on a range of values that will
work across the remaining batches of reads during the automated process.

To assemble a test batch of reads, you want to run the following with
VelvetOptimiser. The commands below assumes you have an 8-core computer (`-t
8`) and a sizable amount of RAM (>= 24 GB). If you have fewer cores or less RAM,
you'll need to adjust the `-t` parameter.  You probably want to initially try
a range of kmer values between 63 and 75.

So, pick a taxon to test with.  You probably want to start with something
having roughly the mean number of reads across all taxa you included in a run.

First, you want to enter a particular taxons `uce-clean` directory of reads and
create an `assembly` directory to hold the assembly results:

.. code-block:: bash

    cd uce-clean/genus_species1/
    mkdir assembly
    cd assembly
    
Now that you're in the directory, you want to run VelvetOptimiser against the
interleaved read data that are one directory up:

.. code-block:: bash

    VelvetOptimiser -s 63 -e 75 -t 8 -a -c ncon \
    -f '-fastq.gz \
    -shortPaired ../interleaved-adapter-quality-trimmed/genus_species-READ1and2-interleaved.fastq.gz \
    -short ../interleaved-adapter-quality-trimmed/genus_species-READ-singleton.fastq.gz'

VelvetOptimiser will assemble the read data in the
interleaved-adapter-quality-trimmed directory and when it is finished, will
output the optimal kmer size for the assembly to stdout.  Write down and/or
remember this number (probably best to write it down).

You also want to ensure that the optimal kmer size is **inside** of the range
that you passed to the program. If it is, conduct several additional assemblies
with other data to ensure that this range works well across all of your samples.

If the "optimal" kmer size is at the **edge** of the range you input, then you
should adjust the range of kmer values to try, and re-assemble the data to see
if the next "optimal" kmer value is within the range you input.

Once you have a good idea of the range size for the data, you may wish to delete
the test `assembly` directories.  You can also keep them, but you will need to
manually symlink over the contigs folder (see below).

Running assemblo.py
*******************

Once you have a good idea of the range that will work for your data,
you can run the automated script to assemble your data:

.. code-block:: bash

    python ~/git/phyluce/bin/assembly/assemblo.py uce-clean 63 75 7
    
This will run VelvetOptimiser across your data in the `uce-clean` directory,
using a range of kmer values from 63-75 and 7 compute cores. As before, you may
need to adjust the number of compute cores for your system and/or amount of RAM.

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

Linking contig assemblies into the `contigs` directory
******************************************************

`assemblo.py` will assemble your contigs, and symlink them into a `contigs`
directory.  These symlinks allow you to reasonably keep track of which assembly
is linked to a particular taxon, maintaining continuity in your data.

If you have manually assembled any contigs outside of assembly, you will also 
need to manually link those into the contigs folder with something like:

.. code-block:: bash

    ln -s uce-clean/genus_species1/genus_species1/assembly/auto_data_75/contigs.fa \
    contigs/genus_species1.contigs.fasta
    
The other downside of symlinks is that they tend not to "travel" very well.
Meaning, if you want to move the contigs folder somewhere else, this will often
break the symlinks created by `assemblo.py`.  This is often painful to fix by hand,
so there's a bit of code in `phyluce_` to help you do that.  If you have reads in::

    /path/to/my/uce-clean/contigs
    
and you want to create a new `contigs` folder elsewhere, you want to replace the path
in the original symlink with the path in the new symlink.  So, you can run::

    python ~/Git/brant/phyluce/bin/share/replace_many_links.py \
        /path/to/my/uce-clean/contigs \
        /path/to/my/uce-clean/clean/ \
        /new/path/to/my/uce-clean/location/
        name-of-new-folder

Determining success or failure of an assembly/enrichment
********************************************************

This is a little bit tricky without having some previous experience. Generally
speaking, if you cannot reasonably easily find an optimal kmer value for a
given assembly (the "optimal" kmer is **always** at the top or bottom of the
range), then the assembly (and data) are likely not as good as they should be.
The causes of this are too many to describe here, but include low-coverage,
poor enrichment, bad libraries, contamination, etc.

You should generally see that the total number of contigs assembled is within
1-3X the number of contigs that you targeted with your enrichment probe set.
Ideally, you also want the assembled contigs to have an n50 > 300-400 bp.  You
typically **do not** want to see only very large (>10 Kbp) contigs assembled or
only very small (< 300 bp) contigs assembled.

Again, these results depend on a variety of factors including starting DNA
quality, library size, the efficiency of the enrichment, the number of
post-enrichment PCR cycles you used, the amount of sequence data collected for
a given library, etc., etc., etc.

.. _velvet: http://www.ebi.ac.uk/~zerbino/velvet/
.. _ABySS: http://www.bcgsc.ca/platform/bioinfo/software/abyss