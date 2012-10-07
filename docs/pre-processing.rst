#####################
Data Pre-processing
#####################

**************
Demultiplexing
**************

Because we are reducing genomes, and beceause we generally reduce genomes to a
fraction of their original size, we can combine DNA libraries from many sources
into a single lane of HiSeq (or even MiSeq) sequencing.  And, because we
are multiplexing these libraries prior to sequencing, we need to "demultiplex"
these libraries after they're sequenced.

Demultiplexing is a tricky process, made trickier by the platform with which you 
are working.  Thankfully, demultiplexing (or "demuxing") Illumina reads is
rather simple. There are several options:

* Casava_ (the Illumina pipeline)
* `fastx-toolkit`_
* bartab_

I have also written splitaake_, which is a program to help demultiplex
Illumina reads.  splitaake_ is particularly useful when you have very many
sequence tags (aka "barcodes") to demultiplex, when your sequence tags are
long, or some combination of the two.

To demultiplex your data with splitaake_, please see the `website`__.


.. _splitaake: https://github.com/faircloth-lab/splitaake/
.. _screen: http://www.gnu.org/software/screen/
.. _tmux: http://tmux.sourceforge.net/
.. _gzip: http://www.gzip.org/
.. _Casava: http://support.illumina.com/sequencing/sequencing_software/casava.ilmn
.. _fastx-toolkit: http://hannonlab.cshl.edu/fastx_toolkit/
.. _bartab: http://www.phyloware.com/Phyloware/XSTK.html
__ splitaake_ 

******************************
Adapter- and quality- trimming
******************************

After you demultiplex your data, you generally end up with a number of files in 
a directory that look like the following, if you used splitaake::

    some-directory-name/
        bfidt-000_R1.fastq.gz
        bfidt-000_R2.fastq.gz
        bfidt-001_R1.fastq.gz
        bfidt-001_R2.fastq.gz
    
or, if you used Casava::

    Project_name/
        Sample_BFIDT-000
            BFIDT-000_AGTATAGCGC_L001_R1_001.fastq.gz
            BFIDT-000_AGTATAGCGC_L001_R2_001.fastq.gz
        Sample_BFIDT-001
            BFIDT-001_TTGTTGGCGG_L001_R1_001.fastq.gz
            BFIDT-001_TTGTTGGCGG_L001_R2_001.fastq.gz

Ideally, what you want to do is to clean these reads of adapter contamination
and trim low-quality bases from each reads (and probably also drop reads
containing "N" (ambiguous) bases.  Then you want to interleave the resulting
data, where read pairs are maintained, and also have an output file of singleton
data, where read pairs are not.

You can do this however you like.  You want a resulting directory structure that
looks like (replace genus_species with your taxon names)::

    genus_species1/
        interleaved-adapter-quality-trimmed/
            genus_species1-READ1and2-interleaved.fastq.gz
            genus_species1-READ-singleton.fastq.gz
    genus_species2/
        interleaved-adapter-quality-trimmed/
            genus_species2-READ1and2-interleaved.fastq.gz
            genus_species2-READ-singleton.fastq.gz
            
The "READ1and2" file should be in interleaved file containig both reads kept,
and the READ file should contain those data where one read of the pair are kept
(AKA singleton reads).


Illumiprocessor
-----------------

I have a package that I use for adapter and quality trimming named 
`illumiprocessor`_.  It generally automates these processes and produces output
if the format we want downstream.

If you used Casava, it will be easiest to place all of the demuliplexed reads
into a single directory.  If you used `splitaake`_, then things should be all
set.

You need to generate a configuration file, that gives details of your reads,
how you want them processed, and what renaming options to use.  This file is
an extension of the INI file used to `splitaake`_ and looks like this::

    [adapters]
    truseq1:AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    truseq2:CAAGCAGAAGACGGCATACGAGAT*GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
    
    [indexes]
    bfidt-000:AACCGAGTTA
    bfidt-001:AATACTTCCG
    bfidt-002:AATTAAGGCC
    bfidt-003:AAGCTTATCC
    
    [params]
    separate reads:True
    read1:{name}_R1.fastq.gz
    read2:{name}_R2.fastq.gz
    
    [combos]
    genus_species1:bfidt-000
    genus_species2:bfidt-001
    genus_species3:bfidt-002
    genus_species4:bfidt-003
    
    [remap]
    bfidt-000:genus_species1
    bfidt-001:genus_species2
    bfidt-002:genus_species3
    bfidt-003:genus_species4
    
This gives the adapter sequences to trim (with the index indicated by an 
asterisk), the indexes we used, whether reads are separate, a name-formatting
convention, a mapping of species names to index, and a mapping of index names
to species.

You can run illumiprocessor against your data (in `demultiplexed`) with:

.. code-block:: bash
    
    mkdir uce-clean
    python ~/git/illumiprocessor/illumiprocessor.py \
        demultiplexed \
        uce-clean \
        malfaro-fish-illumiprocesser.conf \
        --remap \
        --clean \
        --cores 12 \
        --complex
        
The clean data will appear in `uce-clean` with the following structure::

    uce-clean/
        genus_species/
            interleaved-adapter-quality-trimmed/
                genus_species-READ1and2-interleaved.fastq.gz
                genus_species-READ-singleton.fastq.gz
            stats/
                genus_species-READ1.fastq.gz-adapter-contam.txt
                genus_species--READ2.fastq.gz-adapter-contam.txt
                sickle-trim.txt
            untrimmed/
                genus_species-READ1.fastq.gz (symlink)
                genus_species-READ1.fastq.gz (symlink)
        
`interleaved` contains the cleaned read data in interleaved format, with one
file containing "READ1and2" (both reads kept) and another file containing "READ"
data, where only one read of the pair are kept (AKA singleton reads).