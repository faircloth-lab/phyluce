.. include:: global.rst

.. _Tutorial IV:

********************************************************************
Tutorial IV: Identifying UCE Loci and Designing Baits To Target Them
********************************************************************

The first few tutorials have given you a feel for how to perform
phylogenetic/phylogeograohpc analyses using existing probe sets, new data, and
genomes.  However, what if you are working on a set of organisms without a probe
set targeting conserved loci?  How do you identify those loci and design baits
to target them?  This Tutorial shows you how to do that.

In the examples below, we'll follow a UCE identification and probe design wo
rkflow using data from Coleoptera (beetles).  Although you can follow the entire
tutorial from beginning to end, I've also made the BAM files containing mapped
reads available, which lets you skip the computationally exepensive step of
performing read simulation and alignment.

Starting directory structure
============================

To keep things clear, we're going to assume you are working in some directory,
which I'll call `uce-coleoptera`.  We'll be working from the top of this
directory in the steps below:

.. code-block:: bash

    uce-coleoptera

.. _data-download:

Data download and preparation
=============================

Download the genomes
-------------------------------

.. attention:: You do not neccessarily need to do this as part of this tutorial
    for UCE identification and probe design - If you only want to follow the
    steps for locus identification (skipping probe design and in-silico
    testing), you can simply download the prepared `FASTQ/BAM files from
    figshare <https://dx.doi.org/10.6084/m9.figshare.3487349>`_.

Make a directory to hold the genome sequences:

.. code-block:: bash

    > mkdir genomes

Now, get the genome sequences:

.. code-block:: bash

    # Anoplophora glabripennis (Asian longhorned beetle)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000390285.1_Agla_1.0/GCA_000390285.1_Agla_1.0_genomic.fna.gz

    # Agrilus planipennis (emerald ash borer)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000699045.1_Apla_1.0/GCA_000699045.1_Apla_1.0_genomic.fna.gz

    # Leptinotarsa decemlineata (Colorado potato beetle)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000500325.1_Ldec_1.5/GCA_000500325.1_Ldec_1.5_genomic.fna.gz

    # Onthophagus taurus (beetles)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000648695.1_Otau_1.0/GCA_000648695.1_Otau_1.0_genomic.fna.gz

    # Dendroctonus ponderosae (mountain pine beetle)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000355655.1_DendPond_male_1.0/GCA_000355655.1_DendPond_male_1.0_genomic.fna.gz

    # Tribolium castaneum (red flour beetle)
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000002335.2_Tcas_3.0/GCA_000002335.2_Tcas_3.0_genomic.fna.gz

    # Mengenilla moldrzyki (twisted-wing parasites) [Outgroup]
    > wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000281935.1_Memo_1.0/GCA_000281935.1_Memo_1.0_genomic.fna.gz

You need to unzip all of the genome sequences

.. code-block:: bash

    > gunzip *

Finally, your directory structure should look something like:

.. code-block:: bash

    uce-coleoptera
    └── genomes
        ├── GCA_000002335.2_Tcas_3.0_genomic.fna
        ├── GCA_000281935.1_Memo_1.0_genomic.fna
        ├── GCA_000355655.1_DendPond_male_1.0_genomic.fna
        ├── GCA_000390285.1_Agla_1.0_genomic.fna
        ├── GCA_000500325.1_Ldec_1.5_genomic.fna
        ├── GCA_000648695.1_Otau_1.0_genomic.fna
        └── GCA_000699045.1_Apla_1.0_genomic.fna


Cleanup the genome sequences
----------------------------

.. attention:: You do not neccessarily need to do this as part of this tutorial
    for UCE identification and probe design - If you only want to follow the
    steps for locus identification (skipping probe design and in-silico
    testing), you can simply download the prepared `FASTQ/BAM files from
    figshare <https://dx.doi.org/10.6084/m9.figshare.3487349>`_.

When you get genome sequences from NCBI, the FASTA headers of most
scaffold/contigs contain a lot of extra cruft that can cause problems with some
of the steps in the UCE identification and probe design workflow. I usually
remove the extra stuff, maintaining only the accession number information for
each contig / scaffold.  To do that, I use a little python script that looks
like the following:

.. code-block:: python

    from Bio import SeqIO
    with open("Name_of_Genome_File.fna", "rU") as infile:
    with open("outfileName.fasta", "w") as outf:
        for seq in SeqIO.parse(infile, 'fasta'):
            seq.name = ""
            seq.description = ""
            outf.write(seq.format('fasta'))

For these genome assemblies, the important bits are in the FASTA header just
before the space in the name.  The code above basically keeps this information
before the space and discards the remaining FASTA header information.

In the case of the genomes above, here are the commands I ran (note also that
this creates an output file with a different name from the input file):

.. code-block:: python

    from Bio import SeqIO

    # Anoplophora glabripennis (Asian longhorned beetle)
    with open("GCA_000390285.1_Agla_1.0_genomic.fna", "rU") as infile:
        with open("anoGla1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Agrilus planipennis (emerald ash borer)
    with open("GCA_000699045.1_Apla_1.0_genomic.fna", "rU") as infile:
        with open("agrPla1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Dendroctonus ponderosae (mountain pine beetle)
    with open("GCA_000355655.1_DendPond_male_1.0_genomic.fna", "rU") as infile:
        with open("denPon1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Leptinotarsa decemlineata (Colorado potato beetle)
    with open("GCA_000500325.1_Ldec_1.5_genomic.fna", "rU") as infile:
        with open("lepDec1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Mengenilla moldrzyki (twisted-wing parasites) [Outgroup]
    with open("GCA_000281935.1_Memo_1.0_genomic.fna", "rU") as infile:
        with open("menMol1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Onthophagus taurus (beetles)
    with open("GCA_000648695.1_Otau_1.0_genomic.fna", "rU") as infile:
        with open("ontTau1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))

    # Tribolium castaneum (red flour beetle)
    with open("GCA_000002335.2_Tcas_3.0_genomic.fna", "rU") as infile:
        with open("triCas1.fasta", "w") as outf:
            for seq in SeqIO.parse(infile, 'fasta'):
                seq.name = ""
                seq.description = ""
                outf.write(seq.format('fasta'))


Now, you can remove the original files downloaded from NCBI:

.. code-block:: bash

    rm *.fna

And, your directory structure should look something like:

.. code-block:: bash

    uce-coleoptera
    └── genomes
        ├── anoGla1.fasta
        ├── agrPla1.fasta
        ├── denPon1.fasta
        ├── lepDec1.fasta
        ├── menMol1.fasta
        ├── ontTau1.fasta
        └── triCas1.fasta


Put genomes in their own directories
------------------------------------

Because of some historical reasons (and how I organize our lab data), each
genome sequence needs to be in its own directory.  We can do that pretty easily
by running:

.. code-block:: bash

    > cd uce-coleoptera
    > cd genomes
    > for critter in *; do mkdir ${critter%.*}; mv $critter ${critter%.*}; done

Now the directory structure looks like:

.. code-block:: text

    uce-coleoptera
    └── genomes
        ├── agrPla1
        │   └── agrPla1.fasta
        ├── anoGla1
        │   └── anoGla1.fasta
        ├── denPon1
        │   └── denPon1.fasta
        ├── lepDec1
        │   └── lepDec1.fasta
        ├── menMol1
        │   └── menMol1.fasta
        ├── ontTau1
        │   └── ontTau1.fasta
        └── triCas1
            └── triCas1.fasta


Convert genomes to 2bit format
------------------------------

Later during the workflow, we're going to need to have our genomes in `2bit`
format, which is a binary format for genomic data that is easier and faster to
work with than FASTA format.  We'll use a BASH script to convert all of the
sequences, above, to `2bit` format:

.. code-block:: bash

    > cd uce-coleoptera
    > cd genomes
    > for critter in *; do faToTwoBit $critter/$critter.fasta $critter/${critter%.*}.2bit; done

Now the directory structure looks like:

.. code-block:: text

    uce-coleoptera
    └── genomes
        ├── agrPla1
        │   ├── agrPla1.2bit
        │   └── agrPla1.fasta
        ├── anoGla1
        │   ├── anoGla1.2bit
        │   └── anoGla1.fasta
        ├── denPon1
        │   ├── denPon1.2bit
        │   └── denPon1.fasta
        ├── lepDec1
        │   ├── lepDec1.2bit
        │   └── lepDec1.fasta
        ├── menMol1
        │   ├── menMol1.2bit
        │   └── menMol1.fasta
        ├── ontTau1
        │   ├── ontTau1.2bit
        │   └── ontTau1.fasta
        └── triCas1
            ├── triCas1.2bit
            └── triCas1.fasta

Simulate reads from genomes
---------------------------

.. attention:: You do not neccessarily need to take this step as part of the
    tutorial - you can simply download prepared, simulated `FASTQ files from
    figshare <https://ndownloader.figshare.com/files/5513204>`_.

In order to locate UCE loci across a selection of different genomes, we're going
to align reads from each taxon, above, to a reference genome sequence (the
"base" genome sequence) using a permissive raw read aligner.  You can use reads
from low-coverage, genome scans or you can use reads simulated from particular
genomes.  In this tutorial, we're going to use this latter approach and simulate
reads (without sequencing error) from the genomes that we will align to the base
genome. To accomplish this, we'll use `art
<http://www.niehs.nih.gov/research/resources/software/biostatistics/art/>`_,
which is a robust read simulator that is reasonably flexible.

Because we're using simulated reads to locate UCE loci, we want to turn off the
built-in feature of `art` that adds some sequencing error to simulated reads.
This results in a general form of the call to `art` that looks like:

.. code-block:: bash

    > art_illumina \
        --paired \
        --in ~/path/to/input/genome.fasta \
        --out prefix-of-output-file \
        --len 100 \
        --fcov 2 \
        --mflen 200 \
        --sdev 150 \
        -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

This simulates reads from the `--in genome.fasta` that are 100 bp in length,
cover the genome randomly to roughly 2X, have an insert size of 200 bp, and have
an inserts size standard deviation of 150.  The last line turns off the
simulation of sequencing error in each of these reads and also turns off the
creation of an alignment file showing where the reads came from in the genome
sequence.

In the case of the Coleoptera genomes you downloaded, here are the commands I
ran to simulate the read data that we will use in the next step (note also that
this creates an output file with a different name from the input file).  First,
create a directory to hold the simulated read data:

.. code-block:: bash

    > cd uce-coleoptera
    > mkdir reads
    > cd reads

Then, simulate the reads using `art`:

.. code-block:: bash

    > art_illumina \
        --paired \
        --in ../genomes/agrPla1/agrPla1.fasta \
        --out agrPla1-pe100-reads \
        --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

    > art_illumina \
        --paired \
        --in ../genomes/anoGla1/anoGla1.fasta \
        --out anoGla1-pe100-reads \
        --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

    > art_illumina \
        --paired \
        --in ../genomes/denPon1/denPon1.fasta \
        --out denPon1-pe100-reads \
        --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

    > art_illumina \
        --paired \
        --in ../genomes/lepDec1/lepDec1.fasta \
        --out lepDec1-pe100-reads \
        --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

    > art_illumina \
        --paired \
        --in ../genomes/ontTau1/ontTau1.fasta \
        --out ontTau1-pe100-reads \
        --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na


Now, you should see 2 read files for each taxon.  We are going to "break" the
pairs by merging the read information together, then we are going to zip the
resulting file that contains the `R1` and `R2` data.  We can accomplish this
pretty easily using a quick BASH script:

.. attention:: Note that we've dropped menMol1 from the read simulation
    process.  This is largely because it is an outgroup to beetles.  We'll use
    it later, when we're performing *in silico* tests of the UCE bait set.

.. code-block:: bash

    for critter in agrPla1 anoGla1 denPon1 lepDec1 ontTau1;
        do
            echo "working on $critter";
            touch $critter-pe100-reads.fq;
            cat $critter-pe100-reads1.fq > $critter-pe100-reads.fq;
            cat $critter-pe100-reads2.fq >> $critter-pe100-reads.fq;
            rm $critter-pe100-reads1.fq;
            rm $critter-pe100-reads1.fq;
            gzip $critter-pe100-reads.fq;
        done;

If we take a look at our directory structure, it should look like:

.. code-block:: text

    uce-coleoptera
    ├── genomes (collapsed)
    └── reads
        ├── anoGla1-pe100-reads.fq.gz
        ├── agrPla1-pe100-reads.fq.gz
        ├── denPon1-pe100-reads.fq.gz
        ├── lepDec1-pe100-reads.fq.gz
        └── ontTau1-pe100-reads.fq.gz

.. _uce-read-alignment:

Read alignment to the base genome
==================================

.. attention:: You do not neccessarily need to run this step as part of the
    tutorial - you can simply download the prepared, `BAM files from
    figshare <https://ndownloader.figshare.com/files/5513489>`_.

    Because we also provide the BAM files created below, you can choose to just
    start with the BAM files in the :ref:`uce-identification` section.


Prepare the base genome
-----------------------

Now that we have read data representing each of our exemplar taxa, we need to
align these reads to the "base" genome sequence, in this case the genome sequence
of *Tribolium castaneum* (aka `triCas1`).  We selected this assembly as the "base"
genome because of it's age (i.e., better-assembled) and level of annotation.

We will perform the read alignments to `triCas1` using the permissive read
aligner, stampy_, which works well when aligning sequences to a divergent
reference sequence.  Hoever, before running the alignments, we need to prepare
the base genome.  And, before we do that, let's create a direcetory to work in:

.. code-block:: bash

    > cd uce-coleoptera
    > mkdir base
    > cd base

Now, let's copy the base genome to this directory (for simplicity):

.. code-block:: bash

    > cp ../genomes/triCas1.fasta ./

If we take a look at our directory structure, it now looks like:

.. code-block:: text

    uce-coleoptera
    ├── base
    │   └── triCas1.fasta
    ├── genomes
    │   ├── agrPla1
    │   │   ├── agrPla1.2bit
    │   │   └── agrPla1.fasta
    │   ├── anoGla1
    │   │   ├── anoGla1.2bit
    │   │   └── anoGla1.fasta
    │   ├── denPon1
    │   │   ├── denPon1.2bit
    │   │   └── denPon1.fasta
    │   ├── lepDec1
    │   │   ├── lepDec1.2bit
    │   │   └── lepDec1.fasta
    │   ├── menMol1
    │   │   ├── menMol1.2bit
    │   │   └── menMol1.fasta
    │   ├── ontTau1
    │   │   ├── ontTau1.2bit
    │   │   └── ontTau1.fasta
    │   └── triCas1
    │       ├── triCas1.2bit
    │       └── triCas1.fasta
    └── reads
        ├── anoGla1-pe100-reads
        ├── agrPla1-pe100-reads
        ├── denPon1-pe100-reads
        ├── lepDec1-pe100-reads
        └── ontTau1-pe100-reads

Now, we need to run the commands to prepare the `triCas1` genome for use with
stampy:

.. code-block:: bash

    > cd uce-coleoptera
    > cd base
    > stampy.py --species="tribolium-castaneum" --assembly="triCas1" -G triCas1 triCas1.fasta
    > stampy.py -g triCas1 -H triCas1

If we look at our directory structure, it should look like:

.. code-block:: text

    uce-coleoptera
    ├── base
    │   ├── triCas1.fasta
    │   ├── triCas1.sthash
    │   └── triCas1.stidx
    ├── genomes (collapsed)
    └── reads
        ├── anoGla1-pe100-reads
        ├── agrPla1-pe100-reads
        ├── denPon1-pe100-reads
        ├── lepDec1-pe100-reads
        └── ontTau1-pe100-reads

Align reads to the base genome
------------------------------

.. attention:: You do not neccessarily need to run this step as part of the
    tutorial - you can simply download the prepared, `BAM files from
    figshare <https://ndownloader.figshare.com/files/5513489>`_.

    Because we also provide the BAM files created below, you can choose to just
    start with the BAM files in the :ref:`uce-identification` section.

Now that we've prepared our base genome, we need to perform the actual alignment
of the simulated reads to the base genome. And, before we do that, let's create
a directory to hold the resulting alignments:

.. code-block:: bash

    > cd uce-coleoptera
    > mkdir alignments

Our resulting directory structure should now look like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    ├── base
    │   ├── triCas1.fasta
    │   ├── triCas1.sthash
    │   └── triCas1.stidx
    ├── genomes (collapsed)
    └── reads
        ├── anoGla1-pe100-reads
        ├── agrPla1-pe100-reads
        ├── denPon1-pe100-reads
        ├── lepDec1-pe100-reads
        └── ontTau1-pe100-reads

Now, we need to perform the alignments on a taxon-by-taxon basis (and/or you can
run these in parallel using HPC).  To do this easily (and on a local
computer) we can use a BASH script to run the alignments serially:

.. warning:: Note that I am using 16 physical CPU cores (`$cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

.. code-block:: bash

    export cores=16
    export base=triCas1
    export base_dir=$HOME/uce-coloeptera/alignments
    for critter in agrPla1 anoGla1 denPon1 lepDec1 ontTau1;
        do
            export reads=$critter-pe100-reads.fq.gz;
            mkdir -p $base_dir/$critter;
            cd $base_dir/$critter;
            stampy.py --maxbasequal 93 -g ../../base/$base -h ../../base/$base
            --substitutionrate=0.05 -t$cores --insertsize=400 -M
            ../../reads/$reads | samtools view -Sb - > $critter-to-$base.bam;
        done;

This code basically loops over each exemplar genomes, aligns the reads to the
base genome sequence, and converts the resulting output to BAM_ format (which
is a binary, compressed version of SAM_ format).

When the alignments have completed, your directory structure should look
something like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │    ├── anoGla1
    │    │   └── anoGla1-to-triCas1.bam
    │    ├── agrPla1
    │    │   └── agrPla1-to-triCas1.bam
    │    ├── denPon1
    │    │   └── denPon1-to-triCas1.bam
    │    ├── lepDec1
    │    │   └── lepDec1-to-triCas1.bam
    │    └── ontTau1
    │        └── ontTau1-to-triCas1.bam
    ├── base
    │   ├── triCas1.fasta
    │   ├── triCas1.sthash
    │   └── triCas1.stidx
    ├── genomes (collapsed)
    └── reads
        ├── anoGla1-pe100-reads
        ├── agrPla1-pe100-reads
        ├── denPon1-pe100-reads
        ├── lepDec1-pe100-reads
        └── ontTau1-pe100-reads

Now, these BAM_ files are pretty large because they contain all mapped as well
as all unmapped reads. We want to remove those unmapped reads, which will also
reduce file-size.  We can do that using samtools_ `view` and a BASH script.
We're also going to create a directory named `all` and symlink all of the
reduced BAM files to this directory.

.. warning:: The script, as-written, removes the BAM files containing both
    mapped and unmapped reads.  If you don't want to do this, remove the `rm
    $critter/$critter-to-triCas1.bam;` line.

.. code-block:: bash

    > cd uce-coleoptera
    > cd alignments
    > mkdir all
    > for critter in agrPla1 anoGla1 denPon1 lepDec1 ontTau1;
        do
            samtools view -h -F 4 -b $critter/$critter-to-triCas1.bam > $critter/$critter-to-triCas1-MAPPING.bam;
            rm $critter/$critter-to-triCas1.bam;
            ln -s ../$critter/$critter-to-triCas1-MAPPING.bam all/$critter-to-triCas1-MAPPING.bam;
        done;

Now, your directory structure should look something like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │    ├── anoGla1
    │    │   └── anoGla1-to-triCas1-MAPPING.bam
    │    ├── all
    │    │   ├── agrPla1-to-triCas1-MAPPING.bam -> ../agrPla1/agrPla1-to-triCas1-MAPPING.bam
    │    │   ├── anoGla1-to-triCas1-MAPPING.bam -> ../anoGla1/anoGla1-to-triCas1-MAPPING.bam
    │    │   ├── denPon1-to-triCas1-MAPPING.bam -> ../denPon1/denPon1-to-triCas1-MAPPING.bam
    │    │   ├── lepDec1-to-triCas1-MAPPING.bam -> ../lepDec1/lepDec1-to-triCas1-MAPPING.bam
    │    │   └── ontTau1-to-triCas1-MAPPING.bam -> ../ontTau1/ontTau1-to-triCas1-MAPPING.bam
    │    ├── agrPla1
    │    │   └── agrPla1-to-triCas1-MAPPING.bam
    │    ├── denPon1
    │    │   └── denPon1-to-triCas1-MAPPING.bam
    │    ├── lepDec1
    │    │   └── lepDec1-to-triCas1-MAPPING.bam
    │    └── ontTau1
    │        └── ontTau1-to-triCas1-MAPPING.bam
    ├── base
    │   ├── triCas1.fasta
    │   ├── triCas1.sthash
    │   └── triCas1.stidx
    ├── genomes (collapsed)
    └── reads
        ├── anoGla1-pe100-reads
        ├── agrPla1-pe100-reads
        ├── denPon1-pe100-reads
        ├── lepDec1-pe100-reads
        └── ontTau1-pe100-reads

These `*-MAPPING.bam files are available from figshare
<https://ndownloader.figshare.com/files/5513489>`_.

What it all means
-----------------

Basically, because we've now mapped simulated (or actual) sequence data from
several exemplar taxa to a closely-related "base" genome sequence, we've
essentially identified putatively orthologous loci shared between the exemplar
taxa and the base taxon. These conserved regions are where the simulated (or
actual) sequene data mapped with a sequence divergence of < 5%.

Now, we need to filter this large number of conserved regions to remove things
like repetitive regions, but also to find *which* loci are shared among *all*
exemplar taxa - not simply a single exemplar and the base taxon.

.. _uce-identification:

Conserved locus identification
==============================

You can download the `alignment data generated in the steps above from figshare
<https://ndownloader.figshare.com/files/5513489>`_ for use in subsequent steps,
rather than generating them youself (thus, saving you time).

.. attention:: If you are starting the tutorial at this position after
    downloading the `*-MAPPING.bam files from figshare
    <https://ndownloader.figshare.com/files/5513489>`_, then you will need to
    create a directory to work in named `uce-coleoptera` and then place all of
    the `*-MAPPING.bam` files in a subdirectory of `uce-coleoptera` names
    `alignments/all`.  Your resulting directory structure should look like:

    .. code-block:: text

        uce-coleoptera
        ├── alignments
        │   └── all
        │       ├── agrPla1-to-triCas1-MAPPING.bam -> ../agrPla1/agrPla1-to-triCas1-MAPPING.bam
        │       ├── anoGla1-to-triCas1-MAPPING.bam -> ../anoGla1/anoGla1-to-triCas1-MAPPING.bam
        │       ├── denPon1-to-triCas1-MAPPING.bam -> ../denPon1/denPon1-to-triCas1-MAPPING.bam
        │       ├── lepDec1-to-triCas1-MAPPING.bam -> ../lepDec1/lepDec1-to-triCas1-MAPPING.bam
        │       └── ontTau1-to-triCas1-MAPPING.bam -> ../ontTau1/ontTau1-to-triCas1-MAPPING.bam
        └── genomes (collapsed)

    If you want to go beyond conserved locus identification and design probes
    from the target taxa, you will also need to download the appropriate
    genomes.  See the :ref:`data-download` section.

Convert BAMS to BEDS
--------------------

In the steps above, we have generated BAM files representing reads that stampy_
has mapped to the base genome.  Those reads that map align to putatively
conserved sequence regions (that we need to filter), and these alignments of
mapping reads should reside in our `alignments/all` directory.

To begin the filtering process, we're going to convert each BAM_  file to BED_
format, which is an interval-based format that is easy and fast to manipulate
with a software suite known as bedtools_.  But before we do that, we are doing
to create a `bed` directory to hold all of these BED_ format files.

.. code-block:: bash

    > cd uce-coleoptera

    # make a directory to hold the BED data
    > mkdir bed
    > cd bed

    # now, convert our *-MAPPING.bam files to BED format
    > for i in ../alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done

Your directory structure should look something like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all
    │       ├── agrPla1-to-triCas1-MAPPING.bam
    │       ├── anoGla1-to-triCas1-MAPPING.bam
    │       ├── denPon1-to-triCas1-MAPPING.bam
    │       ├── lepDec1-to-triCas1-MAPPING.bam
    │       └── ontTau1-to-triCas1-MAPPING.bam
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   └── ontTau1-to-triCas1-MAPPING.bam.bed
    └── genomes (collapsed)

Sort the converted BEDs
-----------------------

Before moving forward with the `merge` command, below, we need to sort the
resulting BED files, which orders each lines of data in the BED_ file by
chromosome/scaffold/contig and position along that chromosome/scaffold/contig.
Again, we can do this with some bash scripting:

.. code-block:: bash

    > cd uce-coleoptera
    > cd bed
    > for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done

Your directory structure should look something like the following (note that I
have collapsed the directory listing for `all`):

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   └── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    └── genomes (collapsed)


Merge overlapping or nearly-overlapping intervals
-------------------------------------------------

Because there may be small gaps between proximate regions of conservation
(which may result because we're using data that are either from low-coverage,
simulated reads or low-coverage *actual* reads) we need to merge together
putative regions of conservation that are proximate.  Luckily bedtools_ has
a handy tool to do that - the `merge` function.

.. code-block:: bash

    > cd uce-coleoptera
    > cd bed
    > for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done

Your directory structure should look something like the following (note that I
have collapsed the directory listing for `all`):

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    │   └── ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed
    └── genomes (collapsed)

To get some idea of the total number of merged, putatively conserved regions in
each exemplar taxon that are shared with the base genome, we can simply loop
over the files and count the number of lines in each:

.. code-block:: bash

    > cd uce-coleoptera
    > cd bed

    > for i in *.bam.sort.merge.bed; do wc -l $i; done
    19810 agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    48350 anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    21390 denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    33144 lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    25188 ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed


Remove repetitive intervals
---------------------------

At this point, we've mapped reads to the base genome, kept those regions where
reads mapped, converted that to a BED_, and merged intervals that are very close
to one another.

What we have not done is remove any putatively conserved intervals shared
between exemplar taxa and the base genome that are repetitive regions.  To do
this, we're going to run a python program over all of the BED files for each
exemplar taxon.  This program will look at the intervals shared between the
exemplar taxon and the base genome and it will remove intervals where the base
genome is shorter than 80 bp and where > 25 % of the base genome is masked.

.. code-block:: bash

    > cd uce-coleoptera
    > cd bed

    > for i in *.sort.merge.bed;
        do
            phyluce_probe_strip_masked_loci_from_set \
                --bed $i \
                --twobit ../genomes/triCas1/triCas1.2bit \
                --output ${i%.*}.strip.bed \
                --filter-mask 0.25 \
                --min-length 80
        done;

    Screened 19810 sequences from agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 3113 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 16697.
    Screened 48350 sequences from anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 13226 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 35124.
    Screened 21390 sequences from denPon1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 3008 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 18382.
    Screened 33144 sequences from lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 6585 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 26559.
    Screened 9637 sequences from menMol1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 4379 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 5258.
    Screened 25188 sequences from ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed.  Filtered 6505 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 18683.

When this finishes, your directory structure should look like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   └── ontTau1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    └── genomes (collapsed)

Determining locus presence in multiple genomes
----------------------------------------------

Up to this point, we've been processing each file on a taxon-by-taxon basis,
where each taxon had data aligned to the base genome.  Now, we need to determine
which loci are conserved *across* taxa.  To do that, we first need to prepare a
configuration file (named `bed-files.conf`) that gives the paths to each of our
`*.bam.sort.merge.strip.bed` files.  That file needs to be in configuration file
format, like so:

.. code-block:: text

    [beds]
    agrPla1:agrPla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    anoGla1:anoGla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    denPon1:denPon1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    lepDec1:lepDec1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    ontTau1:ontTau1-to-triCas1-MAPPING.bam.sort.merge.strip.bed

The `[beds]` line is the "header" line, and that is followed by each taxon name
(on the left) and the name of the BED file we want to process (on the right).
You should place this file in the `bed` directory.  If you place it elsewhere,
you'll need to use full paths on the right hand side.

Your directory structure should now look like (note new `bed-files.conf`)

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── bed-files.conf
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   └── ontTau1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    └── genomes (collapsed)

Now, we're going to run the following program, that creates a record of which
alignment intervals are shared among taxa.  We need to pass the location of the
`bed-files.conf` to this program, along with the name of our base genome, and a
name for the output database that will be created:

.. code-block:: text

    > phyluce_probe_get_multi_merge_table \
        --conf bed-files.conf \
        --base-taxon triCas1 \
        --output coleoptera-to-triCas1.sqlite

    agrpla1.................
    anogla1....................................
    denpon1...................
    lepdec1...........................
    onttau1...................
    Creating database
    Inserting results

The program shows the results of inserting data for each exemplar taxon that
we've selected.  If we take a look at the table contents (see :ref:`Database`
for more instructions on sqlite_ databases), we see something like the
following:

.. code-block:: text

    sqlite> select * from triCas1 limit 10;
    uce         chromo      start       stop        agrpla1     anogla1     denpon1     lepdec1     onttau1
    ----------  ----------  ----------  ----------  ----------  ----------  ----------  ----------  ----------
    1           GG695505.1  532         632         1           1           0           0           1
    2           GG695547.1  826         926         0           1           0           0           0
    3           GG695547.1  1121        1221        0           1           0           0           0
    4           GG695547.1  1293        1393        0           1           0           0           0
    5           GG694821.1  1002        1102        0           0           0           1           0
    6           GG695519.1  73          193         0           1           1           0           0
    7           GG695519.1  222         380         1           0           1           0           1
    8           GG695519.1  925         1129        1           1           1           0           1
    9           DS497688.1  17907       18022       0           0           1           0           0
    10          DS497688.1  19840       19934       0           1           0           0           0

The first row of this table (which is limited to 10 rows of results by the query
although it is 60699 rows long) shows that for `triCas1` contig `GG695505.1`,
`agrpla1`, `anogla1`, and `onttau1` have reads that overlap at position 532 to
632.

Determining shared, conserved, loci
-----------------------------------

Now that we have our table of results, we can run a quick query (using a
Python_ program) against the table to look at results, more generally.  The
following code queries the database and writes out the number of loci shared by
the base taxon (`triCas1`) and 1, 2, 3, 4, and 5 (all) of the exemplar taxa that
we've aligned to the base genome.  You need to give the program the path to the
database created above and the name of the base taxon:

.. code-block:: python

    > python phyluce_probe_query_multi_merge_table \
            --db coleoptera-to-triCas1.sqlite  \
            --base-taxon triCas1

    Loci shared by triCas1 + 0 taxa:    60,699.0
    Loci shared by triCas1 + 1 taxa:    60,699.0
    Loci shared by triCas1 + 2 taxa:    32,431.0
    Loci shared by triCas1 + 3 taxa:    15,834.0
    Loci shared by triCas1 + 4 taxa:    6,471.0
    Loci shared by triCas1 + 5 taxa:    1,822.0

The output from this program basically says that, if we are interested in only
those loci found in all exemplar taxa that align to the base genome, there are
1,822 of those.  Similarly, if we're willing to be a little less strict about
things, there are 6,471 conserved loci that are shared by triCas1 and 4 of the
exemplar taxa.

.. admonition:: Question: How conservative should I be?
    :class: admonition tip

    Basically, the question boils down to "Should I select only the set of loci
    shared by all exemplars and the base genome or shoul I be more liberal?".
    It's also a hard question to answer.  In most cases, I'm pretty happy
    selecting *n-1* or *n-2* where *n* is the total number of exemplar taxa.  In
    the example below, however, we've selected *n* as the "ideal".  This is
    largely because we have so little information about coleopteran genomes - so
    we want to be pretty darn sure these loci are found in most/all of them.

Now that we have a general sense of the number of conserved loci in each class
of sharing across exemplars (e.g. 5 (all), 4, 3, 2, 1), we need to extract those
loci that fall within one of these classes.  In this case (and as noted in the
box, above), we're going to ouput only those conserved loci that we've
identified as being shared between the base genome and *all* exemplars.  We do
that with a slightly different Python_ script.  This script takes the database
name, the base genome, the count of exemplar taxa shared across, and the name of
an output file as input.  The output file will be BED formatted.

.. code-block:: python

    > phyluce_probe_query_multi_merge_table \
            --db coleoptera-to-triCas1.sqlite \
            --base-taxon triCas1 \
            --output triCas1+5.bed \
            --specific-counts 5

    Counter({'anogla1': 1822, 'lepdec1': 1822, 'agrpla1': 1822, 'denpon1': 1822, 'onttau1': 1822})


Conserved locus validation
==========================

Extract FASTA sequence from base genome for temp bait design
------------------------------------------------------------

Now that we've indentified conserved sequences shared among the base genome and
the exemplar taxa, we need to start designing baits to capture these loci.  The
first step in this process is to extract FASTA sequences from the base genome
that correspond to the loci we've identified.  We do that with a Python_ script
that takes as input the BED_ file we created, above, the `2bit`-formatted base
genomes, a length of sequence we want to extract (160 bp), and an output FASTA
filename.

.. admonition:: Question: Why buffer to 160 bp?
    :class: admonition tip

    We are extracting FASTA regions of 160 bp because that allows us to place 2
    baits right in the center of this region at 3x tiling density which means
    that standard 120 bp baits will overlap by 40 bp and have 80 bp to each side
    (total length 160 bp).

To run the code, we use:

.. code-block:: python

    > phyluce_probe_get_genome_sequences_from_bed \
            --bed triCas1+5.bed  \
            --twobit ../genomes/triCas1/triCas1.2bit \
            --buffer-to 160 \
            --output triCas1+5.fasta

    Screened 1822 sequences.  Filtered 7 < 160 bp or with > 25.0% masked bases or > 0 N-bases. Kept 1815.

That should produce a fasta file whose contents look similar to:

.. code-block:: text

    >slice_0 |DS497688.1:249724-250111
    AAAATCAAAGTCGAATACAAAGGCGAATCTAAGACTTTCTATCCTGAAGAGATCAGTTCC
    ATGGTacttacaaaaatgaaggaaacTGCCGAAGCCTATTTAGGCAAATCGGTCACAAAT
    GCCGTTATCACCGTACCAGCCTATTTCAACGATTCGCAAAGGCAGGCAACTAAAGATGCC
    GGTACTATTTCCGGCTTGCAAGTTTTGCGTATTATTAACGAACCTACGGCTGCTGCCATT
    GCCTACGGTTTGGATAAGAAGGGAACTGGGGAACGTAATGTCTTGATTTTTGATCTGGGT
    GGTGGTACTTTTGATGTGAGCATTTTGACCATTGAGGATGGCATTTTCGAGGTCAAGTCC
    ACCGCTGGTGATACGCATTTGGGTGGC
    >slice_1 |DS497688.1:250513-250673
    CCTGATGAGGCTGTTGCCTATGGAGCTGCCGTCCAAGCCGCCATTTTGCACGGTGATAAG
    TCGGAAGAGGTTCAAGATTTGCTACTTTTGGACGTTACTCCACTTTCATTGGGTATTGAA
    ACAGCAGGCGGTGTGATGACTGCTTTGATCAAGCGTAACA
    >slice_2 |DS497688.1:250682-250991
    CAACCAAACAAACGCAAACTTTCACCACCTACTCTGATAACCAACCCGGTGTATTGATCC
    AAGTGTACGAAGGCGAACGAGCGATGACTAAAGACAATAACCTTTTGGGTAAATTCGAAT
    TGACTGGAATCCCACCGGCACCAAGAGGTGTTCCCCAAATCGAAGTCACCTTTGATATTG
    ACGCCAACGGGATTTTGAACGTCACAGCCATCGAGAAGAGCACCAACAAGGAGAACAAAA
    TCACCATCACCAATGATAAGGGACGTTTGAGCAAGGAAGATATCGAACGGATGGTCAACG
    AAGCCGAGA


Design a temporary bait set from the base taxon
-----------------------------------------------

Now that we've extracted the appropriate loci from the base genome, we need to
design bait sequences targeting these loci.  For that, we use a different
Python_ script.  This program takes as input the FASTA file we just created, and
some design-specific information (`--probe-prefix, --design, --designer`).  The
design options (`--tiling-density, --two-probes, --overlap`) ensure that we
select two baits per locus with 3x tiling that overlap the middle of the
targeted locus.  Finally, we remove (`--masking, --remove-gc`) potentially
problematic baits with >25% repeat content and GC content outside of the range
of 30-70% (30 % > GC > 70%).

.. code-block:: python

    > phyluce_probe_get_tiled_probes \
        --input triCas1+5.fasta \
        --probe-prefix "uce-" \
        --design coleoptera-v1 \
        --designer faircloth \
        --tiling-density 3 \
        --two-probes \
        --overlap middle \
        --masking 0.25 \
        --remove-gc \
        --output triCas1+5.temp.probes

    Probes removed for masking (.) / low GC % (G) / ambiguous bases (N):
    GGGGGGGGGGGGGGGGGGGGGGGGGGGG


    Conserved locus count = 1805
    Probe Count = 3602

Remove duplicates from our temporary bait set
---------------------------------------------

Because we haven't search for duplicates among our loci and because reducing
longer reads to shorter ones (e.g. designing baits from loci) can introduce
duplicate baits, we need to screen the resulting bait set for duplicates.  To do
that, we follow a 2-stage process - first to align all probes to themselves then
to use those alignments to remove potentially duplicates baits/loci.  First we
run a lastz_ search of all baits to themselves.  This program takes as input the
temp probes we just designed (as both `--target` and `--query`), relatively low
values for `--identity` and `--coverage` to make sure we identify as many
duplicates as possible, and the program writes these results to the `--output`
file:

.. code-block:: python

    > phyluce_probe_easy_lastz \
        --target triCas1+5.temp.probes \
        --query triCas1+5.temp.probes \
        --identity 50 --coverage 50 \
        --output triCas1+5.temp.probes-TO-SELF-PROBES.lastz

    Started:  Fri Jun 03, 2016  13:57:54
    Ended:  Fri Jun 03, 2016  13:57:55
    Time for execution:  0.0284410158793 minutes

Now that we've run the alignments, we need to screen them alignments and remove
the duplicate baits from the bait set.  This program takes as input the lastz_
results from above and the temp-probe file, as well as the probe-prefix that we
used during probe design, above.  The results are written to a file that is
equivalent to the probe file name + `DUPE-SCREENED`, so in this case the output
file is named `triCas1+5.temp-DUPE-SCREENED.probes`.

.. code-block:: python

    > phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
        --fasta triCas1+5.temp.probes  \
        --lastz triCas1+5.temp.probes-TO-SELF-PROBES.lastz \
        --probe-prefix=uce-

    Parsing lastz file...
    Screening results...
    Screened 3601 fasta sequences.  Filtered 292 duplicates. Kept 3019.


Align baits against exemplar genomes
------------------------------------

Now that we have a duplicate-free (or putatively duplicate free) set of
temporary baits designed from conserved loci in the base genome, we're going to
use some in-silico alignments to see if we can locate these loci in the several
exemplar genomes.

.. attention:: For the following analyses, you need genome assemblies for each
    of the exemplar taxa, formatted as `2bit` files.

We'll use the results of these alignments to design a bait set that includes
baits designed from the base genome, but also from the exemplar taxa.  This
should allow our bait set to work more consistently across broad groups of
organisms.

In terms of directory structure, things should look pretty similar to the
following:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── bed-files.conf
    │   ├── coleoptera-to-triCas1.sqlite
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── triCas1+5.bed
    │   ├── triCas1+5.bed.missing.matrix
    │   └── triCas1+5.fasta
    └── genomes
        ├── agrPla1
        │   ├── agrPla1.2bit
        │   └── agrPla1.fasta
        ├── anoGla1
        │   ├── anoGla1.2bit
        │   └── anoGla1.fasta
        ├── denPon1
        │   ├── denPon1.2bit
        │   └── denPon1.fasta
        ├── lepDec1
        │   ├── lepDec1.2bit
        │   └── lepDec1.fasta
        ├── menMol1
        │   ├── menMol1.2bit
        │   └── menMol1.fasta
        ├── ontTau1
        │   ├── ontTau1.2bit
        │   └── ontTau1.fasta
        └── triCas1
            ├── triCas1.2bit
            └── triCas1.fasta

Note that we have all the genomes in their directory, in both FASTA and `2bit`
formats. We're also have a new genome sequence in here - that of `menMol1`
(Mengenilla moldrzyki [twisted-wing parasites]), which represents the outgroup
to Coleoptera.  We're adding this taxon because it helps us bridge the base of
the tree - e.g. the divergence between the outgroup and the exemplar taxa that
we're using to design probes.

.. admonition:: Question: What exemplar taxa should I use for bait design?
    :class: admonition tip

    This is a really hard question to answer.  In old, divergent groups with few
    genomic resources, the answer is usually "all the species" with genomic
    data.  Basically, you want to include exemplars that make the divergence
    among baits targeting the same loci something >20% or so.  That said, even
    this number is a bit of a guess - no one has systematically tested how
    "sticky" baits are when they are used to enrich loci across divergent
    groups.  We know they are pretty sticky and in certain cases can enrich loci
    as much as 35%-40% divergent from the bait sequence.  Generally speaking, I
    try to include exemplar taxa during probe design that bridge the known
    diversity of a given group... again, in many cases this is hard (or
    impossible) to know given current data.  So, you may have to take a bit of a
    guess.

So, assuming that you have the appropriate `2bit` files in
`uce-coleoptera/genomes`, we are going to align the temporary probes that we've
designed to the exemplar genomes, and we're going to run these and subsequent
bait design steps in a new directory, named `probe-design`.  So:

.. code-block:: bash

    > cd uce-coleoptera
    > mkdir probe-design

Now, you're directory structure should look like:

.. code-block:: text

    uce-coleoptera
    ├── alignments
    │   └── all (collapsed)
    ├── bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── agrPla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── anoGla1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── bed-files.conf
    │   ├── coleoptera-to-triCas1.sqlite
    │   ├── denPon1-to-triCas1-MAPPING.bam.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── denPon1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── lepDec1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── menMol1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.bed
    │   ├── ontTau1-to-triCas1-MAPPING.bam.sort.merge.strip.bed
    │   ├── triCas1+5.bed
    │   ├── triCas1+5.bed.missing.matrix
    │   └── triCas1+5.fasta
    ├── genomes (collapsed)
    └── probe-design

We need to align the temporary probe sequences to each genome, which we can do
using the following code, which takes as input our temporary probe file, the
list of genomes we want to align the probes against, the path to the genomes,
the minimum sequence identity to accept a match (on the low end of the spectrum
for this step), and number of compute cores to use, and the name of an output
database to create and the output directory in which to store the lastz results.

.. warning:: Note that I am using 16 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

.. code-block:: python

    > phyluce_probe_run_multiple_lastzs_sqlite \
        --probefile ../bed/triCas1+5.temp-DUPE-SCREENED.probes \
        --scaffoldlist agrPla1 anoGla1 denPon1 lepDec1 ontTau1 triCas1 menMol1 \
        --genome-base-path ../genomes \
        --identity 50 \
        --cores 16 \
        --db triCas1+5+menMol1.sqlite \
        --output coleoptera-genome-lastz

    Running against agrPla1.2bit
    Running with the --huge option.  Chunking files into 10000000 bp

    < ... snip ... >

    Cleaning up the chunked files...
    Cleaning /nfs/data1/working/bfaircloth-insects/coleoptera/temp/probe-design/coleoptera-genome-lastz/triCas1+5.temp-DUPE-SCREENED.probes_v_menMol1.lastz
    Creating menMol1 table
    Inserting data to menMol1 table

Extract sequence around conserved loci from exemplar genomes
------------------------------------------------------------

Based on the alignments of the temporary probe set to the exemplar genomes, we
need to extract FASTA data from each of the exemplar sequences so that we can
design baits targeting the conserved loci in each. This is pretty similar to
what we did earlier for the temporary probe set, except that now we're running
the extraction across all the exemplar taxa.

Before we begin, we need to make a configuration file with all the genome
locations in it (again, as before):

.. code-block:: python

    [scaffolds]
    menMol1:/path/to/uce-coleoptera/genomes/menMol1/menMol1.2bit
    agrPla1:/path/to/uce-coleoptera/genomes/agrPla1/agrPla1.2bit
    anoGla1:/path/to/uce-coleoptera/genomes/anoGla1/anoGla1.2bit
    denPon1:/path/to/uce-coleoptera/genomes/denPon1/denPon1.2bit
    lepDec1:/path/to/uce-coleoptera/genomes/lepDec1/lepDec1.2bit
    ontTau1:/path/to/uce-coleoptera/genomes/ontTau1/ontTau1.2bit
    triCas1:/path/to/uce-coleoptera/genomes/triCas1/triCas1.2bit

Using the configuration file, we need to extract the FASTA sequence that we
need from each exemplar taxon.  Here, we're buffering each locus to 180 bp to
give us a little more room to work with during the probe design step.  The
program takes our config file as input, along with the folder of lastz_ results
created above.  The `--name-pattern` argument allows us to match files int the
`--lastz` directory, `--probes` is how we buffer the sequence, and we pass the
name of an output directory to `--output`:

.. code-block:: bash

    > phyluce_probe_slice_sequence_from_genomes \
        --conf coleoptera-genome.conf \
        --lastz coleoptera-genome-lastz \
        --probes 180 \
        --name-pattern "triCas1+5.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
        --output coleoptera-genome-fasta

    2016-06-03 15:07:16,642 - Phyluce - INFO - =================== Starting Phyluce: Slice Sequence ===================
    2016-06-03 15:07:16,644 - Phyluce - INFO - ------------------- Working on menMol1 genome -------------------
    2016-06-03 15:07:16,645 - Phyluce - INFO - Reading menMol1 genome
    2016-06-03 15:07:20,221 - Phyluce - INFO - menMol1: 884 uces, 139 dupes, 745 non-dupes, 0 orient drop, 2 length drop, 738 written
    2016-06-03 15:07:20,222 - Phyluce - INFO - ------------------- Working on agrPla1 genome -------------------
    2016-06-03 15:07:20,223 - Phyluce - INFO - Reading agrPla1 genome
    2016-06-03 15:07:24,759 - Phyluce - INFO - agrPla1: 1410 uces, 184 dupes, 1226 non-dupes, 7 orient drop, 63 length drop, 1156 written
    2016-06-03 15:07:24,760 - Phyluce - INFO - ------------------- Working on anoGla1 genome -------------------
    2016-06-03 15:07:24,761 - Phyluce - INFO - Reading anoGla1 genome
    2016-06-03 15:07:29,926 - Phyluce - INFO - anoGla1: 1474 uces, 224 dupes, 1250 non-dupes, 6 orient drop, 35 length drop, 1209 written
    2016-06-03 15:07:29,926 - Phyluce - INFO - ------------------- Working on denPon1 genome -------------------
    2016-06-03 15:07:29,929 - Phyluce - INFO - Reading denPon1 genome
    2016-06-03 15:07:34,472 - Phyluce - INFO - denPon1: 1361 uces, 305 dupes, 1056 non-dupes, 6 orient drop, 30 length drop, 1020 written
    2016-06-03 15:07:34,472 - Phyluce - INFO - ------------------- Working on lepDec1 genome -------------------
    2016-06-03 15:07:34,473 - Phyluce - INFO - Reading lepDec1 genome
    2016-06-03 15:07:40,020 - Phyluce - INFO - lepDec1: 1436 uces, 259 dupes, 1177 non-dupes, 10 orient drop, 28 length drop, 1139 written
    2016-06-03 15:07:40,021 - Phyluce - INFO - ------------------- Working on ontTau1 genome -------------------
    2016-06-03 15:07:40,022 - Phyluce - INFO - Reading ontTau1 genome
    2016-06-03 15:07:44,350 - Phyluce - INFO - ontTau1: 1361 uces, 206 dupes, 1155 non-dupes, 9 orient drop, 41 length drop, 1105 written
    2016-06-03 15:07:44,350 - Phyluce - INFO - ------------------- Working on triCas1 genome -------------------
    2016-06-03 15:07:44,351 - Phyluce - INFO - Reading triCas1 genome
    2016-06-03 15:07:49,499 - Phyluce - INFO - triCas1: 1513 uces, 199 dupes, 1314 non-dupes, 26 orient drop, 46 length drop, 1242 written

If we look at out directory structure, it looks something like:

.. code-block:: text

    uce-coleoptera
    ├── alignments (collapsed)
    ├── bed (collapsed)
    ├── genomes (collapsed)
    └── probe-design
        ├── coleoptera-genome.conf
        ├── coleoptera-genome-fasta
        │   ├── agrpla1.fasta
        │   ├── anogla1.fasta
        │   ├── denpon1.fasta
        │   ├── lepdec1.fasta
        │   ├── menmol1.fasta
        │   ├── onttau1.fasta
        │   └── tricas1.fasta
        ├── coleoptera-genome-lastz
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_agrPla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_anoGla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_denPon1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_lepDec1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_menMol1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_ontTau1.lastz.clean
        │   └── triCas1+5.temp-DUPE-SCREENED.probes_v_triCas1.lastz.clean
        └── triCas1+5+menMol1.sqlite

The FASTA files we just created are in `uce-coleoptera/coleoptera-genome-fasta`.
The output from the program that you see basically shows you how many UCE loci
we extracted from each of the exemplar genomes.  As expected, the lowest number
we located and extracted are from the `menMol1` (outgroup) genome.

If we have a look in one of these FASTA files, it looks like:

.. code-block:: text

    > less probe-design/coleoptera-genome-fasta/agrpla1.fasta

    >slice_0|contig:KL218988.1|slice:301686-301866|uce:uce-500|match:301721-301831|orient:+|probes:1
    TCGAACTTCTGGTGCTTGTCACCCTTGATGTCCGCACCGAATTCCTTCACGAGACTGTTT
    CTAAAACTTTGGACAAGATGATTGGAGACACAGAAACAGACGAACATCAGAAACGTGTAT
    ATCTTCACTTTCTAGAGCATTCATACAAACTTATTACCAGATGTACTCAGCAGCAGCTTT
    >slice_1|contig:KL218988.1|slice:319624-319804|uce:uce-501|match:319637-319791|orient:+|probes:2
    atgattttttcaaAGGTTACAGCGAAGTCCTCGATTCTACAGCAGATCGAAGAACTAGGA
    GAAGAGACTGGCCTGGTGTGCTGTATTTGTCGCGAGGGATACAAGTATCAACCTGCCAAG
    GTATTGGGAATTTATACGTTTACAAAGAGGTGCAACGTGGACGAGTTCGAAGCAAAACCA
    >slice_2|contig:KL219144.1|slice:184423-184603|uce:uce-503|match:184453-184573|orient:+|probes:1
    agttttaaataatcttACCTAAAGAACTAAAATGAAGAAGCATTTCGTCTGCTCGTAAGT
    CTTGAGCAATGACATCATCAGCGTAGACACATATTATAGCAAGGCATATAAATAGGTGAA
    AGTAGTCTGTTAGATAATTTGCCCAACAAGCTTCCCACAGTCTAAGGGCAACACCTTCGG
    >slice_3|contig:KL219144.1|slice:237138-237318|uce:uce-504|match:237150-237306|orient:+|probes:2
    CTGTGCAAGAGTGAGTGCCATTGATGCAACACTTGAGCGAGATGATCTAAACCTCCATGG
    TGAAAATGAAGAATTTTATATTGAGATTCCCTCGAAGCAACAACCACCTGCCCTGATGTG
    CAGCTTGAGTCGTTAAAGAAAAGCCTTAAAGATCTCATTTGGCTTAGATCAACGCTGAAC

Find which loci we detect consistently
--------------------------------------

As before, we want to determine which loci we are detecting consistently across
all of the exemplar taxa when doing these in-silico searches.  To do that, we'll
run another bit of python code.  Here, we're working in the
`uce-coleoptera/coleoptera-genome-lastz` directory.  This program will create a
relational database that houses detections of loci in the exemplar taxa.  It
takes, as input, the folder of FASTAs we just created, the base genome taxon,
and a name to use as the output database:

.. code-block:: python

    > phyluce_probe_get_multi_fasta_table \
        --fastas coleoptera-genome-fasta \
        --output multifastas.sqlite \
        --base-taxon triCas1

    menmol1.
    agrpla1..
    anogla1..
    denpon1..
    lepdec1..
    onttau1..
    tricas1..
    Creating database
    Inserting results

If we take a look at the table contents in the database (see :ref:`Database`
for more instructions on sqlite_ databases), we see something like the
following:

.. code-block:: text

    locus       menmol1     agrpla1     anogla1     denpon1     lepdec1     onttau1     tricas1
    ----------  ----------  ----------  ----------  ----------  ----------  ----------  ----------
    uce-500     0           1           1           0           1           1           1
    uce-501     1           1           1           1           1           1           1
    uce-503     1           1           1           1           0           1           1
    uce-504     1           1           1           1           1           1           1
    uce-505     1           0           1           0           0           1           1
    uce-506     1           1           1           1           1           1           1
    uce-507     0           1           1           1           1           1           1
    uce-508     1           1           1           1           1           0           1
    uce-509     1           0           1           1           1           1           1
    uce-967     0           0           0           1           1           1           0

Which shows our detection of conserved loci in each of the exemplar taxa when we
search for them using the temporary probes that we designed from the base
genome.  As before, we can get some idea of the distribution of hits among
exemplar taxa (e.g., are loci detected in "all", n-1 taxa, n-2 taxa, etc.).

.. code-block:: python

    > phyluce_probe_query_multi_fasta_table \
        --db multifastas.sqlite \
        --base-taxon triCas1

    Loci shared by 0 taxa:  1,437.0
    Loci shared by 1 taxa:  1,437.0
    Loci shared by 2 taxa:  1,355.0
    Loci shared by 3 taxa:  1,303.0
    Loci shared by 4 taxa:  1,209.0
    Loci shared by 5 taxa:  1,099.0
    Loci shared by 6 taxa:  820.0
    Loci shared by 7 taxa:  386.0

Again, we've got to make a decision here about how conservative we want to be
regarding baits that hit all/some taxa.  We only get 386 loci that we detect
in all exemplars (including the base genome and the `menMol1` outgroup).  That
seems too strict (particularly because this total includes `menMol1`, which is
really divergent from our taxa of interest).  We also have to keep in mind that
we can randomly fail to detect loci that are actually present, either by chance
or do to sequence divergences that are >50% (the value we used in our search).
In the end, I settled on loci we detected in ≥ 4 exemplar taxa.  So we need to
extract those from the database and store them in `triCas1+5-back-to-4.conf`:

.. code-block:: python

    phyluce_probe_query_multi_fasta_table \
        --db multifastas.sqlite \
        --base-taxon triCas1 \
        --output triCas1+5-back-to-4.conf \
        --specific-counts 4

    Counter({'tricas1': 1160, 'anogla1': 1140, 'agrpla1': 1091, 'lepdec1': 1080, 'onttau1': 1043, 'denpon1': 969, 'menmol1': 658})
    Total loci = 1209

The values above show the number of loci detected in each exemplar taxon and the
total number of loci we'll be targeting with the bait set we're about to design.

Your directory should look something like the following:

.. code-block:: text

    uce-coleoptera
    ├── alignments (collapsed)
    ├── bed  (collapsed)
    ├── genomes  (collapsed)
    └── probe-design
        ├── coleoptera-genome.conf
        ├── coleoptera-genome-fasta
        │   ├── agrpla1.fasta
        │   ├── anogla1.fasta
        │   ├── denpon1.fasta
        │   ├── lepdec1.fasta
        │   ├── menmol1.fasta
        │   ├── onttau1.fasta
        │   └── tricas1.fasta
        ├── coleoptera-genome-lastz
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_agrPla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_anoGla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_denPon1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_lepDec1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_menMol1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_ontTau1.lastz.clean
        │   └── triCas1+5.temp-DUPE-SCREENED.probes_v_triCas1.lastz.clean
        ├── multifastas.sqlite
        ├── triCas1+5-back-to-4.conf
        ├── triCas1+5-back-to-4.conf.missing.matrix
        └── triCas1+5+menMol1.sqlite


Final bait set design
=====================

Design a bait set using all exemplar genomes (and the base)
-----------------------------------------------------------

Now that we've settled on the set of loci we'll try to enrich, we want to
design baits to target them.  In contrast to the steps we took before to design
the temporary bait set, we're using all of the exemplar genomes and the base
genome to design probes.  This way, we'll have a heterogeneous bait mix that
contains probes designed from each exemplar but targeting the same locus, which
should make the probe set we're designing more "universal".

To do this, we use a program similar to what we used before, except that this
program has been modified to design probes across many exemplar genomes
(instead of just one).  As input, we give the program the name of the directory
holding all of our fastas and the name of the config file we created in the step
above.  Then, as before, we need to add some metadata that will be incorporated
to the bait set design file, and we tell the program to tile at 3x density, use
a "middle" overlap, remove baits with >25% masking, and to design two probes
targeting each locus.  Finally, we write this probe set to a file named
`coleoptera-v1-master-probe-list.fasta`.

.. code-block:: python

    phyluce_probe_get_tiled_probe_from_multiple_inputs \
        --fastas coleoptera-genome-fasta \
        --multi-fasta-output triCas1+5-back-to-4.conf \
        --probe-prefix "uce-" \
        --designer faircloth \
        --design coleoptera-v1 \
        --tiling-density 3 \
        --overlap middle \
        --masking 0.25 \
        --remove-gc \
        --two-probes \
        --output coleoptera-v1-master-probe-list.fasta

    GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGNNGGGGGGGGGGGGGGGGNGGGGGGGGGGGGNNGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGNGGGGGGGGGGGGGGGGGGGGGGGNGGGGGGNGGGGGGGGGGGGGGGGGGGGGGGGGNNGG

    Conserved locus count = 1209
    Probe Count = 14113

Note that the number of baits that we've designed to target 1209 conserved loci
is quite high - this is because we're including roughly 2 baits for 1209 loci
across 7 exemplar taxa (16926 is the theoretical maximum).

Remove duplicates from our bait set
-----------------------------------

As before, we need to check out bait set for duplicate loci.  This time, the
search is going to take longer, because of the larger number of baits.  We'll
align all the probes to themselves, then read in the alignments, and filter the
probe list to remove putative duplicates.

Align probes to themselves at low stringency to identify duplicates:

.. code-block:: python

    phyluce_probe_easy_lastz \
        --target coleoptera-v1-master-probe-list.fasta \
        --query coleoptera-v1-master-probe-list.fasta \
        --identity 50 \
        --coverage 50 \
        --output coleoptera-v1-master-probe-list-TO-SELF-PROBES.lastz

    Started:  Fri Jun 03, 2016  15:50:52
    Ended:  Fri Jun 03, 2016  15:51:11
    Time for execution:  0.322272149722 minutes

Now, screen the alignements and filter our master probe list to remove
duplicates:

.. code-block:: python

    phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
        --fasta coleoptera-v1-master-probe-list.fasta \
        --lastz coleoptera-v1-master-probe-list-TO-SELF-PROBES.lastz \
        --probe-prefix=uce-

    Parsing lastz file...
    Screening results...
    Screened 14112 fasta sequences.  Filtered 37 duplicates. Kept 13674.

The master probe list that has been filtered of putatively duplicate loci is now
located in `coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta`.

Your directory should look something like the following:

.. code-block:: text

    uce-coleoptera
    ├── alignments (collapsed)
    ├── bed  (collapsed)
    ├── genomes  (collapsed)
    └── probe-design
        ├── coleoptera-genome.conf
        ├── coleoptera-genome-fasta
        │   ├── agrpla1.fasta
        │   ├── anogla1.fasta
        │   ├── denpon1.fasta
        │   ├── lepdec1.fasta
        │   ├── menmol1.fasta
        │   ├── onttau1.fasta
        │   └── tricas1.fasta
        ├── coleoptera-genome-lastz
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_agrPla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_anoGla1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_denPon1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_lepDec1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_menMol1.lastz.clean
        │   ├── triCas1+5.temp-DUPE-SCREENED.probes_v_ontTau1.lastz.clean
        │   └── triCas1+5.temp-DUPE-SCREENED.probes_v_triCas1.lastz.clean
        ├── coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta
        ├── coleoptera-v1-master-probe-list.fasta
        ├── coleoptera-v1-master-probe-list-TO-SELF-PROBES.lastz
        ├── multifastas.sqlite
        ├── triCas1+5-back-to-4.conf
        ├── triCas1+5-back-to-4.conf.missing.matrix
        └── triCas1+5+menMol1.sqlite

The master bait list
====================

What we've created, above, is the master bait list that contains baits targeting
the conserved locus we identified.  Because we've designed probes from multiple
exemplar taxa, the number of overall baits is high (as high as 14 baits
targeting each conserved locus).  This bait set is ready for synthesis and
subsequent enrichment of these conserved loci shared among coleoptera.


Subsetting the master probe list
--------------------------------

Sometimes we might not want to synthesize *all* of the baits for *all* of
the loci.  For instance, we might be enriching loci from species that are nested
within the clade defined by (('Anoplophora glabripennis (Asian longhorned
beetle)':'Leptinotarsa decemlineata (Colorado potato beetle)')'Dendroctonus
ponderosae (mountain pine beetle)'), and because we're only working with these
species, we might want to drop the baits targeting UCE loci in Agrilus
planipennis (emerald ash borer), Tribolium castaneum (red flour beetle),
Onthophagus taurus (taurus scarab), and Mengenilla moldrzyki (Strepsiptera).
This is actually pretty easy to do - we just need to subset the baits to
include those taxa that we *do* want.  Given the example, above, we can run:

.. code-block:: python

    > phyluce_probe_get_subsets_of_tiled_probes \
        --probes coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta \
        --taxa anogla1 lepdec1 denpon1 \
        --output coleoptera-v1-master-probe-list-DUPE-SCREENED-SUBSET-CLADE_1.fasta

    All probes = 13674
    --- Probes by taxon ---
    anogla1 2189
    menmol1 1236
    lepdec1 2083
    denpon1 1867
    onttau1 1970
    agrpla1 2086
    tricas1 2243
    --- Post  filtering ---
    Conserved locus count = 1169
    Probe Count = 6139


In-silico test of the bait design
=================================

Now that we've designed our baits, it's always good to run a sanity check on the
data - if we use the baits to collect data from a selection of available genomes
(or other genetic data), can we reconstruct a phylogeny that is sane, given what
we know about the specific taxa?

How we do that is outlined below.  Many of the steps we've run before, so I'm
not going to explain these quite as much as I have previously.  Several other of
the steps that we're going to run are also outlined in :ref:`Tutorial I`.

First, we need to make a directory to hold our in-silico test results:

.. code-block:: bash

    > cd uce-coleoptera
    > mkdir probe-design-test
    > cd probe-design-test

Now, our directory tree should look something like:

.. code-block:: text

    uce-coleoptera
    ├── alignments (collapsed)
    ├── bed (collapsed)
    ├── genomes  (collapsed)
    ├── probe-design (collapsed)
    └── probe-design-test

Align our bait set to the extant genome sequences
-------------------------------------------------

Here (and if you have them), you may want to include *all* the genomic data you
have access to - particularly if you removed some taxa because they were closely
related to other taxa and designing probes from these closely related groups was
redundant.  To run these alignments (you've seen this before):

.. code-block:: python

    phyluce_probe_run_multiple_lastzs_sqlite \
        --db triCas1+5+strepsiptera-test.sqlite \
        --output coleoptera-genome-lastz \
        --probefile ../probe-design/coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta \
        --scaffoldlist agrPla1 anoGla1 denPon1 lepDec1 ontTau1 triCas1 menMol1 \
        --genome-base-path ../genomes \
        --identity 50 \
        --cores 16

.. warning:: Note that I am using 16 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Now, we need to extract fasta data for each of these loci.  This is effectively
the same as what we've done before, but notice the use of `--flank` in place of
`--probe`.  This tells the program that we want to extract larger chunks of
sequence, in thie base 400 bp to the each side of a given locus (if possible):

.. code-block:: python

    > phyluce_probe_slice_sequence_from_genomes \
        --conf coleoptera-genome.conf \
        --lastz coleoptera-genome-lastz \
        --output coleoptera-genome-fasta \
        --flank 400 \
        --name-pattern "coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta_v_{}.lastz.clean"

    2016-06-03 16:36:47,712 - Phyluce - INFO - =================== Starting Phyluce: Slice Sequence ===================
    2016-06-03 16:36:47,714 - Phyluce - INFO - ------------------- Working on menMol1 genome -------------------
    2016-06-03 16:36:47,715 - Phyluce - INFO - Reading menMol1 genome
    2016-06-03 16:36:57,447 - Phyluce - INFO - menMol1: 766 uces, 73 dupes, 693 non-dupes, 0 orient drop, 2 length drop, 691 written
    2016-06-03 16:36:57,447 - Phyluce - INFO - ------------------- Working on agrPla1 genome -------------------
    2016-06-03 16:36:57,448 - Phyluce - INFO - Reading agrPla1 genome
    2016-06-03 16:37:12,538 - Phyluce - INFO - agrPla1: 1151 uces, 81 dupes, 1070 non-dupes, 0 orient drop, 34 length drop, 1036 written
    2016-06-03 16:37:12,538 - Phyluce - INFO - ------------------- Working on anoGla1 genome -------------------
    2016-06-03 16:37:12,539 - Phyluce - INFO - Reading anoGla1 genome
    2016-06-03 16:37:29,320 - Phyluce - INFO - anoGla1: 1167 uces, 103 dupes, 1064 non-dupes, 0 orient drop, 17 length drop, 1047 written
    2016-06-03 16:37:29,321 - Phyluce - INFO - ------------------- Working on denPon1 genome -------------------
    2016-06-03 16:37:29,322 - Phyluce - INFO - Reading denPon1 genome
    2016-06-03 16:37:44,958 - Phyluce - INFO - denPon1: 1126 uces, 174 dupes, 952 non-dupes, 1 orient drop, 16 length drop, 935 written
    2016-06-03 16:37:44,959 - Phyluce - INFO - ------------------- Working on lepDec1 genome -------------------
    2016-06-03 16:37:44,959 - Phyluce - INFO - Reading lepDec1 genome
    2016-06-03 16:38:01,794 - Phyluce - INFO - lepDec1: 1156 uces, 142 dupes, 1014 non-dupes, 6 orient drop, 22 length drop, 986 written
    2016-06-03 16:38:01,794 - Phyluce - INFO - ------------------- Working on ontTau1 genome -------------------
    2016-06-03 16:38:01,796 - Phyluce - INFO - Reading ontTau1 genome
    2016-06-03 16:38:16,611 - Phyluce - INFO - ontTau1: 1134 uces, 100 dupes, 1034 non-dupes, 3 orient drop, 14 length drop, 1017 written
    2016-06-03 16:38:16,612 - Phyluce - INFO - ------------------- Working on triCas1 genome -------------------
    2016-06-03 16:38:16,613 - Phyluce - INFO - Reading triCas1 genome
    2016-06-03 16:38:32,786 - Phyluce - INFO - triCas1: 1172 uces, 70 dupes, 1102 non-dupes, 13 orient drop, 13 length drop, 1076 written

Match contigs to baits
----------------------

In the step above, we essentially extracted FASTA data for each taxon, and wrote
those out into individual FASTA files. These are the equivalent of the assembled
contigs that we use in the standard phyluce_ pipeline, so now, we're going to
use that workflow.  Note that the filtering in
`phyluce_assembly_match_contigs_to_probes` is more strict that what we used
above to identify contigs.

.. code-block:: python

    > phyluce_assembly_match_contigs_to_probes \
        --contigs coleoptera-genome-fasta \
        --probes ../probe-design/coleoptera-v1-master-probe-list-DUPE-SCREENED.fasta \
        --output in-silico-lastz \
        --min_coverage 67 \
        --log-path log

    2016-06-03 16:40:36,888 - phyluce_assembly_match_contigs_to_probes - INFO - agrpla1: 903 (87.16%) uniques of 1036 contigs, 0 dupe probe matches, 116 UCE loci removed for matching multiple contigs, 117 contigs removed for matching multiple UCE loci
    2016-06-03 16:41:06,688 - phyluce_assembly_match_contigs_to_probes - INFO - anogla1: 927 (88.54%) uniques of 1047 contigs, 0 dupe probe matches, 111 UCE loci removed for matching multiple contigs, 116 contigs removed for matching multiple UCE loci
    2016-06-03 16:41:28,524 - phyluce_assembly_match_contigs_to_probes - INFO - denpon1: 819 (87.59%) uniques of 935 contigs, 0 dupe probe matches, 85 UCE loci removed for matching multiple contigs, 89 contigs removed for matching multiple UCE loci
    2016-06-03 16:41:54,879 - phyluce_assembly_match_contigs_to_probes - INFO - lepdec1: 900 (91.28%) uniques of 986 contigs, 0 dupe probe matches, 80 UCE loci removed for matching multiple contigs, 81 contigs removed for matching multiple UCE loci
    2016-06-03 16:42:05,300 - phyluce_assembly_match_contigs_to_probes - INFO - menmol1: 527 (76.27%) uniques of 691 contigs, 0 dupe probe matches, 58 UCE loci removed for matching multiple contigs, 58 contigs removed for matching multiple UCE loci
    2016-06-03 16:42:29,353 - phyluce_assembly_match_contigs_to_probes - INFO - onttau1: 854 (83.97%) uniques of 1017 contigs, 0 dupe probe matches, 127 UCE loci removed for matching multiple contigs, 130 contigs removed for matching multiple UCE loci
    2016-06-03 16:43:01,303 - phyluce_assembly_match_contigs_to_probes - INFO - tricas1: 934 (86.80%) uniques of 1076 contigs, 0 dupe probe matches, 140 UCE loci removed for matching multiple contigs, 141 contigs removed for matching multiple UCE loci

Get match counts and extract FASTA information
----------------------------------------------

Now, we need to get the count of matches that we recovered to UCE loci in the
probe set, and extract all of the "good" loci to a monolithic FASTA (see
:ref:`Tutorial I` if this is not making sense):

.. code-block:: python

    phyluce_assembly_get_match_counts \
        --locus-db in-silico-lastz/probe.matches.sqlite \
        --taxon-list-config in-silico-coleoptera-taxon-sets.conf \
        --taxon-group 'all' \
        --output taxon-sets/insilico-incomplete/insilico-incomplete.conf \
        --log-path log \
        --incomplete-matrix

    2016-06-03 16:50:02,610 - phyluce_assembly_get_match_counts - INFO - There are 7 taxa in the taxon-group '[all]' in the config file in-silico-coleoptera-taxon-sets.conf
    2016-06-03 16:50:02,610 - phyluce_assembly_get_match_counts - INFO - Getting UCE names from database
    2016-06-03 16:50:02,617 - phyluce_assembly_get_match_counts - INFO - There are 1172 total UCE loci in the database
    2016-06-03 16:50:02,708 - phyluce_assembly_get_match_counts - INFO - Getting UCE matches by organism to generate a INCOMPLETE matrix
    2016-06-03 16:50:02,709 - phyluce_assembly_get_match_counts - INFO - There are 1093 UCE loci in an INCOMPLETE matrix
    2016-06-03 16:50:02,709 - phyluce_assembly_get_match_counts - INFO - Writing the taxa and loci in the data matrix to /nfs/data1/working/bfaircloth-insects/coleoptera/triCas1+5+strepsiptera-test/taxon-sets/insilico-incomplete/insilico-incomplete.conf

Now, extract the FASTA information for each locus into a monolithic FASTA file:

.. code-block:: python

    phyluce_assembly_get_fastas_from_match_counts \
        --contigs ../../coleoptera-genome-fasta \
        --locus-db ../../in-silico-lastz/probe.matches.sqlite \
        --match-count-output insilico-incomplete.conf \
        --output insilico-incomplete.fasta \
        --incomplete-matrix insilico-incomplete.incomplete \
        --log-path log

Align the conserved locus data
------------------------------

Now, we need to align the sequence data for each conserved locus in our data
set.  We'll do this using standard phyluce_ tools (mafft).  First, change into
the working directory:

.. code-block:: bash

    cd taxon-sets/insilico-incomplete

Now, align the sequences:

.. code-block:: python

    phyluce_align_seqcap_align \
        --fasta insilico-incomplete.fasta \
        --output mafft \
        --taxa 7 \
        --incomplete-matrix \
        --cores 12 \
        --no-trim \
        --output-format fasta \
        --log-path log

.. warning:: Note that I am using 12 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Trim the conserved locus alignments
-----------------------------------

Still following the standard phyluce_ workflow, trim the resulting alignments:

.. code-block:: python

    phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
        --alignments mafft \
        --output mafft-gblocks \
        --b1 0.5 \
        --b4 8 \
        --cores 12 \
        --log log

.. warning:: Note that I am using 12 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Remove the locus names from each alignment
------------------------------------------

And, remove the locus names from each of the resulting alignments:

.. code-block:: python

    phyluce_align_remove_locus_name_from_nexus_lines \
        --alignments mafft-gblocks \
        --output mafft-gblocks-clean \
        --cores 12 \
        --log-path log

.. warning:: Note that I am using 12 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Get stats across the aligned loci
---------------------------------

Compute stats across the alignments:

.. code-block:: python

    python ~/git/phyluce/bin/align/phyluce_align_get_align_summary_data \
        --alignments mafft-gblocks-clean \
        --cores 12 \
        --log-path log

    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - ========= Starting phyluce_align_get_align_summary_data =========
    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - Version: git 6ab3a4b
    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - Argument --alignments: triCas1+5+strepsiptera-test/taxon-sets/insilico-incomplete/mafft-gblocks-clean
    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - Argument --cores: 12
    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - Argument --input_format: nexus
    2016-06-03 16:57:11,675 - phyluce_align_get_align_summary_data - INFO - Argument --log_path: triCas1+5+strepsiptera-test/taxon-sets/insilico-incomplete/log
    2016-06-03 16:57:11,676 - phyluce_align_get_align_summary_data - INFO - Argument --output: None
    2016-06-03 16:57:11,676 - phyluce_align_get_align_summary_data - INFO - Argument --show_taxon_counts: False
    2016-06-03 16:57:11,676 - phyluce_align_get_align_summary_data - INFO - Argument --verbosity: INFO
    2016-06-03 16:57:11,676 - phyluce_align_get_align_summary_data - INFO - Getting alignment files
    2016-06-03 16:57:11,710 - phyluce_align_get_align_summary_data - INFO - Computing summary statistics using 12 cores
    2016-06-03 16:57:15,381 - phyluce_align_get_align_summary_data - INFO - ----------------------- Alignment summary -----------------------
    2016-06-03 16:57:15,382 - phyluce_align_get_align_summary_data - INFO - [Alignments] loci:      994
    2016-06-03 16:57:15,382 - phyluce_align_get_align_summary_data - INFO - [Alignments] length:    644,447
    2016-06-03 16:57:15,382 - phyluce_align_get_align_summary_data - INFO - [Alignments] mean:      648.34
    2016-06-03 16:57:15,383 - phyluce_align_get_align_summary_data - INFO - [Alignments] 95% CI:    9.39
    2016-06-03 16:57:15,383 - phyluce_align_get_align_summary_data - INFO - [Alignments] min:       240
    2016-06-03 16:57:15,383 - phyluce_align_get_align_summary_data - INFO - [Alignments] max:       1,444
    2016-06-03 16:57:15,383 - phyluce_align_get_align_summary_data - INFO - ------------------- Informative Sites summary -------------------
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] loci:   994
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] total:  169,024
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] mean:   170.04
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] 95% CI: 4.48
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] min:    0
    2016-06-03 16:57:15,384 - phyluce_align_get_align_summary_data - INFO - [Sites] max:    390
    2016-06-03 16:57:15,386 - phyluce_align_get_align_summary_data - INFO - ------------------------- Taxon summary -------------------------
    2016-06-03 16:57:15,386 - phyluce_align_get_align_summary_data - INFO - [Taxa] mean:            5.76
    2016-06-03 16:57:15,386 - phyluce_align_get_align_summary_data - INFO - [Taxa] 95% CI:          0.07
    2016-06-03 16:57:15,386 - phyluce_align_get_align_summary_data - INFO - [Taxa] min:             3
    2016-06-03 16:57:15,386 - phyluce_align_get_align_summary_data - INFO - [Taxa] max:             7
    2016-06-03 16:57:15,387 - phyluce_align_get_align_summary_data - INFO - ----------------- Missing data from trim summary ----------------
    2016-06-03 16:57:15,387 - phyluce_align_get_align_summary_data - INFO - [Missing] mean: 0.00
    2016-06-03 16:57:15,387 - phyluce_align_get_align_summary_data - INFO - [Missing] 95% CI:       0.00
    2016-06-03 16:57:15,387 - phyluce_align_get_align_summary_data - INFO - [Missing] min:  0.00
    2016-06-03 16:57:15,387 - phyluce_align_get_align_summary_data - INFO - [Missing] max:  0.00
    2016-06-03 16:57:15,399 - phyluce_align_get_align_summary_data - INFO - -------------------- Character count summary --------------------
    2016-06-03 16:57:15,399 - phyluce_align_get_align_summary_data - INFO - [All characters]        3,655,040
    2016-06-03 16:57:15,399 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]           3,518,743
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]            946 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]            946 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]            865 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]            865 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]            865 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]            275 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]            275 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2016-06-03 16:57:15,399 - phyluce_align_get_align_summary_data - INFO - [All characters]        3,655,040
    2016-06-03 16:57:15,399 - phyluce_align_get_align_summary_data - INFO - [Nucleotides]           3,518,743
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - ---------------- Data matrix completeness summary ---------------
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 50%]            946 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 55%]            946 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 60%]            865 alignments
    2016-06-03 16:57:15,400 - phyluce_align_get_align_summary_data - INFO - [Matrix 65%]            865 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 70%]            865 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 75%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 80%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 85%]            657 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 90%]            275 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Matrix 95%]            275 alignments
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - ------------------------ Character counts -----------------------
    2016-06-03 16:57:15,401 - phyluce_align_get_align_summary_data - INFO - [Characters] '-' is present 136,297 times
    2016-06-03 16:57:15,402 - phyluce_align_get_align_summary_data - INFO - [Characters] 'A' is present 1,047,965 times
    2016-06-03 16:57:15,402 - phyluce_align_get_align_summary_data - INFO - [Characters] 'C' is present 708,469 times
    2016-06-03 16:57:15,402 - phyluce_align_get_align_summary_data - INFO - [Characters] 'G' is present 706,567 times
    2016-06-03 16:57:15,402 - phyluce_align_get_align_summary_data - INFO - [Characters] 'T' is present 1,055,742 times
    2016-06-03 16:57:15,402 - phyluce_align_get_align_summary_data - INFO - ========= Completed phyluce_align_get_align_summary_data ========


.. warning:: Note that I am using 12 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Generate an incomplete matrix
-----------------------------

Now, given the alignments that we have, let's generate a 70% complete matrix

.. code-block:: python

    phyluce_align_get_only_loci_with_min_taxa \
        --alignments mafft-gblocks-clean \
        --taxa 7 \
        --output mafft-gblocks-70p \
        --percent 0.75 \
        --cores 12 \
        --log log

    2016-06-03 16:58:19,687 - phyluce_align_get_only_loci_with_min_taxa - INFO - Copied 865 alignments of 994 total containing ≥ 0.75 proportion of taxa (n = 5)

.. warning:: Note that I am using 12 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.


Prep raxml files, run raxml ML searches, and reconcile best tree w/ bootreps
----------------------------------------------------------------------------

Setup the PHYLIP-formatted files for raxml:

.. code-block:: python

    phyluce_align_format_nexus_files_for_raxml \
        --alignments mafft-gblocks-70p \
        --output mafft-gblocks-70p-raxml \
        --log-path log --charsets

Now, run raxml against this phylip file

.. code-block:: python

     raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -N 20 -p 772374015 -n BEST -s mafft-gblocks-70p.phylip -o menmol1 -T 10
     raxml/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -N autoMRE -p 772374015 -b 444353738 -n bootrep -s mafft-gblocks-70p.phylip -o menmol1 -T 10

.. warning:: Note that I am using 10 physical CPU cores (`--cores`) to do this
    work. You need to use the number of physical cores available on *your*
    machine.

Now, reconcile the best ML tree w/ the bootreps:

.. code-block:: python

    raxmlHPC-SSE3 -f b \
        -m GTRGAMMA \
        -t RAxML_bestTree.BEST \
        -z RAxML_bootstrap.bootrep \
        -n FINAL -o menmol1

And rename the tips.  To do this, setup a config file with the old and new
names, like:

.. code-block:: text

    [all]
    agrpla1:Agrilus planipennis (emerald ash borer)
    anogla1:Anoplophora glabripennis (Asian longhorned beetle)
    denpon1:Dendroctonus ponderosae (mountain pine beetle)
    lepdec1:Leptinotarsa decemlineata (Colorado potato beetle)
    menmol1:Mengenilla moldrzyki (Strepsiptera)
    onttau1:Onthophagus taurus (taurus scarab)
    tricas1:Tribolium castaneum (red flour beetle)

And rename the tips:

.. code-block:: python

    phyluce_genetrees_rename_tree_leaves \
        --order left:right \
        --input-format newick \
        --output-format newick \
        --config rename.conf \
        --section all \
        --input RAxML_bipartitions.FINAL \
        --output RAxML_bipartitions.NAME.FINAL.tre
