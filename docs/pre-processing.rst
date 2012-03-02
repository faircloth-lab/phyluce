#####################
Data Pre-processing
#####################

**************
Demultiplexing
**************

Because we are reducing genomes, and beceause we generally reduce genomes to a
fraction of their original size, we can combine DNA libraries from many sources
into a single lane of HiSeq (or even MiSeq) sequencing.  Therefore, because we
are multiplexing these libraries prior to sequencing, we need to "demultiplex"
these libraries after they're sequenced.

Demultiplexing is a tricky process, made tricker by the platform with which you 
are working.  Thankfully, demultiplexing (or "demuxing") Illumina reads is
rather simple. There are several options:

* the Illumina pipeline
* `fastx-tools <>`_

I have also written `splitaake`_, which is a program to help demultiplex
Illumina reads.  `splitaake`_ is particularly useful when you have very many
sequence tags (aka "barcodes") to demultiplex, when your sequence tags are
long, or some combination of the two.

To demultiplex your data with `splitaake`_, you need a couple of things - first
you need to `install splitaake <>`_.  Second, you need to create a configuration
file (e.g. my-test-file.conf) for `splitaake`_ that relates your sequence tags
to their samples.  Here  is an example:

.. code-block:: python

	[L007]
	#name:sequence-tag
	BFIDT-000:AACCGAGTTA
	BFIDT-001:AATACTTCCG
	BFIDT-002:AATTAAGGCC
	BFIDT-003:AAGCTTATCC
	BFIDT-004:CTAACGATCG
	BFIDT-005:CTACAACGGC

Once all that's done, you may want to start up `splitaake`_ in a `screen`_ or 
`tmux`_ session by running:

.. code-block:: bash

	python dmux.py r1.fastq.gz r2.fastq.gz r3.fastq.gz \
			my-test-file.conf --section L007

This can take a rather long time for a single lane of 150-190 million Illumina
reads - on the order of 15-18 hours for 50-60 sequence tags.  These are ways to
speed up this process, but it's not always as easy as it seems (I/O becomes a 
problem very rapidly).

When you're done, regardless of the program that you used, you should have
something similar to a directory containing your reads data in fastq format::

	/Data
		L007/
			r1.fastq.gz
			r2.fastq.gz
			r3.fastq.gz
			dmux/
				BFIDT-000.fastq.gz
				BFIDT-001.fastq.gz
				BFIDT-002.fastq.gz
				BFIDT-003.fastq.gz
				BFIDT-004.fastq.gz
				BFIDT-005.fastq.gz

In the above, the reasd from the original files `r1.fastq.gz`, `r2.fastq.gz`, 
`r2.fastq.gz` have been split into their component parts - with each name in the
configuration file above becoming its own gzip file of results.


.. _splitaake: https://github.com/faircloth-lab/splitaake/
.. _screen: http://www.gnu.org/software/screen/
.. _tmux: http://tmux.sourceforge.net/
.. _gzip: http://www.gzip.org/

******************************
Adapter- and quality- trimming
******************************