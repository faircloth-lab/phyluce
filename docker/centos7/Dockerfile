FROM centos:7
MAINTAINER "Brant Faircloth" <brant _at_ faircloth-lab _dot_ org>
ENV container docker

# update yum
RUN yum -y update
RUN yum clean all
RUN yum -y install bzip2 git ssh

# download conda
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	&& bash /tmp/miniconda.sh -bfp /usr/local/ \
	&& rm -rf /tmp/miniconda.sh

RUN yum clean all
RUN rm -rf /var/cache/yum

# update conda
RUN conda update conda

# set conda REPOS
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# get set of dependencies for phyluce
RUN conda install nomkl abyss=1.5.2 bcftools bedtools biopython bx-python dendropy gblocks lastz mafft muscle pandas picard pysam pyvcf raxml samtools seqtk trimal trinity velvet spades
RUN conda clean --all --yes