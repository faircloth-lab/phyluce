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
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# get set of dependencies for phyluce
RUN conda install nomkl abyss=1.5.2 bcftools bedtools biopython bowtie bwa bx-python dendropy gatk=3.8=py36_0 gblocks jellyfish=2.2.3 lastz mafft muscle pandas picard pysam pytest pyvcf raxml samtools seqtk trimal trinity=2.1.1 velvet=1.2.10 spades=3.12.0
RUN conda clean --all --yes

# get gatk, setup, and cleanup
RUN curl -sSL 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef' -o /usr/local/GATK-archive-3.8-1.tar.bz2 \
    && cd /usr/local \
    && tar -xjf GATK-archive-3.8-1.tar.bz2 \
    && gatk-register /usr/local/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
    && rm -rf /usr/local/GATK-archive-3.8-1.tar.bz2 \
    && rm -rf /usr/local/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef \
    && rm -rf /usr/local/LICENSE.txt