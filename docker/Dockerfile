FROM ubuntu:20.04
MAINTAINER "Brant Faircloth" <brant _at_ faircloth-lab _dot_ org>
ENV container docker
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# update yum
RUN apt-get update -q && \
    apt-get install -q -y \
        bzip2 \
        ca-certificates \
        git \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mercurial \
        subversion \
        wget \
    && apt-get clean


# download conda - code from https://github.com/ContinuumIO/docker-images/blob/master/miniconda3/debian/Dockerfile
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    echo "${CONDA_MD5}  miniconda.sh" > miniconda.md5 && \
    if ! md5sum --status -c miniconda.md5; then exit 1; fi && \
    mkdir -p /opt && \
    sh miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh miniconda.md5 && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    # echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

# update conda
RUN conda update conda

# set conda REPOS
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# get set of dependencies for phyluce
RUN conda clean --all --yes

# get gatk, setup, and cleanup
RUN wget --quiet https://raw.githubusercontent.com/faircloth-lab/phyluce/devel/distrib/phyluce-py36-Linux-conda.yml &&
    conda env create -n phyluce --file phyluce-py36-Linux-conda.yml &&
    echo "conda activate phyluce" >> ~/.bashrc