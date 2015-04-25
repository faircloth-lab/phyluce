FROM centos:6
MAINTAINER "Brant Faircloth" <brant _at_ faircloth-lab _dot_ org>
ENV container docker

# update yum
RUN yum -y update; yum clean all
# install wget
RUN yum -y install wget
RUN yum -y install tar
RUN yum -y install bzip2
RUN yum -y install git

# add test user
RUN useradd -ms /bin/bash test

# switch to test user
USER test
ENV HOME /home/test
WORKDIR /home/test/

# download conda
RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O $HOME/miniconda.sh"]
RUN chmod 0755 $HOME/miniconda.sh
RUN ["/bin/bash", "-c", "$HOME/miniconda.sh -b -p $HOME/conda"]
ENV PATH="$HOME/conda/bin:$PATH"
RUN rm $HOME/miniconda.sh

# update conda
RUN conda update conda

# install phyluce
RUN conda config --add channels http://conda.binstar.org/faircloth-lab
RUN conda install phyluce

# install ipython
RUN conda install ipython

# clone phyluce source to $HOME/git/phyluce
RUN mkdir git
RUN cd $HOME/git && git clone https://github.com/faircloth-lab/phyluce.git
RUN cd $HOME/git/phyluce && git fetch --all && git checkout -b working origin/working
