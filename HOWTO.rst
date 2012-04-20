Install dependencies
####################

- Biopython
- bx-python
- numpy

- Kent tools

  ::

    git clone git://genome-source.cse.ucsc.edu/kent.git
    export MACHTYPE=x86_64
    mkdir ~/bin
    mkdir ~/bin/$MACHTYPE
    cd kent/src/lib
    make
    cd ../jkOwnLib
    make

    # kent tools need to be told the location of libpng in /usr/X11/lib

- scythe:     https://github.com/vsbuffalo/scythe

  ::

    git clone https://github.com/vsbuffalo/scythe.git
    cd scythe/
    make build
    ./scythe    # scythe should produce output here
    mv scythe /usr/local/bin

- sickle:     https://github.com/najoshi/sickle

  ::

    git clone https://github.com/najoshi/sickle.git
    cd sickle/
    make
    ./sickle    # sickle should produce output
    mv sickle /usr/local/bin

- bf-tools:   https://github.com/brantfaircloth/bf-tools

  bf-tools is now https://github.com/faircloth-lab/seqtools

  - If pip is installed

    ::

      pip install seqtools

  - If easy_install is used

    ::

      easy_install seqtools

- BioPerl http://www.bioperl.org/wiki/Installing_BioPerl
   
  Installing BioPerl

  ::

    perl -MCPAN -e shell
    # (in CPAN Shell)
    install Bundle::CPAN
    reload
    install Module::Build
    o conf prefer_installer MB
    o conf commit
    reload
    d /bioperl/
    force install C/CJ/CJFIELDS/BioPerl-1.6.1.tar.gz

  Installing BioPerl from source

  ::

    curl -O http://bioperl.org/DIST/BioPerl-1.6.1.tar.bz2
    tar -xvf BioPerl-1.6.1.tar.bz2
    cd BioPerl-1.6.1/
    perl Build.PL    # may need sudo
    ./Build test     # may need sudo
    ./Build install  # may need sudo

  Testing BioPerl

  ::

    perl -MBio::Perl -e ''   # should run without errors


- VelvetOptimiser

  Installing Velvet

  ::

    git clone git://github.com/dzerbino/velvet.git 
    cd velvet
    make
    mv velvetg /usr/local/bin
    mv velveth /usr/local/bin

  Installing VelvetOptimiser

  ::

    curl -O http://bioinformatics.net.au/VelvetOptimiser-2.2.0.tar.gz
    tar -xvf VelvetOptimiser-2.2.0.tar.gz
    cd VelvetOptimiser-2.2.0/
    mv VelvetOpt/ /Library/Perl/5.10.0/     # install module
    mv VelvetOptimiser.pl /usr/local/bin/VelvetOptimiser






