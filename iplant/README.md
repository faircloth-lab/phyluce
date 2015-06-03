# iPlant phyluce image setup

* start image from iplant

## Setup

### as root

* install zsh and git (git needed for antigen.zsh)

        yum install zsh
        yum install git

* edit /etc/default/useradd to make zsh the default shell for new users

        # useradd defaults file
        GROUP=100
        HOME=/home
        INACTIVE=-1
        EXPIRE=
        SHELL=/bin/zsh
        SKEL=/etc/skel
        CREATE_MAIL_SPOOL=yes

* clone antigen to /usr/local/antigen
* edit /etc/zshrc for all users to add:

        source /usr/local/antigen/antigen.zsh
        # Load the oh-my-zsh's library.
        antigen use oh-my-zsh

        # Bundles from the default repo (robbyrussell's oh-my-zsh).
        antigen bundle command-not-found

        # Syntax highlighting bundle.
        antigen bundle zsh-users/zsh-syntax-highlighting

        # ZSH port of Fish shell's history search feature.
        antigen bundle zsh-users/zsh-history-substring-search

        # Load the theme.
        antigen theme steeef

        # Tell antigen that you're done.
        antigen apply

* installed anaconda into /usr/local/anaconda
* change permissions for anaconda dir

        chown -R root:users /usr/local/anaconda
        chmod -R 0775 /usr/local/anaconda

* edited /etc/profile.d to add phyluce.sh, putting anaconda at front of path for zsh and bash

        export PATH=/usr/local/anaconda/bin:$PATH

* edited /etc/profile.d to add phyluce.csh, putting anaconda at front of path for tcsh

        setenv PATH /usr/local/anaconda/bin:$PATH

* add faircloth-lab conda repo to systemwide condarc:

        conda config --system --add channels https://conda.binstar.org/faircloth-lab

### as user

* install phyluce

        conda install phyluce

* run zsh (optional) - should start automagically once image created
    
        zsh

## Enhancements

### ipython notebook

* start the server on iplant:

        ipython notebook --no-browser --port=8889

* forward over ssh on your local machine, enter password when prompted

        ssh -N -f -L localhost:8888:localhost:8889 username@xx.yy.zz.qq

* browse to (on your local machine):

        http://localhost:8888
