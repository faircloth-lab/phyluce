#!/usr/bin/env python
# encoding: utf-8

import glob
from distutils.core import setup

align = ["{}".format(i) for i in glob.glob("bin/align/*.py")]
assembly = ["{}".format(i) for i in glob.glob("bin/assembly/*.py")]
genetrees = ["{}".format(i) for i in glob.glob("bin/genetrees/*.py")]

scrpt = []
scrpt.extend(align)
scrpt.extend(assembly)
scrpt.extend(genetrees)

setup(
    name='phyluce',
    version='2.0.0',
    description='software for UCE (and general) phylogenomics',
    url='https://github.com/faircloth-lab/phyluce',
    author='Brant C. Faircloth',
    author_email='borg@faircloth-lab.org',
    license='BSD',
    platforms='any',
    packages=[
        'phyluce',
    ],
    data_files=[('config', ['config/phyluce.conf'])],
    scripts=scrpt,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    )
