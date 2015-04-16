#!/usr/bin/env python
# encoding: utf-8

import glob
from distutils.core import setup

align = ["{}".format(i) for i in glob.glob("bin/align/*")]
assembly = ["{}".format(i) for i in glob.glob("bin/assembly/*")]
genetrees = ["{}".format(i) for i in glob.glob("bin/genetrees/*")]
ncbi = ["{}".format(i) for i in glob.glob("bin/ncbi/*")]
probes = ["{}".format(i) for i in glob.glob("bin/probes/*")]
snps = ["{}".format(i) for i in glob.glob("bin/snps/*")]
utilities = ["{}".format(i) for i in glob.glob("bin/utilities/*")]

scrpt = []
scrpt.extend(align)
scrpt.extend(assembly)
scrpt.extend(genetrees)
scrpt.extend(ncbi)
scrpt.extend(probes)
scrpt.extend(snps)
scrpt.extend(utilities)

setup(
    name='phyluce',
    version='1.5.0',
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
