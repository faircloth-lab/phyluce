import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup, find_packages
from glob import glob


if __name__ == '__main__':
    setup(
        name='phyluce',
        version="0.1dev",
        description="phylogenetic inference from ultraconserved elements",
        author="Brant Faircloth",
        url="http://github.com/faircloth-lab/phyluce/",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
             ],
        long_description=open('README.rst').read(),
        packages=['phyluce',],
        scripts=glob('bin/*/*.py')
        )
