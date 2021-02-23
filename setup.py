# Always prefer setuptools over distutils
import os
import glob
import pathlib
from setuptools import setup, find_packages
from collections import defaultdict

from phyluce import __version__

# import pdb

# Copy the executable files to a single place.  Used below.
scrpt = []
for d in [
    "align",
    "assembly",
    "genetrees",
    "ncbi",
    "probes",
    "utilities",
    "workflow",
]:
    scrpt.extend([i for i in glob.glob(os.path.join("bin", d, "*"))])


"""
# the following function and call are used to grab the data files
# that we need for workflows and tests

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        # skip hidden files and dirs
        filenames = [f for f in filenames if not f[0] == "."]
        directories[:] = [d for d in directories if not d[0] == "."]
        for filename in filenames:
            paths.append(os.path.join("..", path, filename))
    return paths


t_files = []
for dir in ["phyluce/tests/test-expected"]:
    t_files.extend(package_files(dir))
"""


# the following function and call are used to grab the data files
# that we need for workflows and tests
def data_files(directory):
    paths = defaultdict(list)
    for (path, directories, filenames) in os.walk(directory):
        # skip hidden files and dirs
        filenames = [f for f in filenames if not f[0] == "."]
        directories[:] = [d for d in directories if not d[0] == "."]
        for filename in filenames:
            # pdb.set_trace()
            # paths.append([path,os.path.join(path, filename)])
            paths[path].append(os.path.join(path, filename))
    return paths


d_files = []
for dir in ["config", "workflows", "phyluce/tests/test-expected"]:
    files = data_files(dir)
    for k, v in files.items():
        if k.startswith("config") or k.startswith("workflows"):
            k = os.path.join("phyluce", k)
        d_files.append((k, v))

setup(
    name="phyluce",  # Required
    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,  # Required
    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description="Software for UCE (and general) phylogenomics",  # Optional
    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    #
    # Often, this is the same as your README, so you can just read it in from
    # that file directly (as we have already done above)
    #
    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    # long_description=long_description,  # Optional
    # Denotes that our long_description is in Markdown; valid values are
    # text/plain, text/x-rst, and text/markdown
    #
    # Optional if long_description is written in reStructuredText (rst) but
    # required for plain-text or Markdown; if unspecified, "applications should
    # attempt to render [the long_description] as text/x-rst; charset=UTF-8 and
    # fall back to text/plain if it is not valid rst" (see link below)
    #
    # This field corresponds to the "Description-Content-Type" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    long_description_content_type="text/x-rst",  # Optional (see note above)
    # This should be a valid link to your project's main homepage.
    #
    # This field corresponds to the "Home-Page" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    url="https://github.com/faircloth-lab/phyluce",  # Optional
    # This should be your name or the name of the organization which owns the
    # project.
    author="Brant C. Faircloth",  # Optional
    # This should be a valid email address corresponding to the author listed
    # above.
    author_email="borg@faircloth-lab.org",  # Optional
    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        # Pick your license as you wish
        "License :: OSI Approved :: BSD License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3 :: Only",
    ],
    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a list of additional keywords, separated
    # by commas, to be used to assist searching for the distribution in a
    # larger catalog.
    keywords="uce, phylogenomics, phylogenetics, genetics",  # Optional
    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    # package_dir={'': 'phyluce',},  # Optional
    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=find_packages(),  # Required
    # Specify which Python versions you support. In contrast to the
    # 'Programming Language' classifiers above, 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires="~=3.6",
    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # install_requires=['peppercorn'],  # Optional
    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    # extras_require={  # Optional
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    # },
    #
    # If there are data files included in your packages that need to be
    # installed, specify them here.
    # package_data={"": t_files},
    data_files=d_files,
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # `pip` to create the appropriate form of executable for the target
    # platform.
    #
    # For example, the following would provide a command called `sample` which
    # executes the function `main` from this package when invoked:
    # entry_points={  # Optional
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    # },
    scripts=scrpt,
    # List additional URLs that are relevant to your project as a dict.
    #
    # This field corresponds to the "Project-URL" metadata fields:
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #
    # Examples listed include a pattern for specifying where the package tracks
    # issues, where the source is hosted, where to say thanks to the package
    # maintainers, and where to support the project financially. The key is
    # what's used to render the link text on PyPI.
    project_urls={  # Optional
        "Bug Reports": "https://github.com/faircloth-lab/phyluce/issues/",
        "Source": "https://github.com/faircloth-lab/phyluce/",
    },
)
