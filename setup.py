#!/usr/bin/env python3

import re
import sys
from setuptools import setup, find_packages


if sys.version_info < (3, 5, 0, 'final', 0):
    raise SystemExit('Python 3.5 or later is required!')

with open('README.md') as fd:
    long_description = fd.read()

with open('icet/__init__.py') as fd:
    lines = '\n'.join(fd.readlines())
version = re.search("__version__ = '(.*)'", lines).group(1)
maintainer = re.search("__maintainer__ = '(.*)'", lines).group(1)
url = re.search("__url__ = '(.*)'", lines).group(1)
description = re.search("__description__ = '(.*)'", lines).group(1)

name = 'icet'  # PyPI name

# Linux-distributions may want to change the name:
if 0:
    name = 'python-icet'

setup(name=name,
      version=version,
      description=description,
      url=url,
      maintainer=maintainer,
      platforms=['unix'],
      install_requires=['numpy', 'ase', 'scipy', 'sklearn', 'pandas'],
      packages=find_packages(),
      long_description=long_description,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Physics'])
