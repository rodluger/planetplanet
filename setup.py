#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from setuptools import setup, find_packages

# Hackishly inject a constant into builtins to enable importing of the
# module in "setup" mode. Stolen from `kplr`
import sys
if sys.version_info[0] < 3:
  import __builtin__ as builtins
else:
  import builtins
builtins.__PLANETPLANET_SETUP__ = True
import planetplanet

long_description = \
"""
Planet-planet occultations.
"""

# Setup!
setup(name = 'planetplanet',
      version = planetplanet.__version__,
      description = 'Planet-planet occultations',
      long_description = long_description,
      classifiers = [
                      'Development Status :: 4 - Beta',
                      'License :: OSI Approved :: MIT License',
                      'Programming Language :: Python',
                      'Programming Language :: Python :: 3',
                      'Topic :: Scientific/Engineering :: Astronomy',
                    ],
      url = 'http://github.com/rodluger/planetplanet',
      author = 'Rodrigo Luger',
      author_email = 'rodluger@uw.edu',
      license = 'MIT',
      packages = ['planetplanet', 'planetplanet.photo', 'planetplanet.detect'],
      install_requires = [
                          'numpy>=1.8',
                          'scipy',
                          'matplotlib',
                          'tqdm',
                          'astropy'
                         ],
      include_package_data = True,
      zip_safe = False,
      test_suite='nose.collector',
      tests_require=['nose']
      )