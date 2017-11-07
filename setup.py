#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from setuptools import setup, find_packages, Extension
import glob
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"
    
# Hackishly inject a constant into builtins to enable importing of the
# module in "setup" mode. Stolen from `kplr`
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__PLANETPLANET_SETUP__ = True
import planetplanet

# PLANETPLANET C EXTENSION. Borrowing heavily from REBOUND here.
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args=['-L/usr/local/lib', 
                     '-Wl,-install_name,@rpath/libppo' + suffix]
else:
    extra_link_args=['-L/usr/local/lib']
libppomodule = Extension('libppo',
                   sources = glob.glob('rebound/src/*.c') + \
                             ['progress/progress.c',
                              'planetplanet/photo/orbit.c',
                              'planetplanet/photo/eyeball.c',
                              'planetplanet/photo/ppo.c',
                             ],
                   include_dirs = ['rebound/src/', 
                                   'progress/', 
                                   'planetplanet/photo/', 
                                   '/usr/local/include'],
                   define_macros=[ ('LIBREBOUND', None) ],
                   extra_compile_args=['-Wall', '-I/usr/local/include', 
                                       '-fstrict-aliasing', '-O3', '-std=c99',
                                       '-Wno-unknown-pragmas', '-DLIBREBOUND', 
                                       '-D_GNU_SOURCE', '-fPIC'],
                   extra_link_args=extra_link_args,
                   libraries=['gsl', 'gslcblas', 'm']
                   )

long_description = \
"A photodynamical code that computes transits, eclipses, phase curves, " + \
"and planet-planet/planet-moon occultations in planetary systems."

# Setup!
setup(name = 'planetplanet',
      version = planetplanet.__version__,
      description = 'Photodynamical code for planet-planet occultations',
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
      license = 'GPL',
      packages = ['planetplanet', 'planetplanet.photo', 'planetplanet.detect'],
      install_requires = [
                          'numpy>=1.8',
                          'scipy',
                          'matplotlib',
                          'six',
                          'tqdm',
                          'astropy',
                          'numba>=0.34',
                          'pandas',
                          'rebound'
                         ],
      include_package_data = True,
      zip_safe = False,
      test_suite='nose.collector',
      tests_require=['nose'],
      ext_modules = [libppomodule],
      )