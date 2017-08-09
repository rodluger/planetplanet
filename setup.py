#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from setuptools import setup, find_packages, Extension
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

# REBOUND C EXTENSION
if sys.platform == 'darwin':
  from distutils import sysconfig
  vars = sysconfig.get_config_vars()
  vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
  extra_link_args=['-Wl,-install_name,@rpath/libppo'+suffix]
else:
  extra_link_args=[]
libppomodule = Extension('libppo',
                   sources = ['planetplanet/photo/rebound/rebound.c',
                              'planetplanet/photo/rebound/integrator_ias15.c',
                              'planetplanet/photo/rebound/integrator_whfast.c',
                              'planetplanet/photo/rebound/integrator_whfasthelio.c',
                              'planetplanet/photo/rebound/integrator_hermes.c',
                              'planetplanet/photo/rebound/integrator_leapfrog.c',
                              'planetplanet/photo/rebound/integrator_janus.c',
                              'planetplanet/photo/rebound/integrator_sei.c',
                              'planetplanet/photo/rebound/integrator.c',
                              'planetplanet/photo/rebound/gravity.c',
                              'planetplanet/photo/rebound/boundary.c',
                              'planetplanet/photo/rebound/display.c',
                              'planetplanet/photo/rebound/collision.c',
                              'planetplanet/photo/rebound/tools.c',
                              'planetplanet/photo/rebound/derivatives.c',
                              'planetplanet/photo/rebound/tree.c',
                              'planetplanet/photo/rebound/particle.c',
                              'planetplanet/photo/rebound/output.c',
                              'planetplanet/photo/rebound/input.c',
                              'planetplanet/photo/rebound/simulationarchive.c',
                              'planetplanet/photo/rebound/transformations.c',
                              'planetplanet/photo/orbit.c',
                              'planetplanet/photo/complex.c',
                              'planetplanet/photo/roots.c',
                              'planetplanet/photo/eyeball.c',
                              'planetplanet/photo/ppo.c',
                              'planetplanet/photo/progress.c'
                              ],
                   include_dirs = ['planetplanet/photo/rebound', 'planetplanet/photo/'],
                   define_macros=[ ('LIBREBOUND', None) ],
                   extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-Wno-unknown-pragmas', '-DLIBREBOUND', '-D_GNU_SOURCE', '-fPIC'],
                   extra_link_args=extra_link_args,
                   )

long_description = \
"""
A photodynamical code that computes transits, eclipses, phase curves, and planet-planet/planet-moon occultations in planetary systems.
"""

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
      tests_require=['nose'],
      ext_modules = [libppomodule],
      )