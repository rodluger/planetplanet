#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import os, subprocess

# Version number
__version__ = "0.0.1"
__author__ = "Rodrigo Luger (rodluger@uw.edu)"
__copyright__ = "Copyright 2017 Rodrigo Luger"

# Was everest imported from setup.py?
try:
  __PLANETPLANET_SETUP__
except NameError:
  __PLANETPLANET_SETUP__ = False

if not __PLANETPLANET_SETUP__:
  
  # Check if the code needs to be compiled
  if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'photo', 'libppo.so')):
    subprocess.call(["make"], cwd = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'photo'))
  
  # Import stuff
  from . import photo, detect
  from .photo import *

