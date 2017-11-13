#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import os, subprocess

# Version number
__version__ = "1.0.0"
__author__ = "Rodrigo Luger (rodluger@uw.edu)"
__copyright__ = "Copyright 2017 Rodrigo Luger"
__url__ = "https://raw.githubusercontent.com/rodluger/planetplanet/master/planetplanet/__init__.py"

# Was everest imported from setup.py?
try:
    __PLANETPLANET_SETUP__
except NameError:
    __PLANETPLANET_SETUP__ = False

if not __PLANETPLANET_SETUP__:

    # Import stuff
    from . import photo, detect, constants, pool
    from .version import VersionCheck
    from .photo import *
    from .photo.trappist1 import *
    from .detect import *

    # Check for updates?
    if not os.environ.get('PLANETPLANET_SKIP_VERSION_CHECK'):
        VersionCheck()
