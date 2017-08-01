#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
orbit_viewer.py
---------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball
import matplotlib.pyplot as pl

eyeball.DrawOrbit(inc = 70., Omega = 30.)
pl.show()