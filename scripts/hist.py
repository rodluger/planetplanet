#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py
-------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np

# Instantiate the Trappist-1 system
system = Trappist1(uncertainty = True, ttvs = False)
system.settings.dt = 0.001
system.scatter_plot(0, 365 * 3)
pl.show()