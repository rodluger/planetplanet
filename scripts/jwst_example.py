#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
jwst_example.py
---------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, oversample = 1, airless = True)

# Get the next occultation
t = system.next_occultation(1000, system.c, occultor = system.b)
time = np.arange(t - 0.1, t + 0.1, 5 / 1440.)
system.compute(time)

system.plot_lightcurve(15.)
system.observe()
pl.show()