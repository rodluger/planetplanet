#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
next_occultation.py
-------------------

Compute the time of the next occultation of a given planet
and plot the light curve.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True)

# Get the next occultation of b
t = system.next_occultation(10000, system.d, occultor = system.b)

# Get the light curve around that point
time = np.linspace(t - 0.05, t + 0.05, 1000)
system.compute(time)
system.plot_lightcurve()
pl.show()