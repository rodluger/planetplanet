#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
broken.py
---------

Numerical issues lead to NANs in the integral
computation. Here's an example. Need to fix this.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, ttvs = False, phasecurve = False, adaptive = True)

# Get the occultation light curve
time = np.linspace(4.95, 4.96, 1000)
system.compute(time)

# Plot all of the occultations
system.plot_lightcurve()
pl.show()