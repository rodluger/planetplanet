#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
lightcurve.py
-------------

Computes a light curve of the TRAPPIST-1 system over ten days, with
orbital parameters drawn at random from their prior distributions.
All transits, secondary eclipses, planet-planet occultations, and mutual
transits are shown. Click on an individual event to see it in detail.
Then click anywhere in the pop-up window to animate the event.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, phasecurve = True, airless = True)

# Get the occultation light curves for the first 10 days
time = np.linspace(0., 10., 10000)
system.compute(time)

# Plot all of the occultations
system.plot_lightcurve()
pl.show()