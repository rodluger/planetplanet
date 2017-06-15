#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
binning.py
----------

Demonstrates how oversampling the transit/occultation is necessary
when the exposure time is a non-negligible fraction of the event duration.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [1.])

# Planet b
b = Planet('b', m = 1, per = 3, inc = 90., r = 3., trn0 = 0, 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# System
system = System(star, b)

# Infinitesimal exposure time
time_u = np.linspace(-0.03, 0.03, 1000)
system.oversample = 1
system.compute(time_u)
flux_u = np.array(system.flux[:,0])
flux_u /= flux_u[0]

# 5 min exposure, 10x binning
time_b = np.arange(-0.03, 0.03, 5. / 1440.)
system.oversample = 10
system.compute(time_b)
flux_b10 = np.array(system.flux[:,0])
flux_b10 /= flux_b10[0]

# 5 min exposure, 30x binning
system.oversample = 30
system.compute(time_b)
flux_b30 = np.array(system.flux[:,0])
flux_b30 /= flux_b30[0]

# 5 min exposure, 100x binning
system.oversample = 100
system.compute(time_b)
flux_b100 = np.array(system.flux[:,0])
flux_b100 /= flux_b100[0]

# Plot!
pl.plot(time_u, flux_u, label = 'unbinned')
pl.plot(time_b, flux_b10, '.', label = '10x binning')
pl.plot(time_b, flux_b30, '.', label = '30x binning')
pl.plot(time_b, flux_b100, '.', label = '100x binning')
pl.legend()
pl.xlabel('Time [days]', fontweight = 'bold', fontsize = 16)
pl.ylabel('Normalized Flux', fontweight = 'bold', fontsize = 16)
pl.show()