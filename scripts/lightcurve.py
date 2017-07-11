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
from matplotlib.ticker import MaxNLocator
import numpy as np

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, phasecurve = True, airless = True, nbody = True, seed = 999)

# Get the occultation light curves over 10 random days
tstart = np.random.random() * 10000
time = np.linspace(tstart, tstart + 10., 10000)
system.compute(time)

# Plot all of the occultations
for body in system.bodies:
  body.time -= body.time[0]
  body.time_hr -= body.time_hr[0]
fig, ax, ann = system.plot_lightcurve()

# Manually tweak the positions of some of the annotations for visibility
ann[18].set_position((7, -5))
ann[19].set_position((7, -5))
ann[29].set_position((0, -20))
ann[38].set_position((-3, -7))
ann[43].set_position((0,-25))
ann[48].set_position((5, -15))
ann[49].set_position((-5, -15))
ann[50].set_position((3, -15))
ann[51].set_position((0, -32))
ann[52].set_position((0, -32))

# Appearance
ax.set_xlim(0, 10.)
ax.set_ylim(0.9819, 1.005)
for body in system.bodies:
  ax.plot([-1.001e3, -1.000e3], [1, 1], color = body.color, label = body.name, lw = 4)
ax.annotate("Occultations by", xy = (0.175, 0.92), ha = 'left', va = 'bottom', fontweight = 'bold', fontsize = 10, xycoords = "axes fraction")
ax.legend(loc = (0.325, 0.9), ncol = 8, frameon = False, handlelength = 1)
ax.set_xlabel("Time [days]", fontweight = 'bold', fontsize = 16)
ax.set_ylabel("Normalized Flux", fontweight = 'bold', fontsize = 16)
for tick in ax.get_xticklabels() + ax.get_yticklabels():
  tick.set_fontsize(12)
ax.get_xaxis().set_major_locator(MaxNLocator(10))

fig.savefig("../img/lightcurve.pdf", bbox_inches = 'tight')
pl.show()
