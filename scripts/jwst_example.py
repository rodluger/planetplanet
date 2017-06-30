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
PARSEC = 3.086e16

def retrograde_bc():
  '''
  A retrograde `c` occults `b` for 180 minutes.
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = False, oversample = 10, airless = True, seed = 1234)

  # Get the next occultation
  t = system.next_occultation(2000, system.b, occultor = system.c)
  system.b.limbdark = [0]
  time = np.arange(t - 0.15, t + 0.15, 5. / 1440.)

  # Compute and plot the light curve
  system.compute(time)
  system.plot_lightcurve(15.)

  # Observe it (one exposure)
  system.observe(stack = 1)
  pl.show()

# Instantiate the Trappist-1 system
system = Trappist1(sample = False, oversample = 10, airless = True, seed = 1)
system.b.Omega = 0
system.c.Omega = 0
system.b.inc = 90
system.c.inc = 90
time = np.arange(10071.6, 10071.85, 10. / 1440.)

# Compute and plot the light curve
system.compute(time)
system.plot_lightcurve(15.)

# Observe it (one exposure)
system.observe(stack = 1)
pl.show()