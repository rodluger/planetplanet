#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
timestep.py
-----------

Testing the minimum timestep we need in the N-body code.
We integrate for one year and look at the errors as a fraction
of the planet radius. Looks like a timestep of 1 hour leads to
negligible (< 1 percent) error over 1 year.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np

fig, ax = pl.subplots(3)
time = np.linspace(0, 365, 10000)

# Tiny timestep (1 minute)
np.random.seed(1234)
system1 = Trappist1(nbody = True, dt = 1. / 1440)
system1.compute(time)

# Large timestep (1 hour)
np.random.seed(1234)
system2 = Trappist1(nbody = True, dt = 1. / 24)
system2.compute(time)

for body1, body2 in zip(system1.bodies[1:], system2.bodies[1:]):
  ax[0].plot(system1.time, (body1.x - body2.x) / body1._r)
  ax[1].plot(system1.time, (body1.y - body2.y) / body1._r)
  ax[2].plot(system1.time, (body1.z - body2.z) / body1._r)

pl.show()