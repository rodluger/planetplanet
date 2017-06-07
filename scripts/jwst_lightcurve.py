'''
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
from planetplanet.detect import jwst
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
np.random.seed(1234)

cadence = 5.0 # mins
d = 12.2      # pc
pc_to_meter = u.pc.in_units(u.m)

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, ttvs = False, phasecurve = True, adaptive = True)

# Get the occultation light curves for the first 10 days
time = np.linspace(0., 10., 10000)
system.compute(time, lambda1 = 5, lambda2 = 30, R = 2000)

# Calculate flux at a distance, d
flux = system.flux / (d * pc_to_meter)**2
lam = system.wavelength

# Get MIRI Filter "wheel"
wheel = jwst.get_miri_filter_wheel()

# Compute MIRI lightcurves
jwst.lightcurves(wheel, flux, time, lam, obscad=cadence, plot=True)
