'''
test_jwst_lightcurves.py
------------------------
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from planetplanet.detect import jwst
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(1234)

cadence = 5.0 # mins
tfact = 0.001 # amplitude of variations
tmin = 0.0    # days
tmax = 4.0    # days
Nthr = 10000
Nlhr = 1000

# Get MIRI Filter "wheel"
wheel = jwst.get_miri_filter_wheel()

# Create fake time-dependent spectra
t, l, flux = jwst.create_fake_data(Nlhr, Nthr, tmin=tmin, tmax=tmax, tfact=tfact)
time = t[0]
lam = l[0]

# Compute MIRI lightcurves
jwst.lightcurves(wheel, flux, time, lam, obscad=cadence)
