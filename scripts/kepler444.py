#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler444.py
------------

Estimates the signal-to-noise on secondary eclipses for the innermost planet (b)
in the Kepler-444 system with JWST
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
from planetplanet.detect import jwst
import matplotlib.pyplot as pl
import numpy as np

import astropy.units as u

# Kepler-444 and simulation parameters
Nocc = 15         # number of occultations observed
nout = 4.0       # obseved out of transit durations [tint]
lammin = 1.0     # min wavelength [um]
lammax = 30.0    # max wavelength [um]
Tstar = 5040.    # stellar temperature [K]
Rs = 0.752       # stellar radius [Solar Radii]
Rp = 0.4         # planet radius [Earth Radii]
d = 36.          # distance to system [pc]

# Additional params for K-444b
r = 0.04178      # semi-major axis [AU]
A = 0.25         # Planet albedo
e = 1.0          # Planet emissivity
i = 88.0         # Orbital inclination [degrees]
P = 3.6          # Orbital period [days]

# Convert some units to km
Rs_km = Rs * u.solRad.in_units(u.km)
Rp_km = Rp * u.earthRad.in_units(u.km)
r_km = r * u.AU.in_units(u.km)
P_mins = P * 24. * 60.

# Planet temperature [K]
Tplan = Tstar * ((1.-A)/e)**0.25 * (0.5*Rs_km/r_km)**0.5

# transit duration
tdur_mins = (P_mins / np.pi) * np.arcsin(np.sqrt((Rp_km+Rs_km)**2 - r_km/(Rs_km*np.cos(i)))/r_km)

# integration time [seconds]
tdur = Nocc * tdur_mins * 60.0

print("Planet Temperature : %.1f K" %Tplan)
print("Transit Duration : %.2f mins" %(tdur_mins))

#jwst.estimate_eclipse_snr(tint = tdur, nout = nout, lammin = lammin, lammax = lammax,
#                          Tstar = Tstar, Tplan = Tplan, Rs = Rs, Rp = Rp, d = d)

# Estimate for OST
jwst.estimate_eclipse_snr(tint = tdur, nout = nout, lammin = lammin, lammax = lammax,
                          Tstar = Tstar, Tplan = Tplan, Rs = Rs, Rp = Rp, d = d,
                          atel = 144., thermal = False)
