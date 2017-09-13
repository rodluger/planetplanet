#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
kepler444.py |github|
---------------------

Estimates the signal-to-noise on 15 stacked secondary eclipses for the innermost planet (b)
in the Kepler-444 system with JWST.

  .. plot::
     :align: center
     
     from scripts import kepler444
     kepler444._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/kepler444.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet.constants import *
from planetplanet import Star, Planet, System, jwst
import matplotlib.pyplot as pl
import numpy as np
import astropy.units as u

def _test():
    '''
    
    '''
    
    fig1, _, fig2, _ = compute()
    pl.close(fig2)
    pl.show()

def compute():
    '''
    
    '''
    
    # Kepler-444 and simulation parameters
    Nocc = 15                   # number of occultations observed
    nout = 4.0                  # obseved out of transit durations [tint]
    lammin = 1.0                # min wavelength [um]
    lammax = 30.0               # max wavelength [um]
    Tstar = 5040.               # stellar temperature [K]
    Rs = 0.752                  # stellar radius [Solar Radii]
    Rp = 0.4                    # planet radius [Earth Radii]
    d = 36.                     # distance to system [pc]

    # Additional params for K-444b
    r = 0.04178                 # semi-major axis [AU]
    A = 0.25                    # Planet albedo
    e = 1.0                     # Planet emissivity
    i = 88.0                    # Orbital inclination [degrees]
    P = 3.6                     # Orbital period [days]

    # Convert some units to km
    Rs_km = Rs * u.solRad.in_units(u.km)
    Rp_km = Rp * u.earthRad.in_units(u.km)
    r_km = r * u.AU.in_units(u.km)
    P_mins = P * 24. * 60.

    # Planet temperature [K]
    Tplan = Tstar * ((1.-A)/e)**0.25 * (0.5*Rs_km/r_km)**0.5

    # transit duration
    tdur_mins = (P_mins / np.pi) * np.arcsin(np.sqrt((Rp_km + Rs_km) ** 2 
                                   - r_km / (Rs_km * np.cos(i))) / r_km)

    # integration time [seconds]
    tdur = Nocc * tdur_mins * 60.0

    print("Planet Temperature : %.1f K" %Tplan)
    print("Transit Duration : %.2f mins" %(tdur_mins))

    # Estimate for JWST
    fig1, ax1 = jwst.estimate_eclipse_snr(tint = tdur, nout = nout, 
                                          lammin = lammin, lammax = lammax,
                                          Tstar = Tstar, Tplan = Tplan, 
                                          Rs = Rs, Rp = Rp, d = d)

    # Estimate for an OST-like telescope, but using the JWST filters
    fig2, ax2 = jwst.estimate_eclipse_snr(tint = tdur, nout = nout, 
                                          lammin = lammin, lammax = lammax,
                                          Tstar = Tstar, Tplan = Tplan, 
                                          Rs = Rs, Rp = Rp, d = d,
                                          atel = 144., thermal = False)

    return fig1, ax1, fig2, ax2

if __name__ == '__main__':
    fig1, _, fig2, _ = compute()
    fig1.savefig("kepler444_jwst.pdf", bbox_inches = 'tight')
    fig2.savefig("kepler444_ost.pdf", bbox_inches = 'tight')