#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_lightcurve.py
------------------

Test the occultation code.

'''

from planetplanet import Star, Planet, System
import numpy as np
TOL = 0.001
RSUN = 6.957e8
PARSEC = 3.086e16
MICRON = 1.e-6
SBOLTZ = 5.670367e-8

def u1(lam):
  '''
  A really silly linear limb darkening law with a linear
  wavelength dependence.

  '''

  lam = np.atleast_1d(lam)
  result = 0.5 * (1 - (lam - 5) / 10) + 0.5
  result[lam < 5] = 1.
  result[lam > 15] = 0.5
  return result

def test_mutual():
  '''
  
  '''

  # Instantiate the star
  star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [u1])

  # Planet b
  b = Planet('b', m = 1, per = 3, inc = 89.8, r = 3., t0 = 0., 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

  # Planet c
  c = Planet('c', m = 1, per = 30, inc = 90., r = 3., t0 = 0., 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

  # System
  system = System(star, b, c)

  # Get the occultation light curves
  time = np.linspace(-0.06, 0.06, 1000)
  system.compute(time)
  
  # Benchmarked values
  truths = [7.6963687220060673e-15,
            1.1495664454458829e-16,
            8.5930679511725155e-12,
            1.247748695405246e-13]
  
  # Computed values
  values = [system.flux[500,0],
            system.flux[500,-1],
            np.sum(system.flux[:,0]),
            np.sum(system.flux[:,-1])]

  # Check!
  assert np.abs(values[0] - truths[0]) / truths[0] < TOL, "Incorrect flux."
  assert np.abs(values[1] - truths[1]) / truths[1] < TOL, "Incorrect flux."
  assert np.abs(values[2] - truths[2]) / truths[2] < TOL, "Incorrect average flux."
  assert np.abs(values[3] - truths[3]) / truths[3] < TOL, "Incorrect average flux."

def test_limbdark():
  '''
  
  '''
  
  # Instantiate the star
  r = 0.1
  teff = 3200
  star = Star('A', m = 0.1, r = r, nz = 31, color = 'k', limbdark = [u1], teff = teff)

  # True luminosity
  truth = 4 * np.pi * (r * RSUN) ** 2 * SBOLTZ * teff ** 4

  # System
  system = System(star)
  system.distance = 12.2
  
  # Compute the light curve
  time = np.arange(0., 1., 10)
  system.compute(time, lambda1 = 0.1, lambda2 = 100, R = 3000)
  
  # Luminosity
  bol = np.trapz(system.flux[0], system.A.wavelength * 1e6)
  lum = 4 * np.pi * bol * (12.2 * PARSEC) ** 2
  
  # Check!
  assert np.abs(lum - truth) / truth < TOL, "Incorrect bolometric flux."

if __name__ == '__main__':
  test_mutual()
  test_limbdark()