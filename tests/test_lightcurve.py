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

def test_mutual():
  '''
  
  '''
  
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
  
  # Instantiate the star
  star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [u1])

  # Planet b
  b = Planet('b', m = 1, per = 3, inc = 89.8, r = 3., trn0 = 0, 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

  # Planet c
  c = Planet('c', m = 1, per = 30, inc = 90., r = 3., trn0 = 0., 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

  # System
  system = System(star, b, c)

  # Get the occultation light curves
  time = np.linspace(-0.06, 0.06, 1000)
  system.compute(time)
  
  # Check some benchmarked values
  assert np.abs(system.flux[500,0] - 5.11448257092e-15) / 5.11448257092e-15 < TOL, "Incorrect flux."
  assert np.abs(system.flux[500,-1] - 9.55103324782e-17) / 9.55103324782e-17 < TOL, "Incorrect flux."
  assert np.abs(np.sum(system.flux[:,0]) - 5.7110429019e-12) / 5.7110429019e-12 < TOL, "Incorrect average flux."
  assert np.abs(np.sum(system.flux[:,-1]) - 1.03675617599e-13) / 1.03675617599e-13 < TOL, "Incorrect average flux."