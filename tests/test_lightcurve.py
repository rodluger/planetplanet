#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test_lightcurve.py
------------------

Test the occultation code.

'''

from planetplanet import Star, Planet, System
from planetplanet.constants import *
import numpy as np

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

def test_mutual(tol = 1e-10):
  '''
  Test the code's ability to compute mutual transit light curves.
  
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
  truths = [7.69636872134e-15,
            1.14956644064e-16,
            8.59304867308e-12,
            1.24774758808e-13]
  
  # Computed values
  values = [system.flux[500,0],
            system.flux[500,-1],
            np.sum(system.flux[:,0]),
            np.sum(system.flux[:,-1])]

  # Check!
  assert np.abs(values[0] - truths[0]) / truths[0] < tol, "Incorrect flux: %.10e != %.10e." % (values[0], truths[0])
  assert np.abs(values[1] - truths[1]) / truths[1] < tol, "Incorrect flux: %.10e != %.10e." % (values[1], truths[1])
  assert np.abs(values[2] - truths[2]) / truths[2] < tol, "Incorrect average flux: %.10e != %.10e." % (values[2], truths[2])
  assert np.abs(values[3] - truths[3]) / truths[3] < tol, "Incorrect average flux: %.10e != %.10e." % (values[3], truths[3])

def test_compare_methods(tol = 1e-6):
  '''
  Compare the different optimization methods for limb-darkened transits, including the `batman` algorithm.
  
  '''

  # Instantiate the star
  star = Star('star', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [0.4, 0.26])

  # Planet b
  planet = Planet('planet', m = 1, per = 0.5, inc = 90.4, r = 2., t0 = 0, 
                  nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

  # Time array
  time = np.arange(-0.025, 0.025, 0.1 * MINUTE)
  
  # Compute the light curve, no optimization
  system = System(star, planet, batmanopt = False, circleopt = False)
  system.compute(time, lambda1 = 0.5, lambda2 = 2.)
  flux1 = system.star.flux[:,0] / system.star.flux[0,0]
  
  # Compute the light curve w/ circle optimization
  system = System(star, planet, batmanopt = False, circleopt = True)
  system.compute(time, lambda1 = 0.5, lambda2 = 2.)
  flux2 = system.star.flux[:,0] / system.star.flux[0,0]
  
  # Compute the light curve w/ batman optimization
  system = System(star, planet, batmanopt = True, circleopt = False)
  system.compute(time, lambda1 = 0.5, lambda2 = 2.)
  flux3 = system.star.flux[:,0] / system.star.flux[0,0]
  
  # Flux minima
  flux1min = flux1.min()
  flux2min = flux2.min()
  flux3min = flux3.min()
  
  # Flux averages
  flux1mean = flux1.mean()
  flux2mean = flux2.mean()
  flux3mean = flux3.mean()
  
  # Check!
  assert (flux1min - flux2min) / flux1min < tol, "Incorrect flux: %.10e != %.10e" % (flux1min, flux2min)
  assert (flux1min - flux3min) / flux1min < tol, "Incorrect flux: %.10e != %.10e" % (flux1min, flux3min)
  assert (flux1mean - flux2mean) / flux1mean < tol, "Incorrect average flux: %.10e != %.10e" % (flux1mean, flux2mean)
  assert (flux1mean - flux3mean) / flux1mean < tol, "Incorrect average flux: %.10e != %.10e" % (flux1mean, flux3mean)
  
def test_limbdark(tol = 1e-4):
  '''
  Test the limb darkening normalization by comparing the total flux of the star
  to what you get with the Stefan-Boltzmann law.
  
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
  system.compute(time, lambda1 = 0.01, lambda2 = 1000, R = 3000)
  
  # Luminosity
  bol = np.trapz(system.flux[0], system.A.wavelength * 1e6)
  lum = 4 * np.pi * bol * (12.2 * PARSEC) ** 2
  
  # Check!
  assert np.abs(lum - truth) / truth < tol, "Incorrect bolometric flux: %.10e != %.10e." % (lum, truth)

if __name__ == '__main__':
  test_mutual()
  test_limbdark()
  test_compare_methods()