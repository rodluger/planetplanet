#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
custom_map.py |github|
----------------------

Illustrates how to specify a custom surface radiance map and 
use it to visualize the planet and compute occultation light curves.
As an example, I specify a planet with a hotspot given by a Gaussian
in zenith angle, but this can be modified to generate any surface map
that is symmetric about the hotspot.

  .. plot::
     :align: center
     
     from scripts import custom_map
     custom_map._test()
     
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/custom_map.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import
import matplotlib.pyplot as pl
from planetplanet import DrawEyeball, Planet, Star, System
from planetplanet.constants import *
from numba import cfunc
import numpy as np

def _test():
  '''
  
  '''
  
  view_planet()
  secondary_eclipse()
  pl.show()

def CustomMap(fwhm = 30, temp = 500):
  '''
  This is a function generator that returns a surface map. 
  It can have an arbitrary number of args/kwargs, but for
  simplicity we allow the user to tune two parameters: the FWHM of the
  hotspot in degrees, which for this map we'll assume is a Gaussian in 
  zenith angle, and the temperature of the hotspot in K.

  '''
  
  # Note that you can do all sorts of calculations within this function to
  # define variables that can be accessed by the radiance map below. Let's
  # compute the standard deviation (in radians) of the Gaussian from the FWHM:
  std = (np.pi / 180) * fwhm / 2.35482

  # Note that since :py:obj:`func` is a compiled C function, you'll get nasty
  # errors if you try to access functions/classes/dictionaries or other fancy
  # Python objects from inside it. Stick to floats and ints and (maybe) numpy
  # arrays.

  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    The actual function that returns the radiance at a given wavelength and zenith angle.
    
    '''
    
    # Let's compute the radiance of the hotspot from Planck's law.
    # This will be the amplitude of our Gaussian.
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
    amp = a / (np.exp(b) - 1.)

    # Evaluate the Gaussian and return the radiance at the requested
    # zenith angle.
    return amp * np.exp(-(z / (2 * std)) ** 2)

  # This line is important: we need to tell :py:obj:`planetplanet` what kind
  # of radiance map this is. There are two choices: a radially symmetric map,
  # whose axis of symmetry is always the center of the projected disk, and an
  # elliptically symmetric (eyeball) map, whose axis of symmetry rotates as the planet
  # orbits the star, always pointing towards the substellar point (or with
  # an offset, if :py:obj:`Lambda` and :py:obj:`Phi` are set). Let's make this
  # map elliptically symmetric:

  func.maptype = MAP_ELLIPTICAL_CUSTOM

  # We return the actual radiance map when this function is called
  return func

def view_planet(theta = 0.87, gamma = 4.02):
  '''
  View an image of the planet w/ the custom map at a given phase angle.
  Note that I tuned `theta` and `gamma` to match the hotspot offset specified in
  the `secondary_eclipse` function below.
  
  '''
  
  fig = pl.figure(figsize = (3,3))
  DrawEyeball(fig = fig, radiancemap = CustomMap(), theta = theta, gamma = gamma, nz = 51)

def secondary_eclipse():
  '''
  Compute and plot the planet's secondary eclipse light curve.
  
  '''

  # Instantiate the star
  star = Star('A', r = 0.1)

  # Planet b. Let's make it quite big
  b = Planet('b', m = 1, per = 3, inc = 89.8, r = 5., t0 = -1.5, 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = True,
             radiancemap = CustomMap())
  
  # Give the planet a hotspot offset just for fun
  b.Phi = 30
  b.Lambda = 30
  
  # Compute and plot the light curve
  time = np.linspace(-0.1, 0.1, 1000)
  system = System(star, b)
  system.compute(time)
  system.plot_occultation('b', 0.)
  
if __name__ == '__main__':
  view_planet()
  secondary_eclipse()
  pl.show()