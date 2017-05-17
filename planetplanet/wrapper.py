#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`wrapper.py` - C wrapper
--------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
np.seterr(invalid = 'ignore')
import os
from numpy.ctypeslib import ndpointer, as_ctypes
import matplotlib.pyplot as pl
cmap = pl.get_cmap('RdBu_r')
AUREARTH = 23454.9271
MDFAST = 0
NEWTON = 1

# Load the C library
try:
  ppo = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")

class Planet(ctypes.Structure):
  '''
  The class containing all the input planet parameters
  
  '''
  
  _fields_ = [("per", ctypes.c_double),
              ("inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("w", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("r", ctypes.c_double),
              ("noon", ctypes.c_double),
              ("midnight", ctypes.c_double),
              ("nlat", ctypes.c_int),
              ("x", ctypes.c_double),
              ("y", ctypes.c_double),
              ("z", ctypes.c_double)]
        
  def __init__(self, **kwargs):
    self.per = kwargs.pop('per', 1.51087081)
    inc = kwargs.pop('inc', 89.65)
    self.inc = inc * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.)
    self.a = kwargs.pop('a', 0.01111 * AUREARTH)
    self.t0 = kwargs.pop('t0', 7322.51736)
    self.r = kwargs.pop('r', 1.086)
    self.noon = kwargs.pop('noon', 1.)
    self.midnight = kwargs.pop('midnight', 0.1)
    self.nlat = kwargs.pop('nlat', 11)
    self.x = np.nan
    self.y = np.nan
    self.z = np.nan
    
class Settings(ctypes.Structure):
  '''
  The class that contains the model settings
  
  '''
  
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int),
              ("phasecurve", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', NEWTON)
    self.phasecurve = int(kwargs.pop('phasecurve', False))
    
def Image(ax, occultor, occulted, pad = 0.5):
  '''
  
  '''
  
  # Plot the occultor
  r = occultor.r
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', zorder = 99)
  ax.plot(x, -y, color = 'k', zorder = 99)
  ax.fill_between(x, -y, y, color = 'lightgray', zorder = 99)

  # Plot the occulted planet
  r = occulted.r
  x0 = occulted.x - occultor.x
  y0 = occulted.y - occultor.y
  theta = np.arctan(occulted.z / np.abs(occulted.x))
  x = np.linspace(x0 - r, x0 + r, 1000)
  y = np.sqrt(r ** 2 - (x - x0) ** 2)
  ax.plot(x, y0 + y, color = 'k', zorder = 98)
  ax.plot(x, y0 - y, color = 'k', zorder = 98)
  
  # Plot the latitude ellipses
  for lat in np.linspace(0, np.pi, occulted.nlat + 2)[1:-1]:
    a = occulted.r * np.abs(np.sin(lat))
    b = a * np.abs(np.sin(theta))
    x0 = occulted.x - occulted.r * np.cos(lat) * np.cos(theta)
    y0 = occulted.y
    xlimb = occulted.r * np.cos(lat) * np.sin(theta) * np.tan(theta)
    if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
      xmin = x0 - b
    else:
      xmin = x0 - xlimb
    if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
      xmax = x0 + b
    else:
      xmax = x0 - xlimb
    x = np.linspace(x0 - b, x0 + b, 1000)
    if theta > 0:
      x[x < x0 - xlimb] = np.nan
    else:
      x[x > x0 - xlimb] = np.nan
    A = b ** 2 - (x - x0) ** 2
    A[A<0] = 0
    y = (a / b) * np.sqrt(A)
    if np.abs(np.cos(lat)) < 1e-5:
      style = dict(color = 'k', ls = '--')
    else:
      style = dict(color = cmap(0.5 * (np.cos(lat) + 1)), ls = '-')
    ax.plot(x - occultor.x, occulted.y - occultor.y + y, **style)
    ax.plot(x - occultor.x, occulted.y - occultor.y - y, **style)

  # Limits
  ax.set_xlim(occulted.x - occultor.x - (pad + 1) * r, occulted.x - occultor.x + (pad + 1) * r)
  ax.set_ylim(occulted.y - occultor.y - (pad + 1) * r, occulted.y - occultor.y + (pad + 1) * r)

def Lightcurve(time, planets, **kwargs):
  '''
  
  '''
  
  # Initialize
  Flux = ppo.Flux
  Flux.restype = ctypes.c_double
  Flux.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.POINTER(Planet), 
                   Settings, ctypes.ARRAY(ctypes.c_double, len(planets))]
  if type(planets) is Planet:
    planets = [planet]
  n = len(planets)
  p = (Planet * n)(*planets)
  f = ctypes.ARRAY(ctypes.c_double, n)()
  s = Settings(**kwargs)
  
  # Call the C function in a loop
  flux = np.zeros((len(time), n))
  for i, t in enumerate(time):
    Flux(t, n, p, s, f)
    flux[i] = np.array(f)

  return flux


b = Planet(per = 1.51087081, inc = 89.65, a = 0.01111 * AUREARTH, r = 1.086, t0 = 7322.51736, midnight = 1, n = 31)
c = Planet(per = 2.4218233, inc = 89.67, a = 0.01521 * AUREARTH, r = 1.056, t0 = 7282.80728, midnight = 1, n = 31)
d = Planet(per = 4.049610, inc = 89.75, a = 0.02144 * AUREARTH, r = 0.772, t0 = 7670.14165, midnight = 1, n = 31)
e = Planet(per = 6.099615, inc = 89.86, a = 0.02817 * AUREARTH, r = 0.918, t0 = 7660.37859, midnight = 1, n = 31)
f = Planet(per = 9.206690, inc = 89.68, a = 0.0371 * AUREARTH, r = 1.045, t0 = 7671.39767, midnight = 1, n = 31)
g = Planet(per = 12.35294, inc = 89.71, a = 0.0451 * AUREARTH, r = 1.127, t0 = 7665.34937, midnight = 1, n = 31)
h = Planet(per = 18.766, inc = 89.80, a = 0.06 * AUREARTH, r = 0.755, t0 = 7662.55463, midnight = 1, n = 31)
planets = [b, c, d, e, f, g, h]

time = np.linspace(0, 100, 10000)
flux = Lightcurve(time, planets, phasecurve = False)
fig, ax = pl.subplots(7)
for i in range(7):
  ax[i].plot(time, flux[:,i], lw = 1)
  #ax[i].set_xticklabels([])
  #ax[i].set_yticklabels([])
pl.show()