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
  lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")
Flux = lib.Flux
Flux.restype = ctypes.c_double
Flux.argtypes = [ctypes.c_double, ctypes.POINTER(PLANET), ctypes.POINTER(PLANET), ctypes.POINTER(SETTINGS)]

class PLANET(ctypes.Structure):
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
    
class SETTINGS(ctypes.Structure):
  '''
  The class that contains the model settings
  
  '''
  
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', NEWTON)

def Plot(ax, occultor, occulted, pad = 0.5):
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




settings = SETTINGS()
t1b = PLANET(per = 1.51087081, inc = 89.65, a = 0.01111 * AUREARTH, r = 1.086, t0 = 7322.51736, midnight = 0, n = 31)
t1c = PLANET(per = 2.4218233, inc = 89.67, a = 0.01521 * AUREARTH, r = 1.056, t0 = 7282.80728, midnight = 0, n = 31)
time = np.linspace(0, 50, 100000)
flux = [Flux(t, t1b, t1c, settings) for t in time]