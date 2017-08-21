#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
maps.py |github|
----------------

A library of surface radiance map functions.

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/maps.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import
from ..constants import *
import numpy as np
from numba import cfunc

__all__ = ['RadiativeEquilibriumMap', 'LimbDarkenedMap', 'BandedCloudsMap', 'UniformMap', 'NullMap']

def Blackbody(lam, temp):
  '''
  
  '''
  
  a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
  b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
  return a / (np.exp(b) - 1.)

def NullMap():
  '''
  
  '''
  
  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''
    
    return np.nan
  
  return func

def UniformMap(albedo = 0.3, irrad = SEARTH):
  '''
  
  '''
  
  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''
    
    # Compute the temperature
    temp = ((irrad * (1 - albedo)) / (4 * SBOLTZ)) ** 0.25
  
    # Compute the radiance
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
    return a / (np.exp(b) - 1.)
  
  return func

def LimbDarkenedMap(teff = 2500, u = [1.], lambda1 = 5., lambda2 = 15., R = 100):
  '''
  
  '''
  
  # Compute the wavelength grid
  wav = [lambda1]
  while(wav[-1] < lambda2):
    wav.append(wav[-1] + wav[-1] / R)
  wavelength = np.array(wav)
  
  # Ensure limb darkening coefficients are a 2d array in wavelength and mu
  ulam = [None for ui in u]
  for i in range(len(u)):
    if callable(u[i]):
      ulam[i] = u[i](wavelength)
    elif not hasattr(u[i], '__len__'):
      ulam[i] = u[i] * np.ones_like(wavelength)
    else:
      raise Exception("Limb darkening coefficients must be provided as a list of scalars or as a list of functions.")
  u = np.array(ulam)
  nu, nw = u.shape
  B0 = np.zeros_like(wavelength)
  
  # Compute the normalization term, Equation (E5)
  for j, lam in enumerate(wavelength):
    norm = 0
    for i in range(nu):
      norm += u[i, j] / ((i + 2) * (i + 3))
    norm = 1 - 2 * norm
    B0[j] = Blackbody(lam * 1e-6, teff) / norm;
  
  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''

    # Get the index in the wavelength array 
    # (TODO: speed this up using `bsearch`?)
    for j in range(nw):
      if wavelength[j] > lam - 1e-10:
        break

    # Initialize
    flux = B0[j]
    cosz = np.cos(z)
  
    # Loop over the coefficient order
    for i in range(nu):
  
      # The Taylor expansion is in (1 - mu)
      x = (1 - cosz) ** (i + 1)
    
      # Compute the wavelength-dependent intensity
      flux -= u[i, j] * B0[j] * x
    
    return flux
  
  return func
  
def RadiativeEquilibriumMap(albedo = 0.3, tnight = 40., irrad = SEARTH):
  '''
  
  '''

  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''
    
    # Compute the temperature
    if (z < np.pi / 2):
      temp = ((irrad * np.cos(z) * (1 - albedo)) / SBOLTZ) ** 0.25
      if (temp < tnight):
        temp = tnight
    else:
      temp = tnight
  
    # Compute the radiance
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
    return a / (np.exp(b) - 1.)
  
  return func

def BandedCloudsMap(albedo = 0.3, irrad = SEARTH, bands = 5):
  '''
  
  '''

  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
  
    '''
    
    # Compute the temperature
    temp = ((irrad * np.cos(bands * z) ** 2 * (1 - albedo)) / SBOLTZ) ** 0.25
  
    # Compute the radiance
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
    return a / (np.exp(b) - 1.)
  
  return func