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

__all__ = ['LimbDarkenedMap', 'RadiativeEquilibriumMap', 'BandedCloudsMap', 'UniformMap']

def LimbDarkenedMap():
  '''
  
  '''
  
  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''
    
    return np.nan
  
  # Specify that this is the default radially-symmetric function
  func.maptype = MAP_RADIAL_DEFAULT
  
  return func

def RadiativeEquilibriumMap():
  '''
  
  '''
  
  @cfunc("float64(float64, float64)")
  def func(lam, z):
    '''
    
    '''
    
    return np.nan
  
  # Specify that this is the default elliptically-symmetric function
  func.maptype = MAP_ELLIPTICAL_DEFAULT
  
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
  
  # Specify that this is a custom radially-symmetric function
  func.maptype = MAP_RADIAL_CUSTOM
  
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
  
  # Specify that this is a custom elliptically-symmetric function
  func.maptype = MAP_ELLIPTICAL_CUSTOM
  
  return func