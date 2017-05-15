#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`orbit.py` - Orbital solution
-------------------------------------

A :py:mod:`ctypes` wrapper around a C implementation of a Kepler solver. 

.. warning:: The longitude of pericenter `w` may be defined differently here than in other transit codes; \
             watch out for a possible offset of :math:`\pi` from what you're used to.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
import os
from numpy.ctypeslib import ndpointer, as_ctypes

__all__ = ['Orbit']

# Define errors
_ERR_NONE             =   0                                                           # We're good!
_ERR_NOT_IMPLEMENTED  =   1                                                           # Function/option not yet implemented
_ERR_KEPLER           =   3                                                           # Error in the Kepler solver; probably didn't converge
_ERR_BAD_ECC          =   5                                                           # Bad value for eccentricity
_ERR_RADIUS           =   9                                                           # Bad input radius
_ERR_INC              =   10                                                          # Bad inclination
_ERR_STAR_CROSS       =   12                                                          # Star-crossing orbit
_ERR_PER              =   13                                                          # Bad period
_ERR_RHOS_ARS         =   14                                                          # Must specify either rhos or aRs!
_ERR_RHOS             =   15                                                          # Bad rhos
_ERR_ECC_W            =   16                                                          # Bad eccentricity/omega

# Define models
MDFAST     =              9
NEWTON     =              10

# Array IDs
_ARR_M       =             2
_ARR_E       =             3
_ARR_F       =             4
_ARR_R       =             5
_ARR_X       =             6
_ARR_Y       =             7
_ARR_Z       =             8
_ARR_B       =             9

# Other
G           =             6.672e-8
DAYSEC      =             86400.

class PLANET(ctypes.Structure):
      '''
      The class containing all the input orbital parameters
      
      '''
      
      _fields_ = [("inc", ctypes.c_double),
                  ("rhos", ctypes.c_double),
                  ("MpMs", ctypes.c_double),
                  ("esw", ctypes.c_double),
                  ("ecw", ctypes.c_double),
                  ("per", ctypes.c_double),
                  ("RpRs", ctypes.c_double),
                  ("ecc", ctypes.c_double),
                  ("w", ctypes.c_double),
                  ("aRs", ctypes.c_double)]
      
      def __init__(self, **kwargs):   
        self.MpMs = kwargs.pop('MpMs', 0.)
        self.per = kwargs.pop('per', 10.)
        self.RpRs = kwargs.pop('RpRs', 0.1)
        self.bcirc = kwargs.pop('bcirc', 0.)
        inc = kwargs.pop('inc', 90.)
        self.inc = inc * np.pi / 180.
        self.rhos = kwargs.pop('rhos', 1.4)                                           # User may specify either ``rhos`` or ``aRs``
        self.aRs = kwargs.pop('aRs', np.nan)
        if not np.isnan(self.aRs):
          self.rhos = np.nan
        self.ecc = kwargs.pop('ecc', 0.)                                              # User may specify (``esw`` and ``ecw``) or (``ecc`` and ``w``)
        self.w = kwargs.pop('w', 0.)
        self.esw = kwargs.pop('esw', np.nan)
        self.ecw = kwargs.pop('ecw', np.nan)
        if (not np.isnan(self.esw)) and (not np.isnan(self.ecw)):
          self.ecc = np.nan
          self.w = np.nan
         
class ARRAYS(ctypes.Structure):
      '''
      The class that stores the output arrays
      
      '''
      
      _fields_ = [("npts", ctypes.c_int),
                  ("_M", ctypes.POINTER(ctypes.c_double)),
                  ("_E", ctypes.POINTER(ctypes.c_double)),
                  ("_f", ctypes.POINTER(ctypes.c_double)),
                  ("_r", ctypes.POINTER(ctypes.c_double)),
                  ("_x", ctypes.POINTER(ctypes.c_double)),
                  ("_y", ctypes.POINTER(ctypes.c_double)),
                  ("_z", ctypes.POINTER(ctypes.c_double)),
                  ("_b", ctypes.POINTER(ctypes.c_double))]
                  
      def __init__(self, **kwargs):                
        pass
        
      @property
      def M(self):
        return np.array([self._M[i] for i in range(self.npts)])
        
      @property
      def E(self):
        return np.array([self._E[i] for i in range(self.npts)])
        
      @property
      def f(self):
        return np.array([self._f[i] for i in range(self.npts)])
        
      @property
      def r(self):
        return np.array([self._r[i] for i in range(self.npts)])
      
      @property
      def x(self):
        return np.array([self._x[i] for i in range(self.npts)])
        
      @property
      def y(self):
        return np.array([self._y[i] for i in range(self.npts)])
      
      @property
      def z(self):
        return np.array([self._z[i] for i in range(self.npts)])
      
      @property
      def b(self):
        return np.array([self._b[i] for i in range(self.npts)])
             
class SETTINGS(ctypes.Structure):
      '''
      The class that contains the light curve settings
      
      '''
      _fields_ = [("keptol", ctypes.c_double),
                  ("maxkepiter", ctypes.c_int),
                  ("kepsolver", ctypes.c_int)]
      
      def __init__(self, **kwargs):
        self.keptol = kwargs.pop('keptol', 1.e-15)                                    # Kepler solver tolerance
        self.maxkepiter = kwargs.pop('maxkepiter', 100)                               # Maximum number of iterations in Kepler solver
        self.kepsolver = kwargs.pop('kepsolver', NEWTON)                              # Newton solver or fast M&D solver?

# Load the C library
try:
  lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'orbitlib.so'))
except:
  raise Exception("Can't find `orbitlib.so`; please run `make` to compile it.")

# Declare the C functions; user should access these through the Orbit() class below
_Compute = lib.Compute
_Compute.restype = ctypes.c_int
_Compute.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(PLANET), ctypes.POINTER(SETTINGS), ctypes.POINTER(ARRAYS)]

_dbl_free = lib.dbl_free
_dbl_free.argtypes = [ctypes.POINTER(ctypes.c_double)]

# Error handling
def RaiseError(err):
  if (err == _ERR_NONE):
    return
  elif (err == _ERR_NOT_IMPLEMENTED):
    raise Exception("Option not implemented.")
  elif (err == _ERR_BAD_ECC):
    raise Exception("Bad value for ``ecc``.")  
  elif (err == _ERR_RADIUS):
    raise Exception("Bad value for ``RpRs``.")
  elif (err == _ERR_INC):
    raise Exception("Bad value for ``inc``.")
  elif (err == _ERR_STAR_CROSS):
    raise Exception("Star-crossing orbit.") 
  elif (err == _ERR_RHOS_ARS):
    raise Exception("Must specify one of ``rhos`` or ``aRs``.") 
  elif (err == _ERR_RHOS):
    raise Exception("Bad value for ``per``.") 
  elif (err == _ERR_ECC_W):
    raise Exception("Bad value for ``esw`` or ``ecw``.") 
  elif (err == _ERR_KEPLER):
    raise Exception("Error in Kepler solver.")
  else:
    raise Exception("Error in orbital computation (%d)." % err)

class Orbit(object):
  '''
  A user-friendly wrapper around the :py:class:`ctypes` routines.
  
  :param time: The input time array.
  
  :param kwargs: The keyword arguments. These are:
  
    - **inc** - The orbital inclination in degrees. Default `90.`
    - **MpMs** - The planet-star mass ratio. Default `0.`
    - **per** - The planet orbital period in days. Default `10.`
    - **RpRs** - The planet-star radius ratio. Default `0.1`
    - **rhos** or **aRs** - The stellar density in `g/cm^3` or the semi-major axis-stellar radius ratio. \
                            Default is `rhos = 1.4`, the density of the Sun
    - **ecc** and **w** or **esw** and **ecw** - The eccentricity and the longitude of pericenter in radians, \
                                                 or the two eccentricity vectors. Default is `ecc = 0.` and `w = 0.`  

    - **keptol** - The tolerance of the Kepler solver. Default `1.e-15`
    - **maxkepiter** - Maximum number of iterations in the Kepler solver. Default `100`
    - **kepsolver** - The Kepler solver to use. Default `ps.NEWTON` (recommended)
  
  '''
  
  def __init__(self, time, **kwargs):
    self.time = time
    self.planet = PLANET(**kwargs)
    self.settings = SETTINGS(**kwargs)
    self.arrays = ARRAYS() 
    err = _Compute(len(self.time), as_ctypes(self.time), self.planet, self.settings, self.arrays)
    if err != _ERR_NONE: 
      RaiseError(err)
  
  @property
  def M(self):
    return self.arrays.M

  @property
  def E(self):
    return self.arrays.E
  
  @property
  def f(self):
    return self.arrays.f

  @property
  def r(self):
    return self.arrays.r

  @property
  def x(self):
    return self.arrays.x

  @property
  def y(self):
    return self.arrays.y

  @property
  def z(self):
    return self.arrays.z

  @property
  def b(self):
    return self.arrays.b
 
  def __del__(self):
    '''
    Free the C arrays when the last reference to the class goes out of scope.
    
    '''
    
    _dbl_free(self.arrays._M)
    _dbl_free(self.arrays._E)
    _dbl_free(self.arrays._f)
    _dbl_free(self.arrays._r)
    _dbl_free(self.arrays._x)
    _dbl_free(self.arrays._y)
    _dbl_free(self.arrays._z)