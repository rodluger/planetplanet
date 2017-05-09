#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
simple.py
---------

Simple planet-planet occultation flux for a constant
brightness dayside and a constant brightness nightside.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
np.seterr(divide = 'ignore', invalid = 'ignore')
import matplotlib.pyplot as pl

__all__ = ['dFlux']

def EllipseUpper(x, r, xE, yE, a, b):
  '''
  The function describing the upper half of the half-ellipse
  
  '''
  
  A = b ** 2 - (x - xE) ** 2
  if np.abs(A) < 1e-15:
    A = 0
  return yE + a / b * np.sqrt(A)

def EllipseLower(x, r, xE, yE, a, b):
  '''
  The function describing the lower half of the half-ellipse
  
  '''
  
  return yE

def CircleUpper(x, r, xE, yE, a, b):
  '''
  The function describing the upper half of the circle
  
  '''
  
  A = r ** 2 - x ** 2
  if np.abs(A) < 1e-15:
    A = 0
  return np.sqrt(A)

def CircleLower(x, r, xE, yE, a, b):
  '''
  The function describing the lower half of the circle
  
  '''
  
  return -CircleUpper(x, r, xE, yE, a, b)

def IntCircleUpper(x0, x1, r, xE, yE, a, b):
  '''
  The definite integral of the area under the upper
  curve of the circle between `x0` and `x1`
  
  '''
  
  F = [0, 0]
  for i, x in enumerate((x0, x1)):
    z = (r - x) * (r + x)
    if np.abs(b) < 1e-15:
      z = 0
    z = np.sqrt(z)
    F[i] = 0.5 * (x * z + r ** 2 * np.arctan(x / z))
  return F[1] - F[0]

def IntCircleLower(x0, x1, r, xE, yE, a, b):
  '''
  The definite integral of the area under the lower
  curve of the circle between `x0` and `x1`
  
  '''
  
  return -IntCircleUpper(x0, x1, r, xE, yE, a, b)
  
def IntEllipseUpper(x0, x1, r, xE, yE, a, b):
  '''
  The definite integral of the area under the upper
  curve of the half-ellipse between `x0` and `x1`
  
  '''
  
  F = [0, 0]
  for i, x in enumerate((x0, x1)):
    y = x - xE
    z = (b - y) * (b + y)
    if np.abs(b) < 1e-15:
      z = 0
    z = np.sqrt(z)
    F[i] = (a / (2 * b)) * (y * z + b ** 2 * np.arctan(y / z)) + yE * x
  return F[1] - F[0]

def IntEllipseLower(x0, x1, r, xE, yE, a, b):
  '''
  The definite integral of the area under the lower
  curve of the half-ellipse between `x0` and `x1`
  
  '''
  
  return yE * (x1 - x0)
  
def Roots(r, xE, yE, a, b):
  '''
  Returns the points of intersection of the top section of a half-ellipse
  centered at (`xE`, `yE`) with semi-major axis `b` and semi-minor
  axes `a` (bad notation, sorry!) and a circle of radius `r` centered at 
  the origin
  
  '''
  
  # Get the coefficients
  A = a ** 2 / b ** 2 - 1
  B = -2 * a ** 2 * xE / b ** 2
  C = r ** 2 - yE ** 2 - a ** 2 + a ** 2 / b ** 2 * xE ** 2
  D = 4 * yE ** 2 * a ** 2 / b ** 2
  c4 = A ** 2
  c3 = 2 * A * B
  c2 = 2 * A * C + B ** 2 + D
  c1 = 2 * B * C - 2 * D * xE
  c0 = C ** 2 - (b ** 2 - xE ** 2) * D
  
  # Solve numerically
  roots = np.roots([c4, c3, c2, c1, c0])
  
  # Filter out imaginary roots and roots corresponding
  # to the lower half of the ellipse
  good_roots = []
  for x in roots:
    if (x.imag == 0):
      x = x.real
      args = (x, r, xE, yE, a, b)
      if (np.abs(EllipseUpper(*args) - CircleUpper(*args)) < 1e-7) or \
         (np.abs(EllipseUpper(*args) - CircleLower(*args)) < 1e-7):
        good_roots.append(x)
      
  return good_roots

def Limits(r, xE, yE, a, b):
  '''
  Returns tuples of the limits of integration for
  the integrals describing the area of overlap
  between the half-ellipse and the circle
  
  '''
  
  # The integration limits
  limits = []
  
  # Check for ellipse vertices inside circle, including edge cases
  for x, y in [(xE - b, yE), (xE + b, yE)]:
    if (x ** 2) + (y ** 2) <= r ** 2:
      limits.append(x)
  
  # Check for circle vertices inside ellipse, including edge cases
  for x, y in [(-r, 0), (r, 0)]:
    if (y >= yE) and (x - xE) ** 2 / b ** 2 + (y - yE) ** 2 / a ** 2 <= 1:
      limits.append(x)
  
  # Check for circle-line intersections, excluding edge cases
  if r ** 2 >= yE ** 2:
    x = np.sqrt(r ** 2 - yE ** 2)
    if (x > xE - b) and (x < xE + b):
      limits.append(x)
    if (-x > xE - b) and (-x < xE + b):
      limits.append(-x)
  
  # Now get the circle-ellipse intersections, excluding edge cases
  for x in Roots(r, xE, yE, a, b):
    limits.append(x)
  
  # Sort the limits and transform to tuples
  limits = sorted(limits)
  limits = list(zip(limits[:-1], limits[1:]))
  
  return limits

def Area(r, xE, yE, a, b):
  '''
  The area of overlap between the circle and the
  half-ellipse
  
  '''
  
  # Get the integration limits
  limits = Limits(r, xE, yE, a, b)
  area = np.zeros(len(limits))
  
  # Loop over all the sub-regions
  for i, _ in enumerate(area):
    
    # The integration limits
    x0, x1 = limits[i]
    
    # Bisect the limits and sort the boundary functions
    x = (x0 + x1) / 2
    args = (x, r, xE, yE, a, b)
    order = np.argsort(np.array([EllipseUpper(*args), EllipseLower(*args), CircleUpper(*args), CircleLower(*args)]))
    
    # Grab the interior two boundary functions and integrate them!
    funcs = [IntEllipseUpper, IntEllipseLower, IntCircleUpper, IntCircleLower]
    intupper = funcs[order[2]](x0, x1, r, xE, yE, a, b)
    intlower = funcs[order[1]](x0, x1, r, xE, yE, a, b)
    
    # The net area for this region is just their difference
    area[i] = intupper - intlower
    
  return np.sum(area)

def dFlux(xp, yp, rp, xo, yo, ro, theta = np.pi / 8, fday = 1, fnight = 0.5):
  '''
  Returns the flux difference when the planet (`p`) is occulted by the occultor (`o`).
  The planet is centered at (`xp`, `yp`) and has radius `rp`, and the occultor is 
  centered at (`xo`, `yo`) and has radius `ro`. The star is some distance `d` to the 
  left at (`xp` - `d`, `yp`), and the phase angle of the planet is `theta`, such that
  -pi/2, 0, and pi/2 correspond to new phase, half phase, and full phase, respectively.
  
  '''
    
  # Is the terminator ellipse dark or bright?
  assert (theta >= -np.pi / 2) and (theta <= np.pi / 2), "Invalid value for theta."
  tdark = theta < 0
  theta = np.abs(theta)
  
  # The position of the planet when the occultor is at the origin
  x = xp - xo
  y = yp - yo
  d = np.sqrt(x ** 2 + y ** 2)
  
  # Are the circles touching at all?
  if d >= (ro + rp):
    return 0
  
  # The areas of the day and night half-disks
  DHALF = Area(ro, y, -x, rp, rp)
  NHALF = Area(ro, -y, x, rp, rp)
  
  # The area of the terminator half-ellipse
  # and the day and night regions
  if tdark:
    fterm = fnight
    T = Area(ro, y, -x, rp * np.sin(theta), rp)
    D = DHALF - T
    N = NHALF
  else:
    fterm = fday
    T = Area(ro, -y, x, rp * np.sin(theta), rp)
    N = NHALF - T
    D = DHALF

  # Combine everything for the net flux change
  dF = -(fday * D + fterm * T + fnight * N)

  return dF