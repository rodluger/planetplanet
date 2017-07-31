#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
phases.py
---------

*Still working on this.*

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball, Star, Planet, System
import numpy as np
import matplotlib.pyplot as pl

# Get orbital elements
star = Star('A')
b = Planet('b', per = 10., inc = 60., Omega = 0., t0 = 0, ecc = 0, w = 0., dlambda = 0, dpsi = 0)
system = System(star, b)
time = np.linspace(-5, 5, 1000)
system.compute_orbits(time)

# Plot the orbit
fig, ax = pl.subplots(1, figsize = (8,8))
ax.plot(b.x, b.y, 'k-', alpha = 0.5)

# Adjust the plot range
xmin = min(b.x.min(), b.y.min())
xmax = max(b.x.max(), b.y.max())
dx = xmax - xmin
xmin -= 0.1 * dx
xmax += 0.1 * dx
ax.set_xlim(xmin, xmax)
ax.set_ylim(xmin, xmax)
ax.axis('off')

# Get the indices of the images we'll plot, sorted by zorder
inds = np.array(list(range(0, 1000, 50)), dtype = int)
inds = inds[np.argsort([-b.z[i] for i in inds])]

# Plot images at different phases
for i in inds:
  
  # The x, y, z, r position of the planet
  x = b.x[i]
  y = b.y[i]
  z = b.z[i]
  r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    
  # Rotate to a frame XYZ in which the orbital plane is the xz plane
  # The orbit is now edge-on and horizontal on the sky
  # TODO: Account for nonzero Omega
  X = x
  Y = y * np.sin(b._inc) - z * np.cos(b._inc) # = 0
  Z = z * np.sin(b._inc) + y * np.cos(b._inc)
  
  # Get the coordinates of the hotspot in the XYZ frame
  Xh = X - X * (b._r / r)
  Yh = Y - Y * (b._r / r) # = 0
  Zh = Z - Z * (b._r / r)
  
  # Add an offset?
  # TODO

  # Rotate back to the original plane
  xh = Xh
  yh = Yh * np.sin(b._inc) + Zh * np.cos(b._inc)
  zh = Zh * np.sin(b._inc) - Yh * np.cos(b._inc)

  # Compute the effective rotation angle
  rotation = -np.arctan((yh - y) / (xh - x))
  
  # Find the x coordinate of the hotspot in the rotated frame,
  # relative to the planet center
  xhrc = (xh - x) * np.cos(rotation) - (yh - y) * np.sin(rotation)
  
  # Find the effective orbital phase
  if zh >= 0:
    theta = np.arccos(-xhrc / b._r)
  else:
    theta = np.pi + np.arccos(xhrc / b._r)
  
  # Plot the radial vector
  ax.plot([0, x], [0, y], 'k-', alpha = 0.5, lw = 1)
  
  # Get the figure coordinates of the point
  disp_coords = ax.transData.transform((x, y))
  xf, yf = fig.transFigure.inverted().transform(disp_coords)
  
  # Draw the planet
  _, ax_eye, _, _ = eyeball.Draw(theta = theta, rotation = rotation, fig = fig, pos = [xf - 0.015, yf - 0.015, 0.03, 0.03])
  
  # DEBUG
  #ax_eye.plot([-1,1], [0,0], 'k-', lw = 1, alpha = 0.5)
  
# Show
pl.show()

'''
# -*-*- EXPERIMENTAL -*-*-
#
# Correct for a latitudinal hotspot offset. We are effectively
# moving the star up a distance `d` along the axis perpendicular 
# to the orbital plane.
dlat = r * np.tan(b._dlambda)
  
# Rotate to a frame in which the orbital plane is the xy plane
yp = y * np.cos(-b._inc) - z * np.sin(-b._inc)
zp = z * np.cos(-b._inc) + y * np.sin(-b._inc)
xp = x

# Rotate to a frame in which the orbital plane is the xz plane
yp = y * np.sin(b._inc) - z * np.cos(b._inc)
zp = z * np.sin(b._inc) + y * np.cos(b._inc)
xp = x

# TODO: Go from here. The above expressions are ~mostly~ correct.

# DEBUG
#x = xp
#y = yp
#z = zp
  
# -*-*-              -*-*-
'''

'''
# Find the rotation angle that would put it on the x axis
rotation = -np.arctan(yss / xss)

# Find the x coordinate in the rotated frame
xss_r = xss * np.cos(rotation) - yss * np.sin(rotation)

# Find the y coordinate in a frame where the orbital plane
# runs parallel to the x axis
yprime = y * np.cos(-b._Omega) + x * np.sin(-b._Omega)

# This is the effective orbital phase
if yprime > 0:
  theta = np.arccos(-xss_r / b._r)
else:
  theta = np.pi + np.arccos(xss_r / b._r) 
'''

'''
-*-* BEST ONE, no offset *-*-
# Get the coordinates of the hotspot relative to the planet center
xh = -x * (b._r / r)
yh = -y * (b._r / r)
zh = -z * (b._r / r)

# Find the rotation angle that would put it on the x axis
rotation = -np.arctan(yh / xh)

# Find the x coordinate in the rotated frame
xh_r = xh * np.cos(rotation) - yh * np.sin(rotation)

# Find the effective orbital phase
if zh < 0:
  theta = np.arccos(-xh_r / b._r)
else:
  theta = np.pi + np.arccos(xh_r / b._r)
'''