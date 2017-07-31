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
b = Planet('b', per = 10., inc = 80., Omega = 30., t0 = 0, ecc = 0., w = 0., dlambda = 0, dpsi = 0)
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
  
  # The position of the planet in the original frame (frame #0)
  x0 = b.x[i]
  y0 = b.y[i]
  z0 = b.z[i]
  r = np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
    
  # Rotate so that the orbital plane is parallel to the x axis (frame #1)
  x1 = x0 * np.cos(b._Omega) + y0 * np.sin(b._Omega)
  y1 = y0 * np.cos(b._Omega) - x0 * np.sin(b._Omega)
  z1 = z0
  
  # Now rotate so that it's also parallel to the z axis (frame #2)
  x2 = x1
  y2 = 0
  z2 = z1 * np.sin(b._inc) + y1 * np.cos(b._inc)
  
  # Finally, rotate so that the substellar point faces the observer (frame #3)
  # (TODO)
  x3 = x2
  y3 = y2
  z3 = z2
  
  # Get the coordinates of the hotspot in this frame
  xh3 = x3 - x3 * (b._r / r)
  yh3 = 0
  zh3 = z3 - z3 * (b._r / r)
  
  # Get the coordinates relative to the planet center
  dxh3 = xh3 - x3
  dyh3 = yh3 - y3
  dzh3 = zh3 - z3
  
  # Add the offset in longitude. This is a counterclockwise rotation
  # in the xz plane by an angle `dpsi` about the planet center
  dxh3, dzh3 = dxh3 * np.cos(b._dpsi) - dzh3 * np.sin(b._dpsi), dzh3 * np.cos(b._dpsi) + dxh3 * np.sin(b._dpsi)
  
  # Add the offset in latitude. This is trickier (TODO)
  
  # Get the coordinates relative to the star
  xh3 = x3 + dxh3
  yh3 = y3 + dyh3
  zh3 = z3 + dzh3
  
  # Rotate back to frame #2 (TODO)
  xh2 = xh3
  yh2 = yh3
  zh2 = zh3
  
  # Rotate back to frame #1
  xh1 = xh2
  yh1 = yh2 * np.sin(b._inc) + zh2 * np.cos(b._inc)
  zh1 = zh2 * np.sin(b._inc) - yh2 * np.cos(b._inc)
  
  # Finally, rotate back to the original frame (frame #0)
  xh0 = xh1 * np.cos(b._Omega) - yh1 * np.sin(b._Omega)
  yh0 = yh1 * np.cos(b._Omega) + xh1 * np.sin(b._Omega)
  zh0 = zh1

  # Compute the effective rotation angle
  rotation = -np.arctan((yh0 - y0) / (xh0 - x0))
  
  # Find the x coordinate of the hotspot in the rotated frame (r)
  # relative to the planet center
  dxhr = (xh0 - x0) * np.cos(rotation) - (yh0 - y0) * np.sin(rotation)
  
  # Prevent tiny numerical errors from yielding NaNs in the arccos below
  if dxhr < -1: 
    dxhr = -1
  elif dxhr > 1: 
    dxhr = 1
  
  # Find the effective phase angle
  if zh0 >= 0:
    theta = np.arccos(-dxhr / b._r)
  else:
    theta = np.pi + np.arccos(dxhr / b._r)
  
  # Plot the radial vector
  ax.plot([0, x0], [0, y0], 'k-', alpha = 0.5, lw = 1)
  
  # Get the figure coordinates of the point
  disp_coords = ax.transData.transform((x0, y0))
  xf, yf = fig.transFigure.inverted().transform(disp_coords)
  
  # Draw the planet
  _, ax_eye, _, _ = eyeball.Draw(theta = theta, rotation = rotation, fig = fig, pos = [xf - 0.015, yf - 0.015, 0.03, 0.03])
  
# Show
pl.show()