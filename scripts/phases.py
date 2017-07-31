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
b = Planet('b', per = 10., inc = 80., Omega = 30., t0 = 0, ecc = 0., w = 0., dlambda = 0, dpsi = 30)
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
  
  # The position of the planet in the xyz frame
  x = b.x[i]
  y = b.y[i]
  z = b.z[i]
  r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    
  # Rotate to a frame XYZ in which the orbital plane is the xz plane.
  # First rotate so that the orbital plane is parallel to the x axis.
  # This is the x_y_z_ frame
  x_ = x * np.cos(b._Omega) + y * np.sin(b._Omega)
  y_ = y * np.cos(b._Omega) - x * np.sin(b._Omega)
  z_ = z
  
  # Now rotate so that it's also parallel to the z axis to get to the
  # XYZ frame, where we will perform our hotspot shift(s)
  X = x_
  Y = 0
  Z = z_ * np.sin(b._inc) + y_ * np.cos(b._inc)
  
  # Get the coordinates of the hotspot in the XYZ frame
  Xh = X - X * (b._r / r)
  Yh = 0
  Zh = Z - Z * (b._r / r)
  
  # Get the coordinates relative to the planet center
  dXh = Xh - X
  dYh = Yh - Y
  dZh = Zh - Z
  
  # Add the offset in longitude. This is a counterclockwise rotation
  # in the xz plane by an angle `dpsi` about the planet center
  dXh, dZh = dXh * np.cos(b._dpsi) - dZh * np.sin(b._dpsi), dZh * np.cos(b._dpsi) + dXh * np.sin(b._dpsi)
  
  # Add the offset in latitude. This is trickier. TODO
  
  # Get the coordinates relative to the star
  Xh = X + dXh
  Yh = Y + dYh
  Zh = Z + dZh

  # Rotate back to the x_y_z_ frame
  xh_ = Xh
  yh_ = Yh * np.sin(b._inc) + Zh * np.cos(b._inc)
  zh_ = Zh * np.sin(b._inc) - Yh * np.cos(b._inc)
  
  # Now rotate back to the original xyz frame
  xh = xh_ * np.cos(b._Omega) - yh_ * np.sin(b._Omega)
  yh = yh_ * np.cos(b._Omega) + xh_ * np.sin(b._Omega)
  zh = zh_

  # Compute the effective rotation angle
  rotation = -np.arctan((yh - y) / (xh - x))
  
  # Find the x coordinate of the hotspot in the rotated frame,
  # relative to the planet center
  xhrc = (xh - x) * np.cos(rotation) - (yh - y) * np.sin(rotation)
  
  # Prevent tiny numerical errors from yielding NaNs in the arccos below
  if xhrc < -1: 
    xhrc = -1
  elif xhrc > 1: 
    xhrc = 1
  
  # Find the effective phase angle
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
  
# Show
pl.show()