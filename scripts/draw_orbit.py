#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
draw_orbit.py |github|
----------------------

Draws the orbit of TRAPPIST-1 c on the sky. The inclination is set to 60 degrees
for better visualization.

  .. plot::
     :align: center
     
     from planetplanet.photo import Trappist1
     import matplotlib.pyplot as pl
     system = Trappist1(sample = True, phasecurve = False, nbody = True, seed = 999)
     system.c.inc = 60
     system.c.draw_orbit(draw_outline = False, size = 1.5, nz = 51)
     pl.show()
  
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/draw_orbit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
from planetplanet import maps
import matplotlib.pyplot as pl

if __name__ == '__main__':
  
  # Instantiate the TRAPPIST-1 system
  system = Trappist1(sample = True, phasecurve = False, nbody = True, seed = 999)
  
  # Give it an inclination for visualization purposes
  system.c.inc = 0
  
  system.c.Phi = 30
  system.c.radiancemap = maps.BandedCloudsMap()
  
  # Draw the orbit
  system.c.draw_orbit(draw_outline = False, size = 1.5, nz = 51)
  
  # Show!
  pl.show()

