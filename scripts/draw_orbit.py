#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
draw_orbit.py |github|
----------------------

Draws the orbit of TRAPPIST-1 c on the sky. The inclination is set to 60 degrees
for better visualization.

  .. plot::
     :align: center
     
     from scripts import draw_orbit
     draw_orbit._test()
  
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/draw_orbit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
from planetplanet import maps
import matplotlib.pyplot as pl

def _test():
  '''
  
  '''
  
  plot()

def plot():
  
  # Instantiate the TRAPPIST-1 system
  system = Trappist1(sample = True, phasecurve = True, nbody = True, seed = 999)
  
  # Give it a non-edge-on inclination
  system.c.inc = 60
  
  # Draw the orbit
  system.c.draw_orbit(draw_outline = True, size = 2, nz = 11, plot_phasecurve = False)
  
  # Show!
  pl.show()

if __name__ == '__main__':
  plot()