#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
occultation.py |github|
-----------------------

A simple occultation of TRAPPIST-1c by TRAPPIST-1b.
Planet c has a latitudinal hotspot offset, just for fun.


  .. plot::
     :align: center
     
     from scripts import occultation
     occultation._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/occultation.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
from planetplanet.constants import *
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np

def _test():
  '''
  
  '''
  
  plot()

def plot():
  '''
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True, phasecurve = True, nbody = True, seed = 999)
  
  # Fudge: Let's make this a nice, near-full occultation
  system.c.Omega = -0.15
  
  # Give `c` a large latitudinal offset in its hotspot just for fun
  system.c.Phi = 30
  
  # Compute an occultation by `b`
  # This would be on December 3, 2021
  time = np.linspace(9552.9364, 9552.9664, 100)
  system.compute(time)

  # Plot the occultation
  fig, axlc, axxz, axim = system.plot_occultation('c', 9552.95)
  pl.show()

if __name__ == '__main__':
  plot()
  