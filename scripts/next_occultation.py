#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
next_occultation.py |github|
----------------------------

Compute the time of the next occultation of a given planet
and plot the light curve.

  .. plot::
     :align: center
     
     from scripts import next_occultation
     import matplotlib.pyplot as pl
     next_occultation.plot()
     pl.show()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/next_occultation.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(42)

def plot():
  '''
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True)

  # Get the next occultation of b
  t = system.next_occultation(999, system.c, occultor = system.b)

  # Get the light curve around that point
  time = np.linspace(t - 0.05, t + 0.05, 1000)

  system.compute(time)
  fig, ax = system.plot_lightcurve()
  return fig, ax
  
if __name__ == '__main__':
  fig, ax = plot()
  pl.show()