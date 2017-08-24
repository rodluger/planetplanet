#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
plot_orbits.py |github|
-----------------------

Plot the TRAPPIST-1 planet orbits from different viewing angles.

  .. plot::
     :align: center
     
     from scripts import plot_orbits
     import matplotlib.pyplot as pl
     plot_orbits.plot()
     pl.show()
     
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/plot_orbits.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import planetplanet as pp
import matplotlib.pyplot as pl
import numpy as np

def plot():
  '''
  
  '''
  
  # Instantiate the TRAPPIST-1 system
  time = np.linspace(0, 20, 10000)
  system = pp.Trappist1(sample = True)
  
  # Compute and plot the orbits
  system.compute(time)
  fig, ax = system.plot_orbits(cmap = 'jet_r')
  
  return fig, ax

if __name__ == '__main__':
  fig, ax = plot()
  pl.show()