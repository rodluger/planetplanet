#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
scatter.py |github|
-------------------

Computes all occultations that occur in the TRAPPIST-1 system over a 
3 year time period for a random draw from the prior. Plots each
occultation as a circle in a top-view of the system; the circle size, 
transparency, and color indicate the duration, impact parameter, and 
occulting body, respectively.

  .. plot::
     :align: center
     
     from planetplanet import Trappist1
     import matplotlib.pyplot as pl
     system = Trappist1(sample = True, nbody = True)
     system.scatter_plot(0, 365 * 3)
     pl.show()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/scatter.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Trappist1
import matplotlib.pyplot as pl
import numpy as np

if __name__ == '__main__':

  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True, nbody = True)
  fig = system.scatter_plot(0, 365 * 3)
  fig.savefig("scatter.pdf", bbox_inches = 'tight')
  pl.show()