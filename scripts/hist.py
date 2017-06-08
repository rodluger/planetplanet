#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py
-------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import matplotlib.pyplot as pl
import numpy as np
import corner
from tqdm import tqdm

# Number of prior samples
nsamp = 100

# Minimum duration (minutes)
mind = 10.

# Maximum impact parameter
maxb = 0.5

# Draw samples from the prior
hist = [[] for k in range(7)]
count = [np.zeros(nsamp) for k in range(7)]
for n in tqdm(range(nsamp)):

  # Instantiate the Trappist-1 system
  system = Trappist1(uncertainty = True, ttvs = False, quiet = True)
  system.settings.dt = 0.0001
  h = system.histogram(0, 365)
  
  # Loop over the planets
  for k in range(7):
  
    # Count the number of events longer than `mind`
    # and with impact parameter below `minb`
    if len(h[k]):
      count[k][n] = len(np.where((h[k][:,1] < maxb) & (h[k][:,2] > np.log10(mind)))[0])
    
    # Append to cumulative histogram
    hist[k].extend(h[k])
    
# Plot!
for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):
  
  # Frequency histogram
  fig = corner.corner(count[k])
  fig.subplots_adjust(left = 0.15, right = 0.9, bottom = 0.2, top = 0.9)
  ax = fig.axes[0]
  ax.set_xlabel('Occultations per year', fontsize = 12, fontweight = 'bold')
  ax.set_ylabel('Probability', fontsize = 12, fontweight = 'bold')
  fig.savefig('%s.hist.png' % planet, bbox_inches = 'tight')
  pl.close()
  
  # Corner plot
  samples = np.array(hist[k])
  fig = corner.corner(samples, range = [(-190,190), (-0.1,1.1), (-0.1, 3.1)], labels = ["Phase [deg]", "Impact parameter", "Duration [min]"])
  for i, ax in enumerate(fig.axes):
    ax.set_xlabel(ax.get_xlabel(), fontsize = 10, fontweight = 'bold')
    ax.set_ylabel(ax.get_ylabel(), fontsize = 10, fontweight = 'bold')
  fig.axes[0].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[3].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[6].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[3].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
  fig.axes[7].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
  fig.axes[8].set_xticks([0, 1, 2, 3])
  fig.axes[8].set_xticklabels([1, 10, 100, 1000])
  fig.axes[6].set_yticks([0, 1, 2, 3])
  fig.axes[6].set_yticklabels([1, 10, 100, 1000])
  fig.savefig('%s.corner.png' % planet, bbox_inches = 'tight')
  pl.close()