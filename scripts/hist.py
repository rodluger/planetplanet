#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
hist.py
-------

Histograms of the occultation events as a function of phase, duration, and
impact parameter for each of the seven TRAPPIST-1 planets.

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
nsamp = 3000

# Minimum duration (minutes)
mind = 10.

# Maximum impact parameter
maxb = 0.5

# Compute (or load saved)?
compute = True

if compute:

  # Draw samples from the prior
  hist = [[] for k in range(7)]
  count = [np.zeros(nsamp) for k in range(7)]
  for n in tqdm(range(nsamp)):

    # Instantiate the Trappist-1 system
    system = Trappist1(sample = True, nbody = False, quiet = True)
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

  # Save
  np.savez('hist.npz', hist = hist, count = count)

else:

  data = np.load('hist.npz')
  hist = data['hist']
  count = data['count']

# Corner plot
for k, planet in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):  
  samples = np.array(hist[k])
  fig = corner.corner(samples, data_kwargs = {'alpha': 0.01}, range = [(-180,180), (0,1), (0, 3)], labels = ["Phase [deg]", "Impact parameter", "Duration [min]"])
  for i, ax in enumerate(fig.axes):
    ax.set_xlabel(ax.get_xlabel(), fontsize = 10, fontweight = 'bold')
    ax.set_ylabel(ax.get_ylabel(), fontsize = 10, fontweight = 'bold')
  fig.axes[0].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[3].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[6].set_xticks([-180, -90, 0, 90, 180])
  fig.axes[3].set_yticks([0.2, 0.4, 0.6, 0.8])
  fig.axes[7].set_xticks([0.2, 0.4, 0.6, 0.8])
  fig.axes[8].set_xticks([0, 1, 2, 3])
  fig.axes[8].set_xticklabels([1, 10, 100, 1000])
  fig.axes[6].set_yticks([0, 1, 2, 3])
  fig.axes[6].set_yticklabels([1, 10, 100, 1000])
  fig.savefig('%s.corner.png' % planet, bbox_inches = 'tight')
  pl.close()
  
# Frequency histogram
fig, ax = pl.subplots(1, figsize = (8, 6))
system = Trappist1(sample = True, nbody = False, quiet = True)
for k, planet in enumerate(system.bodies[1:]): 
  ax.hist(count[k], color = planet.color, edgecolor = 'none', alpha = 0.25, histtype = 'stepfilled', normed = True, range = (0,50), zorder = 0, label = planet.name)
  ax.hist(count[k], color = planet.color, histtype = 'step', normed = True, range = (0,50), zorder = 1, lw = 2)
leg = ax.legend(loc = 'upper right', title = "Planet")
ax.get_legend().get_title().set_fontweight('bold')
ax.set_xlabel('Occultations per year', fontsize = 12, fontweight = 'bold')
ax.set_ylabel('Probability', fontsize = 12, fontweight = 'bold')
fig.savefig('hist.png', bbox_inches = 'tight')
pl.close()