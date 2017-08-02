#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball_orbit.py
----------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import eyeball
import matplotlib.pyplot as pl

# Plot it!
fig1, ax1, fig2, ax2 = eyeball.DrawOrbit(inc = 60., Omega = 0., ecc = 0.5, w = 0, size = 2, 
                                         plot_phasecurve = True, label_phases = False, rasterize = True)

# Tweak the appearance a little for the PDFs
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.margins(0.01, 0.2)
fig1.savefig('../img/eyeball_orbit.pdf', bbox_inches = 'tight', dpi = 600)
fig2.savefig('../img/eyeball_phasecurve.pdf', bbox_inches = 'tight')

# Show
pl.show()