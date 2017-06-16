'''
'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
from planetplanet.detect import jwst
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.units as u
import numpy as np

mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rcParams['font.size'] = 25.0
mpl.rc('text', usetex=True)

#np.random.seed(1213)
np.random.seed(1234)

cadence = 1.0 # mins
Tstart = 100.0  # days
Tend = 120.0   # days
d = 12.2      # pc
saveplot = False
savetxt = False
pc_to_meter = u.pc.in_units(u.m)
oversample = True
Nos = 1

# Instantiate the Trappist-1 system
system = Trappist1(sample = True, ttvs = False, phasecurve = False, adaptive = True,
                   oversample = Nos, distance = d)

# Compute cadence in days
dt = cadence / 60. / 24.

# Define time grid
if oversample:
    time = np.arange(Tstart, Tend, dt)
else:
    time = np.linspace(Tstart, Tend, 10000)

# Get the occultation light curves
system.compute(time, lambda1 = 4, lambda2 = 30, R = 3000)

# Calculate flux at a distance, d
flux = system.flux # / (d * pc_to_meter)**2
lam = system.wavelength

# Get MIRI Filter "wheel"
wheel = jwst.get_miri_filter_wheel()

# Compute MIRI lightcurves
#jwst.lightcurves(wheel, flux, time, lam, obscad=cadence, plot=True)

# Loop over MIRI Filters
for filt in wheel:

    # Compute lightcurve in filter
    filt.compute_lightcurve(flux, time, lam, obscad=cadence, oversample=oversample)

    # Setup plot
    fig, ax = plt.subplots(figsize=(16,6))
    ax.set_title(r"%s" %filt.name)
    ax.set_ylabel("Relative Flux")
    ax.set_xlabel("Time [days]")

    # Plot lightcurve
    filt.lightcurve.plot(ax0=ax)

    # Save or show plot
    if saveplot:
        fig.savefig("../img/jwst_lc_%s.png" %filt.name, bbox_inches="tight")
    else:
        fig.subplots_adjust(bottom=0.2)

    # Save data file
    if savetxt:
        # Compose data array to save
        data = np.array([filt.lightcurve.time, filt.lightcurve.obs, filt.lightcurve.sig]).T

        # Save txt file
        np.savetxt("jwst_lc_%s_%imin.txt" %(filt.name, cadence), data, fmt=str("%.6e"),
                   header="time [days]      flux         error", comments="")

system.plot_lightcurve()

plt.show()
