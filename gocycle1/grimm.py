"""Test using the new Grimm et al. (2018) data."""
from planetplanet import Trappist1, TRAPPIST_T0
from planetplanet.constants import MINUTE
import numpy as np
import matplotlib.pyplot as pl


def expected(actual, t0, P):
    """Return the expected time of transit."""
    periodic = np.array(t0 + np.arange(-1000, 1000) * P)
    actual = np.atleast_1d(actual)
    e = np.zeros_like(actual)
    for i, t in enumerate(actual):
        e[i] = periodic[np.argmin(np.abs(periodic - t))]
    return e


# Plot the data
fig, ax = pl.subplots(2, 4)
ax = ax.flatten()
periods = [1.51087081, 2.4218233, 4.049610,
           6.099615, 9.206690, 12.35294, 18.767]
t0 = [None for i in range(7)]
P = [None for i in range(7)]
for i, name in enumerate(['b', 'c', 'd', 'e', 'f', 'g', 'h']):

    # Setup
    ax[i].set_title(name)

    # Load
    data_times, data_err = np.loadtxt(name + ".ttv", unpack=True)
    data_times -= 2450000

    # Time of first transit and period, taken from the *data*
    t0[i] = data_times[0]
    data_tnum = np.zeros_like(data_times)
    for j in range(1, len(data_times)):
        data_tnum[j] = data_tnum[j - 1] + \
            np.round((data_times[j] - data_times[j - 1]) /
                     periods[i])
    P[i], _ = np.polyfit(data_tnum, data_times, 1)

    # Plot it
    data_ttv = data_times - expected(data_times, t0[i], P[i])
    ax[i].errorbar(data_times,
                   data_ttv / MINUTE,
                   data_err / MINUTE, fmt='o', zorder=1)


# N-body
cadence = 1 * MINUTE
time = np.arange(TRAPPIST_T0, TRAPPIST_T0 + 800, cadence)

# Integrate
system = Trappist1(seed=12345, nbody=True, flat=True, timestep=MINUTE)
for body in system.bodies[1:]:
    body.host = system.bodies[0]
system.compute_orbits(time)
for i, planet in enumerate(system.bodies[1:]):

    # Get transit times
    inds = np.where((np.array(planet.x[:-1] > 0, dtype=int) -
                    (np.array(planet.x[1:] < 0, dtype=int) == 0))
                    & (planet.z[:-1] > 0))[0]
    nbody_times = planet.time[inds]

    # Plot
    nbody_ttv = nbody_times - expected(nbody_times, t0[i], P[i])
    ax[i].fill_between(nbody_times,
                       (nbody_ttv - cadence) / MINUTE,
                       (nbody_ttv + cadence) / MINUTE,
                       alpha=0.5, color='lightgray', zorder=0)


# Plot the TTV curves
pl.show()
