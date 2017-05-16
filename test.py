import planetplanet as pp
import numpy as np
import matplotlib.pyplot as pl
AUREARTH = 23454.9271
color = pl.get_cmap('RdBu_r')
TOL = 1e-10

time = np.linspace(0., 30., 1000000)


T1b = pp.Orbit(time = time, inc = 89.65, a = 0.01111 * AUREARTH, per = 1.51087081, r = 1.086, ecc = 0.5, w = np.pi / 8, t0 = 7322.51736)
T1c = pp.Orbit(time = time, inc = 89.67, a = 0.01521 * AUREARTH, per = 2.4218233, r = 1.056, ecc = 0.5, w = np.pi / 4, t0 = 7282.80728)

# We need to flip the problem when the right side of the planet
# is being illuminated; this is just a reflection about the y axis.
# NOTE: for non-transiting planets, we need to actually rotate the planet,
# since theta is no longer equal to the orbital angle 
neginds = np.where(T1b.x < 0)[0]
T1b.x[neginds] *= -1
T1c.x[neginds] *= -1

# Compute the orbital angle (for an edge-on orbit this is really 
# just the true anomaly, but with a phase offset)
theta = np.arctan(T1b.y / (TOL + T1b.x))

# We compute our light curves with the occultor at the origin, so
# let's find the relative position of the occulted planet and
# the radial separation of the two
x0 = T1b.x - T1c.x
y0 = T1b.y - T1c.y
sep = np.sqrt(x0 ** 2 + y0 ** 2)

# We only care about times when the planets are overlapping; get
# those indices
inds = np.where(sep <= (T1b.r + T1c.r))[0]
lims = np.where(np.diff(inds) > 1)[0]
lims = np.concatenate(([-1], lims, [len(inds) - 1]))
inds = [inds[a:b] for a, b in zip(lims[:-1] + 1, lims[1:] + 1)]

# Instantiate the eyeballs
occultor = pp.eyeball.Occultor(T1c.r)
occulted = pp.eyeball.Occulted(0, 0, T1b.r, 0, occultor, midnight = 0)

f = [np.empty(len(n)) for n in inds]
t = [np.empty(len(n)) for n in inds]
for i, n in enumerate(inds):
  for j, m in enumerate(n):

    occulted.x0 = x0[m]
    occulted.y0 = y0[m]
    occulted.theta = theta[m]
    f[i][j] = occulted.delta_flux
    t[i][j] = time[m]

dt = 0
dx = 0.001
for ti, fi in zip(t, f):
  x = ti - ti[0] + dt
  x = np.concatenate(([x[0] - 2 * dx, x[0] - dx], x, [x[-1] + dx, x[-1] + 2 * dx]))
  fi = np.concatenate(([0, 0], fi, [0, 0]))
  pl.plot(x, fi, 'k-')
  dt = x[-1] + 2 * dx
pl.show()

