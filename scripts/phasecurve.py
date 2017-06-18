import numpy as np
from scipy.integrate import dblquad
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from tqdm import tqdm
cmap = pl.get_cmap('RdBu_r')
SBOLTZ = 5.670367e-8
KBOLTZ = 1.38064852e-23    
HPLANCK = 6.62607004e-34   
CLIGHT = 2.998e8           
MICRON = 1.e-6

def Temperature(x, y, theta, mult, tnight):
  '''
  Compute the temperature.
  
  '''
  
  # Normalize
  x2 = x * x
  y2 = y * y
  
  # This is a solution to a quadratic equation in z = sin(za) **  2 
  z = 0.5 * ((1 - 2 * x2 - y2) * np.cos(2 * theta) + 2 * x * np.sqrt(1 - x2 - y2) * np.sin(2 * theta) + y2 + 1)
  
  # The temperature
  T = mult * (1 - np.atleast_1d(z)) ** 0.125
  
  # Where are we relative to the terminator?
  xterm = np.sin(theta) * np.sqrt(np.abs(1 - y2))
  
  # Correct the nightside temperature
  T[(x > xterm) | (T < tnight)] = tnight

  return T
  
def Radiance(T, w):
  '''
  
  '''
  
  a = 2 * HPLANCK * CLIGHT * CLIGHT / (w * w * w * w * w)
  b = HPLANCK * CLIGHT / (w * KBOLTZ * T)
  return a / (np.exp(b) - 1)
  
def Interact(A = 0.3, Fstar = 1361.):
  '''
  
  '''
  
  mult = (Fstar * (1 - A) / SBOLTZ) ** 0.25

  fig, ax = pl.subplots(1)
  ax.axis('off')
  fig.subplots_adjust(bottom = 0.2)

  z = np.zeros((100, 100)) * np.nan
  img = pl.imshow(z, cmap = cmap, vmin = 0, vmax = 400., extent = (-1, 1, -1, 1))
  x = np.linspace(-0.99,0.99,1000)
  ax.plot(x, np.sqrt(0.99 ** 2 - x ** 2), 'k-', lw = 2)
  ax.plot(x, -np.sqrt(0.99 ** 2 - x ** 2), 'k-', lw = 2)
  ax.set_xlim(-1.1,1.1)
  ax.set_ylim(-1.1,1.1)

  axslider = pl.axes([0.3, 0.05, 0.44, 0.03])
  slider = Slider(axslider, r'$\theta$', -90., 90., valinit = 45.)

  def update(val):
    theta = slider.val
    for i, x in enumerate(np.linspace(-1,1,100)):
      for j, y in enumerate(np.linspace(-1,1,100)):
        if (x ** 2 + y ** 2 <= 1):
          z[j,i] = Temperature(x, y, theta * np.pi / 180, mult, 100)
    img.set_data(z)
    fig.canvas.draw_idle()
  
  slider.on_changed(update)
  update(45.)

  pl.show()

def Phasecurve(ntheta = 30, tnight = 100, A = 0.3, Fstar = 1361., nx = 100, ny = 50,
               lambda1 = 5., lambda2 = 15., nlam = 10, **kwargs):
  '''
  
  '''
  
  mult = (Fstar * (1 - A) / SBOLTZ) ** 0.25
  theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta, endpoint = True)
  w = np.linspace(lambda1, lambda2, nlam, endpoint = True) * MICRON
  phasecurve = np.zeros((ntheta, nlam))
  
  # Compute
  for n in tqdm(range(ntheta)):
    rad = np.zeros((nx, ny, nlam))
    for i, x in enumerate(np.linspace(-1, 1, nx)):
      for j, y in enumerate(np.linspace(0, 1, ny)):
        if (x ** 2 + y ** 2 <= 1):
          T = Temperature(x, y, theta[n], mult, tnight)
          rad[i,j] = Radiance(T, w)
    phasecurve[n] = 2 * np.sum(rad, axis = (0,1))

  # Plot
  for k in range(nlam):
    pl.plot(theta * 180 / np.pi, phasecurve[:,k])
  pl.xlabel("Orbital angle [deg]", fontweight = 'bold', fontsize = 16)
  pl.ylabel("Spectral Radiance", fontweight = 'bold', fontsize = 16)
  pl.show()
  
Phasecurve()
   