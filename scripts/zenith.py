import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
cmap = pl.get_cmap('RdBu')

def ZenithAngle(x, y, r, theta):
  '''
  Compute the zenith angle.
  
  '''
  
  # Normalize
  x = x / r
  y = y / r
  x2 = x * x
  y2 = y * y
  
  # This is a solution to a quadratic equation in z = sin(za) **  2 
  z = 0.5 * ((1 - 2 * x2 - y2) * np.cos(2 * theta) + 2 * x * np.sqrt(1 - x2 - y2) * np.sin(2 * theta) + y2 + 1)
  
  # Where are we relative to the terminator?
  xterm = np.sin(theta) * np.sqrt(np.abs(1 - y2))
  if (x <= xterm):
    return np.arcsin(np.sqrt(z))
  else:
    return np.pi - np.arcsin(np.sqrt(z))

fig, ax = pl.subplots(1)
ax.axis('off')
fig.subplots_adjust(bottom = 0.2)

z = np.zeros((100, 100)) * np.nan
img = pl.imshow(z, cmap = cmap, vmin = 0, vmax = 180., extent = (-1, 1, -1, 1))
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
        z[j,i] = ZenithAngle(x, y, 1, theta * np.pi / 180) * 180 / np.pi
  img.set_data(z)
  fig.canvas.draw_idle()
  
slider.on_changed(update)
update(45.)

pl.show()