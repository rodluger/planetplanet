import numpy as np
import matplotlib.pyplot as pl


semis = np.array([11.11, 15.21, 21.44, 28.17, 37.1, 45.1, 60])
radii = np.array([1.086, 1.056, 0.772, 0.918, 1.045, 1.127, 0.755]) / (23.455)
incs = np.array([89.65, 89.67, 89.75, 89.86, 89.680, 89.710, 89.80]) * np.pi / 180
colors = ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6', '#999999']


fig = pl.figure(figsize = (5,5))

for i in range(7):
  x = np.linspace(-semis[i], semis[i], 1000)
  a = semis[i]
  b = semis[i] * np.cos(incs[i])
  y = (b / a) * np.sqrt(a ** 2 - x ** 2)
  
  pl.plot(x, y, lw = 1, color = colors[i])
  pl.plot(x, -y, lw = 1, color = colors[i])
  
  pl.fill_between(x, y - radii[i], y + radii[i], color = colors[i], alpha = 0.3)
  pl.fill_between(x, -y - radii[i], -y + radii[i], color = colors[i], alpha = 0.3)
pl.show()