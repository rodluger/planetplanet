import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
import numpy as np
color = pl.get_cmap('plasma_r')

def EllipseUpper(x, x0, y0, a, b):
  '''
  The function describing the upper half of an ellipse
  
  '''
  
  A = b ** 2 - (x - x0) ** 2
  if np.abs(A) < 1e-15:
    A = 0
  return y0 + a / b * np.sqrt(A)

def EllipseLower(x, x0, y0, a, b):
  '''
  The function describing the lower half of an ellipse
  
  '''
  
  A = b ** 2 - (x - x0) ** 2
  if np.abs(A) < 1e-15:
    A = 0
  return y0 - a / b * np.sqrt(A)

def Roots(r, xE, yE, a, b):
  '''
  Returns the points of intersection of an ellipse centered at (`xE`, `yE`) 
  with radii `b` and `a` and a circle of radius `r` centered at the origin
  
  '''
  
  # Get the coefficients
  A = a ** 2 / b ** 2 - 1
  B = -2 * a ** 2 * xE / b ** 2
  C = r ** 2 - yE ** 2 - a ** 2 + a ** 2 / b ** 2 * xE ** 2
  D = 4 * yE ** 2 * a ** 2 / b ** 2
  c4 = A ** 2
  c3 = 2 * A * B
  c2 = 2 * A * C + B ** 2 + D
  c1 = 2 * B * C - 2 * D * xE
  c0 = C ** 2 - (b ** 2 - xE ** 2) * D
  
  # Solve numerically
  roots = np.roots([c4, c3, c2, c1, c0])
  
  # Filter out imaginary roots
  good_roots = []
  for x in roots:
    if (x.imag == 0):
      good_roots.append(x.real)
      
  return good_roots

def Ellipse(l, rp = 1, xp = 0, yp = 0, theta = np.pi / 8, ro = 1, n = 1000):
  '''
  
  '''
  
  # Ellipse parameters
  a = rp * np.abs(np.sin(l))
  b = a * np.abs(np.sin(theta))
  x0 = xp - rp * np.cos(l) * np.cos(theta)
  y0 = yp
  
  # Identify and remove points behind the limb
  d = rp * np.cos(l) * np.sin(theta) * np.tan(theta)
  x = np.linspace(x0 - b, x0 + b, n)
  if theta > 0:
    x[x < x0 - d] = np.nan
  else:
    x[x > x0 - d] = np.nan
    
  # The ellipse equation
  B = b ** 2 - (x - x0) ** 2
  B[B<0] = 0
  ybot = y0 - a / b * np.sqrt(B)
  ytop = y0 + a / b * np.sqrt(B)
  
  # Compute important vertices
  vertices = []
  
  # Ellipse-circle intersections
  if b ** 2 >= d ** 2:
    vertices.append((x0 - d, y0 - a / b * np.sqrt(b ** 2 - d ** 2)))
    vertices.append((x0 - d, y0 + a / b * np.sqrt(b ** 2 - d ** 2)))
  
  # Ellipse x minimum
  if (theta > 0) and (x0 - b > x0 - d):
    vertices.append((x0 - b, y0))
  elif (theta <= 0) and (x0 - b < x0 - d):
    vertices.append((x0 - b, y0))
  
  # Ellipse x maximum
  if (theta > 0) and (x0 + b > x0 - d):
    vertices.append((x0 + b, y0))
  elif (theta <= 0) and (x0 + b < x0 - d):
    vertices.append((x0 + b, y0))
  
  # Compute the ellipse-occultor intersection points
  roots = Roots(ro, x0, y0, a, b)
  for root in roots:
    
    if ((theta > 0) and (root > x0 - d)) or ((theta <= 0) and (root < x0 - d)):

      eup = EllipseUpper(root, x0, y0, a, b)
      elo = EllipseLower(root, x0, y0, a, b)
      oup = EllipseUpper(root, 0, 0, ro, ro)
      olo = EllipseLower(root, 0, 0, ro, ro)
    
      if (np.abs(eup - oup) < 1e-7) or (np.abs(eup - olo) < 1e-7):
        vertices.append((root, eup))
      elif (np.abs(elo - oup) < 1e-7) or (np.abs(elo - olo) < 1e-7):
        vertices.append((root, elo))
    
  return (x, ybot, ytop, vertices)

# Set up the figure
fig, ax = pl.subplots(1, figsize = (7,7))
ax.set_xlim(-3,3)
ax.set_ylim(-3,3)
axslider = pl.axes([0.15, 0.05, 0.7, 0.03])
slider = Slider(axslider, r'$\theta$', -np.pi / 2, np.pi / 2, valinit = 0)

# The latitude grid
latitude = np.linspace(0, np.pi / 2, 5)[1:]



# Occultor
xo = 0
yo = 0
ro = 1
vo = [(xo - ro, yo), (xo + ro, yo)]

# Planet (occulted)
xp = -1
yp = -0.5
rp = 1.25
vp = [(xp - rp, yp), (xp + rp, yp)]

# Arrays for plotting
xparr = np.linspace(xp - rp, xp + rp, 1000)
yparrbot = yp - np.sqrt(rp ** 2 - (xparr - xp) ** 2)
yparrtop = yp + np.sqrt(rp ** 2 - (xparr - xp) ** 2)
xoarr = np.linspace(xo - ro, xo + ro, 1000)
yoarrbot = yo - np.sqrt(ro ** 2 - (xoarr - xo) ** 2)
yoarrtop = yo + np.sqrt(ro ** 2 - (xoarr - xo) ** 2)

def update(val):
  '''
  
  '''
  
  # Set up plot
  ax.clear()
  ax.set_xlim(xp - 1.25 * rp, xp + 1.25 * rp)
  ax.set_ylim(yp - 1.25 * rp, yp + 1.25 * rp)
  
  # Plot the planet and its vertices
  ax.plot(xparr, yparrbot, color = 'k')
  ax.plot(xparr, yparrtop, color = 'k')
  
  # Compute and plot lines of constant surface brightness
  vertices = []
  theta = slider.val
  if theta == 0:
    theta = 1e-7
  for l in latitude:
  
    # Compute the ellipse
    xearr, yearrbot, yearrtop, ve = Ellipse(l, rp = rp, xp = xp, yp = yp, theta = theta, ro = ro)
    ax.plot(xearr, yearrbot, lw = 2, color = color(l / np.pi))
    ax.plot(xearr, yearrtop, lw = 2, color = color(l / np.pi))
    vertices.extend(ve)
    
  # Plot the occultor
  ax.plot(xoarr, yoarrbot, color = 'k')
  ax.plot(xoarr, yoarrtop, color = 'k')
  
  # Compute the planet-occultor intersection points, if they exist
  d = np.sqrt((xo - xp) ** 2 + (yo - yp) ** 2)
  if d <= (ro + rp):
    y = np.sqrt((-d + rp - ro) * (-d - rp + ro) * (-d + rp + ro) * (d + rp + ro)) / (2 * d)
    x = (d ** 2 - rp ** 2 + ro ** 2) / (2 * d)
    cost = (xo - xp) / d
    sint = (yo - yp) / d
    x1 = -x * cost + y * sint
    y1 = -x * sint - y * cost
    x2 = -x * cost - y * sint
    y2 = -x * sint + y * cost
    vertices.append((x1, y1))
    vertices.append((x2, y2))



  # Combine, sort, and plot the vertices
  vertices.extend(vp)
  vertices.extend(vo)
  vertices = np.array(sorted(list(set(vertices))))
  for x, y in vertices:
    ax.plot(x, y, 'ro', ms = 4)
      
  fig.canvas.draw_idle()
slider.on_changed(update)
update(0)

pl.show()
