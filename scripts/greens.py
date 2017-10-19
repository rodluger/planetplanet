import numpy as np
import matplotlib.pyplot as pl

def Circle(x, r, xC, yC):
    '''
    
    '''
    
    A = r ** 2 - (x - xC) ** 2
    if hasattr(A, '__len__'):
        A[A < 0] = 0
        A = np.sqrt(A)
    else:
        if A < 0:
            A = 0
        else:
            A = np.sqrt(A)
    return yC + np.sign(r) * A

def Intersection(xo, yo, ro, tol = 1e-5):
    '''
    Points of intersection between occulted and occultor.
    
    '''
    
    ccroots = []
    d = np.sqrt(xo ** 2 + yo ** 2)
    if (d < (1 + ro)):
        A = (-d + 1 - ro) * (-d - 1 + ro) * (-d + 1 + ro) * (d + 1 + ro)        
        if (A >= 0):
            y = np.sqrt(A) / (2 * d)
            x = -np.sqrt(1 - y ** 2)
            if np.abs(xo) < tol:
                cost = 0
                sint = 1
            else:
                frac = yo / xo
                cost = 1 / np.sqrt(frac ** 2 + 1)
                sint = frac * cost
            if (xo < 0):
                cost *= -1
                sint *= -1
            x1 = x * cost + y * sint
            x2 = x * cost - y * sint
            x3 = -x * cost + y * sint
            x4 = -x * cost - y * sint
            for x in np.unique([x1, x2, x3, x4]):
                if (np.abs(Circle(x, 1, 0, 0) - Circle(x, ro, xo, yo)) < tol):
                    y = Circle(x, 1, 0, 0)
                elif (np.abs(Circle(x, 1, 0, 0) - Circle(x, -ro, xo, yo)) < tol):
                    y = Circle(x, 1, 0, 0)
                elif (np.abs(Circle(x, -1, 0, 0) - Circle(x, ro, xo, yo)) < tol):
                    y = Circle(x, -1, 0, 0)
                elif (np.abs(Circle(x, -1, 0, 0) - Circle(x, -ro, xo, yo)) < tol):
                    y = Circle(x, -1, 0, 0)
                else:
                    continue
                ccroots.append((x, y))
    return ccroots

def Analytic(xo, yo, ro):
    '''
    Analytic flux during an occultation assuming a uniform body.
    
    '''
    
    d = np.sqrt(xo ** 2 + yo ** 2)
    if (d < (1 + ro)):
        A = (-d + 1 - ro) * (-d - 1 + ro) * (-d + 1 + ro) * (d + 1 + ro)
    else:
        F = 0.
        return np.pi
    if d + ro <= 1:
        F = np.pi * ro ** 2
    elif d + 1 <= ro:
        F = np.pi
    else:
        F = ro ** 2 * np.arccos((d ** 2 + ro ** 2 - 1) / (2 * d * ro)) \
           + np.arccos((d ** 2 + 1 - ro ** 2) / (2 * d)) \
           - 0.5 * np.sqrt(A)
    
    return np.pi - F
    
def FUnif(x0, y0, r, phi1, phi2):
    '''
    
    '''
    
    F0 = 0.5 * r * (r * (phi2 - phi1)
                  + x0 * (np.sin(phi2) - np.sin(phi1))
                  + y0 * (np.cos(phi1) - np.cos(phi2)))
    
    return F0

def Greens(F, xo, yo, ro, roots):
    '''
    
    '''
    
    # Check for special cases
    if len(roots) != 2:
        d = np.sqrt(xo ** 2 + yo ** 2) 
        if d > (1 + ro):
            # No occultation
            return F(0, 0, 1, 0, 2 * np.pi)
        else:
            # Occultor completely within occulted planet disk
            F0 = F(0, 0, 1, 0, 2 * np.pi)
            F1 = F(xo, yo, ro, 0, 2 * np.pi)
            return F0 - F1

    # Get the roots of intersection
    x1, y1 = roots[0]
    x2, y2 = roots[1]
    
    # First integrate over the boundary of the occulted body.
    # Get the angular position of the roots relative to the origin.
    phi1 = np.arctan2(y1, x1)
    phi2 = np.arctan2(y2, x2)
    
    # We integrate counter-clockwise along the arc that is *outside*
    # the occultor. The starting point is the one that is to the *right*
    # of the occultor. A little trig shows that we must check whether:
    if (-x1 * yo + y1 * xo) < 0:
        phi1, phi2 = phi2, phi1
    if phi2 < phi1:
        phi2 += 2 * np.pi

    # Compute this segment of the line integral
    F0 = F(0, 0, 1, phi1, phi2)
    
    # Now integrate over the boundary of the occultor.
    # Get the angular position of the roots relative to the occultor center.
    # Sort them so we integrate in the counter-clockwise direction
    phi1 = np.arctan2(y1 - yo, x1 - xo)
    phi2 = np.arctan2(y2 - yo, x2 - xo)
    
    # We integrate counter-clockwise along the arc that is *inside*
    # the occulted planet. As before:
    if (-x1 * yo + y1 * xo) < 0:
        phi1, phi2 = phi2, phi1
    if phi2 < phi1:
        phi2 += 2 * np.pi

    # Compute this segment of the line integral
    F1 = F(xo, yo, ro, phi1, phi2)
    
    # The total flux is just the difference
    return F0 - F1
    
def Plot(xo, yo, ro, roots):
    '''
    
    '''
    
    # Set up the figure
    fig, ax = pl.subplots(1, figsize = (6, 6))
    ax.set_aspect(1)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)

    # Plot the planet
    x = np.linspace(-1,1,1000)
    ax.plot(x, Circle(x, 1, 0, 0), 'k-', lw = 1)
    ax.plot(x, Circle(x, -1, 0, 0), 'k-', lw = 1)

    # Plot the occultor
    x = np.linspace(xo - ro, xo + ro, 1000)
    ax.plot(x, Circle(x, ro, xo, yo), 'r-', lw = 1)
    ax.plot(x, Circle(x, -ro, xo, yo), 'r-', lw = 1)
    
    # Plot the intersection points
    for xp, yp in roots:
        pl.plot(xp, yp, 'ro')
    
    # Show!
    pl.show()

def Compute(xo, yo, ro):
    '''
    Occulted body at origin, unit radius.
    Occultor at xo, yo, radius ro.
    
    '''
    
    # Get points of intersection
    roots = Intersection(xo, yo, ro)
    
    # Get non-occulted area analytically
    FA = Analytic(xo, yo, ro)
    
    # Get it using Green's theorem
    FG = Greens(FUnif, xo, yo, ro, roots)
    
    print(FA, FG)
    
    # Plot!
    Plot(xo, yo, ro, roots)

Compute(-0.5, 0.25, 0.75)