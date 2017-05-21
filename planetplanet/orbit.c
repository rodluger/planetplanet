#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
 
double modulus(double x, double y) {
  /*
      The arithmetic modulus, x mod y
  */
  return x - y * floor(x / y);
} 
 
double sgn(double x) {
  /* 
      Returns the sign of x
  */ 
  return (x > 0) - (x < 0);
}
 
double TrueAnomaly(double E, double ecc) {
  /* 
      The true anomaly as a function of the eccentric anomaly E
      and the eccentricity ecc
  */ 
  if (ecc == 0.) return E;
  else return 2. * atan2(pow(1. + ecc, 0.5) * sin(E / 2.), 
                         pow(1. - ecc, 0.5) * cos(E / 2.));
}

double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter) {
  /* 
      Adapted from Rory Barnes, based on Murray & Dermott 
  */
  
  double dEccA;
  double di_1, di_2, di_3 = 1.0, fi, fi_1, fi_2, fi_3;
  double lo = -2 * PI;
  double up = 2 * PI;
  double next;
  int iter;
  
  if (dEcc == 0.) return dMeanA;                                                      // The trivial circular case
  dEccA = dMeanA + sgn(sin(dMeanA))*0.85*dEcc;

  for (iter = 1; iter <= maxiter; iter++) {
    fi = dEccA - dEcc*sin(dEccA) - dMeanA;
    if (fi > 0)
      up = dEccA;
    else
      lo = dEccA;
    fi_1 = 1.0 - dEcc*cos(dEccA);
    fi_2 = dEcc*sin(dEccA);
    fi_3 = dEcc*cos(dEccA);
    di_1 = -fi / fi_1;
    di_2 = -fi / (fi_1 + 0.5*di_1*fi_2);
    di_3 = -fi / (fi_1 + 0.5*di_2*fi_2 + 1./6.*di_2*di_2*fi_3);
    next = dEccA + di_3;
    
    if (fabs(dEccA - next) < tol) 
      break;
      
    if ((next > lo) && (next < up)) 
      dEccA = next;
    else 
      dEccA = (lo + up) / 2.;
      
    if ((fabs(dEccA - lo) < tol) || (fabs(dEccA - up) < tol))
      break;
  }
  
  if (iter >= maxiter) 
    return -1.;                                                                       // Solver didn't converge
  else
    return dEccA;
}

double EccentricAnomaly(double M, double e, double tol, int maxiter) {
  /*  
      A simpler version of the Kepler solver, borrowed from
      https://github.com/lkreidberg/batman/blob/master/c_src/_rsky.c
  */
  
  double E = M, eps = tol;                                                            // Kreidberg: eps = 1.0e-7;
  
  if (e == 0.) return M;                                                              // The trivial circular case
  
	while(fabs(E - e*sin(E) - M) > eps) E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	return E;
	
}

int OrbitXYZ(int np, PLANET **planet, SETTINGS settings){
  /*
      Compute the orbital parameters
  */
  
  double M, E, f, r, b;
  double tmp;
  int iErr = ERR_NONE;
  int p;
  
  // Loop over all planets
  for (p = 0; p < np; p++) {
  
    // Mean anomaly
    M = 2. * PI / planet[p]->per * modulus(planet[p]->time[planet[p]->t] - planet[p]->t0, planet[p]->per);                  

    // Eccentric anomaly
    if (settings.kepsolver == MDFAST)
      E = EccentricAnomalyFast(M, planet[p]->ecc, settings.keptol, settings.maxkepiter);
    else
      E = EccentricAnomaly(M, planet[p]->ecc, settings.keptol, settings.maxkepiter);
    if (E == -1) return ERR_KEPLER;
  
    // True anomaly
    f = TrueAnomaly(E, planet[p]->ecc);        
  
    // Star-planet[p] separation                                  
    r = planet[p]->a * (1. - planet[p]->ecc * planet[p]->ecc)/(1. + planet[p]->ecc * cos(f));                     
  
    // Instantaneous impact parameter  
    b = r * sqrt(1. - pow(sin(planet[p]->w - PI + f) * sin(planet[p]->inc), 2.));                                         

    // Cartesian sky-projected coordinates
    planet[p]->x[planet[p]->t] = r * cos(planet[p]->w - PI + f);                                       
    planet[p]->z[planet[p]->t] = r * sin(planet[p]->w - PI + f);
    if (b * b - planet[p]->x[planet[p]->t] * planet[p]->x[planet[p]->t] < 1.e-15) 
      // Prevent numerical errors
      planet[p]->y[planet[p]->t] = 0.;                                                                 
    else {
      tmp = modulus(f + planet[p]->w - PI, 2 * PI);
      planet[p]->y[planet[p]->t] = sqrt(b * b - planet[p]->x[planet[p]->t] * planet[p]->x[planet[p]->t]);
      if (!((0 < tmp) && (tmp < PI))) 
        planet[p]->y[planet[p]->t] *= -1;

    }
  
  }
  
	return iErr;
	
}