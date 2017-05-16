#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "orbit.h"
 
void dbl_free(double *ptr){
  /* 
      Called by python to free a double pointer
  */ 
  free(ptr);
} 
 
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

int OrbitXYZ(double time, PLANET *planet, SETTINGS *settings){
  /*
      Compute the orbital parameters
  */
  
  double per, a, inc, w, ecc, fi, tperi0, t0;
  double M, E, f, r, b;
  double tmp;
  int iErr = ERR_NONE;

  // Orbital period
  per = planet->per;                                                                  
  if (!(per > 0.)) return ERR_INPUT;
  
  // Orbital inclination
  inc = planet->inc;                                                                  
  if (!((inc >= 0.) && (inc <= PI / 2.))) return ERR_INPUT;
  
  // Semi-major axis
  a = planet->a;
  if (!(a > 0.)) return ERR_INPUT;
  
  // Time of transit/inferior conjunction
  t0 = planet->t0;
  
  // Longitude of pericenter
  w = planet->w;
  if (!((w >= 0.) && (w <= 2 * PI))) return ERR_INPUT;
  
  // NOTE: Adding a phase shift to be consistent with BATMAN
  w = w - PI;
  
  // Eccentricity
  ecc = planet->ecc;
  if (!((ecc >= 0.) && (ecc < 1))) return ERR_INPUT;
  
  // True anomaly at planet center (Shields et al. 2015)
  fi = (3. * PI / 2.) - w;                
  
  // Time of pericenter passage (Shields et al. 2015)                                            
  tperi0 = per * sqrt(1. - ecc * ecc) / (2. * PI) * (ecc * sin(fi) / 
           (1. + ecc * cos(fi)) - 2. / sqrt(1. - ecc * ecc) * 
           atan2(sqrt(1. - ecc * ecc) * tan(fi/2.), 1. + ecc));                       

  // Mean anomaly
  M = 2. * PI / per * modulus(time - tperi0 - t0, per);                  
  
  // Eccentric anomaly
  if (settings->kepsolver == MDFAST)
    E = EccentricAnomalyFast(M, ecc, settings->keptol, settings->maxkepiter);
  else
    E = EccentricAnomaly(M, ecc, settings->keptol, settings->maxkepiter);
  if (E == -1) return ERR_KEPLER;
  
  // True anomaly
  f = TrueAnomaly(E, ecc);        
  
  // Star-planet separation                                  
  r = a * (1. - ecc * ecc)/(1. + ecc * cos(f));                     
  
  // Instantaneous impact parameter  
  b = r * sqrt(1. - pow(sin(w + f) * sin(inc), 2.));                                         
  
  // Cartesian sky-projected coordinates
  planet->x = r * cos(w + f);                                       
  planet->z = r * sin(w + f);
  if (b * b - planet->x * planet->x < 1.e-15) 
    // Prevent numerical errors
    planet->y = 0.;                                                                 
  else {
    tmp = modulus(f + w, 2 * PI);
    planet->y = sqrt(b * b - planet->x * planet->x);
    if (!((0 < tmp) && (tmp < PI))) 
      planet->y *= -1;
  }

	return iErr;
	
}