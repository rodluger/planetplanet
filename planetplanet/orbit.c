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

int Compute(int npts, double *time, PLANET *planet, SETTINGS *settings, ARRAYS *arr){
  /*
      Compute the orbital parameters
  */    
  double omega, per, a, inc, w, ecc, fi, tperi0, t0;
  double dt, tmp;
  int i;
  int iErr = ERR_NONE;

  arr->npts = npts;
  arr->M = malloc(npts*sizeof(double)); 
  arr->E = malloc(npts*sizeof(double)); 
  arr->f = malloc(npts*sizeof(double)); 
  arr->r = malloc(npts*sizeof(double)); 
  arr->x = malloc(npts*sizeof(double)); 
  arr->y = malloc(npts*sizeof(double)); 
  arr->z = malloc(npts*sizeof(double));
  arr->b = malloc(npts*sizeof(double));

  per = planet->per;                                                                  // Orbital period
  if (!(per > 0.)) return ERR_PER;

  inc = planet->inc;                                                                  // Orbital inclination
  if (!((inc >= 0.) && (inc <= PI / 2.))) return ERR_INC;
  
  a = planet->a;
  t0 = planet->t0;
  
  if (isnan(planet->esw) || isnan(planet->ecw)) {                                     // Eccentricity and longitude of pericenter
    if (isnan(planet->ecc)) return ERR_ECC_W;
    if ((planet->ecc != 0) && isnan(planet->w)) 
      return ERR_ECC_W;
    else if (planet->ecc == 0)
      planet->w = 0;
    else if (isnan(planet->w))
      return ERR_ECC_W;                           
    if ((planet->ecc < 0) || (planet->ecc >= 1)) return ERR_ECC_W;
    if ((planet->w < 0) || (planet->w >= 2 * PI)) return ERR_ECC_W;
    w = planet->w;
    ecc = planet->ecc;
  } else {
    w = atan2(planet->esw, planet->ecw);
    ecc = sqrt(planet->esw * planet->esw + planet->ecw * planet->ecw);
    if ((ecc < 0.) || (ecc >= 1.)) return ERR_BAD_ECC;
    planet->ecc = ecc;
    planet->w = w;
  }
  
  // HACK: My definition of omega in the equations below is apparently
  // off by 180 degrees from Laura Kreidberg's in BATMAN. This isn't elegant,
  // but the two models agree now that I added the following line:
  w = w - PI;
  
  fi = (3. * PI / 2.) - w;                                                            // True anomaly at planet center (Shields et al. 2015)
  tperi0 = per * sqrt(1. - ecc * ecc) / (2. * PI) * (ecc * sin(fi) / 
           (1. + ecc * cos(fi)) - 2. / sqrt(1. - ecc * ecc) * 
           atan2(sqrt(1. - ecc * ecc) * tan(fi/2.), 1. + ecc));                       // Time of pericenter passage (Shields et al. 2015)

  /*
  --- ORBITAL SOLUTION ---
  */
  
  for (i = 0; i < npts ; i++) {
       
    arr->M[i] = 2. * PI / per * modulus(time[i] - tperi0 - t0, per);                  // Mean anomaly
    if (settings->kepsolver == MDFAST)
      arr->E[i] = EccentricAnomalyFast(arr->M[i], ecc, settings->keptol, 
                                       settings->maxkepiter);                         // Eccentric anomaly
    else
      arr->E[i] = EccentricAnomaly(arr->M[i], ecc, settings->keptol, 
                                   settings->maxkepiter);
    if (arr->E[i] == -1) return ERR_KEPLER;
    arr->f[i] = TrueAnomaly(arr->E[i], ecc);                                          // True anomaly
    arr->r[i] = a * (1. - ecc * ecc)/(1. + ecc * cos(arr->f[i]));                     // Star-planet separation
    arr->b[i] = arr->r[i] * sqrt(1. - pow(sin(w + arr->f[i]) * sin(inc), 2.));        // Instantaneous impact parameter                                   
    arr->x[i] = arr->r[i] * cos(w + arr->f[i]);                                       // Cartesian sky-projected coordinates
    arr->z[i] = arr->r[i] * sin(w + arr->f[i]);
    if (arr->b[i] * arr->b[i] - arr->x[i] * arr->x[i] < 1.e-15) 
      arr->y[i] = 0.;                                                                 // Prevent numerical errors
    else {
      tmp = modulus(arr->f[i] + w, 2 * PI);                                           // TODO: Verify this modulus
      arr->y[i] = sqrt(arr->b[i] * arr->b[i] - arr->x[i] * arr->x[i]);
      if (!((0 < tmp) && (tmp < PI))) arr->y[i] *= -1;
    }

  }

	return iErr;
}