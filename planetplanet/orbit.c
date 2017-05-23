#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "rebound.h"
 
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

int Kepler(int np, BODY **body, SETTINGS settings){
  /*
      Compute the instantaneous x, y, and z positions of all 
      planets with a simple Keplerian solver
  */
  
  double M, E, f, r;
  double  cwf, swf, co, so, ci, si;
  int iErr = ERR_NONE;
  int t, p;
  
  // In this simple solver, the star is fixed at the origin (massless planets)
  for (t = 0; t < body[0]->nt; t++) {
    body[0]->x[t] = 0;
    body[0]->y[t] = 0;
    body[0]->z[t] = 0;
  }
  
  // Loop over all planets
  for (p = 1; p < np; p++) {
  
    // Loop over the time array
    for (t = 0; t < body[p]->nt; t++) {
      
      // Mean anomaly
      M = 2. * PI / body[p]->per * modulus(body[p]->time[t] - body[p]->t0, body[p]->per);                  

      // Eccentric anomaly
      if (settings.kepsolver == MDFAST)
        E = EccentricAnomalyFast(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
      else
        E = EccentricAnomaly(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
      if (E == -1) return ERR_KEPLER;
  
      // True anomaly
      f = TrueAnomaly(E, body[p]->ecc);        
  
      // Orbital radius                                  
      r = body[p]->a * (1. - body[p]->ecc * body[p]->ecc)/(1. + body[p]->ecc * cos(f));                     
      
      /*
      // OLD PYSYZYGY CODE
      double b, tmp;
      
      // Instantaneous impact parameter  
      b = r * sqrt(1. - pow(sin(body[p]->w - PI + f) * sin(body[p]->inc), 2.));                                         

      // Cartesian sky-projected coordinates
      body[p]->x[t] = r * cos(body[p]->w - PI + f);                                       
      body[p]->z[t] = r * sin(body[p]->w - PI + f);
      if (b * b - body[p]->x[t] * body[p]->x[t] < 1.e-15) 
        // Prevent numerical errors
        body[p]->y[t] = 0.;                                                                 
      else {
        tmp = modulus(f + body[p]->w - PI, 2 * PI);
        body[p]->y[t] = sqrt(b * b - body[p]->x[t] * body[p]->x[t]);
        if (!((0 < tmp) && (tmp < PI))) 
          body[p]->y[t] *= -1;
      }
      */
      
      // Murray and Dermott p. 51
      cwf = cos(body[p]->w + f);
      swf = sin(body[p]->w + f);    // Negative in pysyzygy
      co = cos(body[p]->Omega);     // Negative in pysyzygy
      so = sin(body[p]->Omega);
      ci = cos(body[p]->inc);
      si = sin(body[p]->inc);
      body[p]->x[t] = r * (co * cwf - so * swf * ci);
      body[p]->y[t] = r * (so * cwf + co * swf * ci);
      body[p]->z[t] = r * swf * si;
      
    }
  
  }
  
	return iErr;
	
}

void heartbeat(struct reb_simulation* r){
  /*
  
  */
  
  // TODO
	
}

int NBody(int np, BODY **body, SETTINGS settings) {
  /*
  
  */
  
  int p, t;
  double tmax;
  double M, E, f;
  struct reb_simulation* r = reb_create_simulation();
  
	// Timestep in days
	r->dt = 0.01;
	
	// Max time in days
	tmax = 10.;
	
	// G in REARTH^3 / MEARTH / day^2
	r->G = 11466.9811868;
	
	// 11th order symplectic corrector
	r->ri_whfast.safe_mode 	= 0;
	r->ri_whfast.corrector 	= 11;		
	r->integrator = REB_INTEGRATOR_WHFAST;
	r->heartbeat = heartbeat;
	r->exact_finish_time = 1;

  // Initialize the star
  struct reb_particle primary = {0};
  primary.m = body[0]->m;
  t = 0;
  
	// Initialize the planets
	for (p = 1; p < np; p++) {
	
	  // Get the true anomaly at the first timestep
    M = 2. * PI / body[p]->per * modulus(body[p]->time[0] - body[p]->t0, body[p]->per);                  
    if (settings.kepsolver == MDFAST)
      E = EccentricAnomalyFast(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
    else
      E = EccentricAnomaly(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
    if (E == -1) return ERR_KEPLER;
    f = TrueAnomaly(E, body[p]->ecc);  
	  
	  
	  struct reb_particle planet = reb_tools_orbit_to_particle(r->G, primary, body[p]->m, 
	                               body[p]->a, body[p]->ecc, body[p]->inc, 
	                               body[p]->Omega, body[p]->w, f);
		reb_add(r, planet);
	}
	
	// Move to center of mass frame and integrate!
	reb_move_to_com(r);
	reb_integrate(r, tmax);
  
  return 0;
}