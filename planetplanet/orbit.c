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

int Kepler(int np, PLANET **planet, SETTINGS settings){
  /*
      Compute the instantaneous x, y, and z positions of all 
      planets with a simple Keplerian solver
  */
  
  double M, E, f, r, b;
  double tmp;
  int iErr = ERR_NONE;
  int t, p;
  
  // Loop over all planets
  for (p = 0; p < np; p++) {
  
    // Loop over the time array
    for (t = 0; t < planet[p]->nt; t++) {
      
      // Mean anomaly
      M = 2. * PI / planet[p]->per * modulus(planet[p]->time[t] - planet[p]->t0, planet[p]->per);                  

      // Eccentric anomaly
      if (settings.kepsolver == MDFAST)
        E = EccentricAnomalyFast(M, planet[p]->ecc, settings.keptol, settings.maxkepiter);
      else
        E = EccentricAnomaly(M, planet[p]->ecc, settings.keptol, settings.maxkepiter);
      if (E == -1) return ERR_KEPLER;
  
      // True anomaly
      f = TrueAnomaly(E, planet[p]->ecc);        
  
      // Star-planet separation                                  
      r = planet[p]->a * (1. - planet[p]->ecc * planet[p]->ecc)/(1. + planet[p]->ecc * cos(f));                     
  
      // Instantaneous impact parameter  
      b = r * sqrt(1. - pow(sin(planet[p]->w - PI + f) * sin(planet[p]->inc), 2.));                                         

      // Cartesian sky-projected coordinates
      planet[p]->x[t] = r * cos(planet[p]->w - PI + f);                                       
      planet[p]->z[t] = r * sin(planet[p]->w - PI + f);
      if (b * b - planet[p]->x[t] * planet[p]->x[t] < 1.e-15) 
        // Prevent numerical errors
        planet[p]->y[t] = 0.;                                                                 
      else {
        tmp = modulus(f + planet[p]->w - PI, 2 * PI);
        planet[p]->y[t] = sqrt(b * b - planet[p]->x[t] * planet[p]->x[t]);
        if (!((0 < tmp) && (tmp < PI))) 
          planet[p]->y[t] *= -1;

      }
  
    }
  
  }
  
	return iErr;
	
}

void heartbeat(struct reb_simulation* r){
  /*
  
  */
  
  // TODO
	
}

int NBody(int np, PLANET **planet, SETTINGS settings) {
  /*
  
  */
  
  int p, t;
  double tmax;
  double vx[np], vy[np], vz[np];
  struct reb_simulation* r = reb_create_simulation();
  
  // Velocities. TODO: Compute analytically
  int stupid = 1;
  if (stupid) {
    double time[6];
    for (t = 0; t < 6; t++)
      time[t] = planet[0]->time[t];
    t = planet[0]->nt;
    for (p = 0; p < np; p++) {
      planet[p]->time[0] = time[0] - 3e-7;
      planet[p]->time[1] = time[0] - 2e-7;
      planet[p]->time[2] = time[0] - 1e-7;
      planet[p]->time[3] = time[0] + 1e-7;
      planet[p]->time[4] = time[0] + 2e-7;
      planet[p]->time[5] = time[0] + 3e-7;
      planet[p]->nt = 6;
    }
    Kepler(np, planet, settings); 
    for (p = 0; p < np; p++) {
      planet[p]->time[0] = time[0];
      planet[p]->time[1] = time[1];
      planet[p]->time[2] = time[2];
      planet[p]->time[3] = time[3];
      planet[p]->time[4] = time[4];
      planet[p]->time[5] = time[5];
      planet[p]->nt = t;
    }
    for (p = 0; p < np; p++) {
      vx[p] = (-(1/60.) * planet[p]->x[0] + (3/20.) * planet[p]->x[1] - (3/4.) * planet[p]->x[2] + 
              (3/4.) * planet[p]->x[3] - (3/20.) * planet[p]->x[4] + (1/60.) * planet[p]->x[5]) / 1e-7;
      vy[p] = (-(1/60.) * planet[p]->y[0] + (3/20.) * planet[p]->y[1] - (3/4.) * planet[p]->y[2] + 
              (3/4.) * planet[p]->y[3] - (3/20.) * planet[p]->y[4] + (1/60.) * planet[p]->y[5]) / 1e-7;
      vz[p] = (-(1/60.) * planet[p]->z[0] + (3/20.) * planet[p]->z[1] - (3/4.) * planet[p]->z[2] + 
              (3/4.) * planet[p]->z[3] - (3/20.) * planet[p]->z[4] + (1/60.) * planet[p]->z[5]) / 1e-7;
    }
  }
  
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

	// Initialize the planets
	for (p = 0; p < np; p++) {
		struct reb_particle body = {0};
		body.x  = planet[p]->x[0]; 		
		body.y  = planet[p]->y[0];	 	
		body.z  = planet[p]->z[0];
		body.vx = vx[p]; 		
		body.vy = vy[p];	 	
		body.vz = vz[p];
		body.m  = planet[p]->m;
		reb_add(r, body);
	}
	
	// Move to center of mass frame and integrate!
	reb_move_to_com(r);
	reb_integrate(r, tmax);
  
  return 0;
}