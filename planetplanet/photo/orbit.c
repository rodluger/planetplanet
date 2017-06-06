#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "rebound.h"
#include "progress.h"

void on_progress(progress_data_t *data) {
  /*
  
  */
  
  progress_write(data->holder);
  
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

int Kepler(int np, BODY **body, SETTINGS settings){
  /*
      Compute the instantaneous x, y, and z positions of all 
      planets with a simple Keplerian solver
  */
  
  double M, E, f, r;
  double  cwf, swf, co, so, ci, si;
  int iErr = ERR_NONE;
  int t, p;
  progress_t *progress = progress_new(body[0]->nt, 50);
  
  if (!settings.quiet) {
    progress->fmt = "[:bar] :percent :elapsed";
    progress_on(progress, PROGRESS_EVENT_PROGRESS, on_progress);
    printf("Computing orbits with the Kepler solver...\n");
  }
  
  // In this simple solver, the star is fixed at the origin (massless planets)
  for (t = 0; t < body[0]->nt; t++) {
    body[0]->x[t] = 0;
    body[0]->y[t] = 0;
    body[0]->z[t] = 0;
  }
  
  // Loop over the time array
  for (t = 0; t < body[0]->nt; t++) {
  
    // Loop over all planets
    for (p = 1; p < np; p++) {
      
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
            
      // Murray and Dermott p. 51
      cwf = cos(body[p]->w + f);
      swf = sin(body[p]->w + f);    // NOTE: This is negative in pysyzygy
      co = cos(body[p]->Omega);     // NOTE: This is negative in pysyzygy
      so = sin(body[p]->Omega);
      ci = cos(body[p]->inc);
      si = sin(body[p]->inc);
      body[p]->x[t] = r * (co * cwf - so * swf * ci);
      body[p]->y[t] = r * (so * cwf + co * swf * ci);
      body[p]->z[t] = r * swf * si;
      
    }
  
    // Display the progress
    if (!settings.quiet) {
	    if (t % 1000 == 0) 
	      progress_tick(progress, 1000);
    }
    
  }
  
  if (!settings.quiet)
    printf("\n");
  
  progress_free(progress);
	return iErr;
	
}

void heartbeat(struct reb_simulation* r){
  /*
  
  */
  
  // Nothing!
	
}

int NBody(int np, BODY **body, SETTINGS settings) {
  /*
  
  */
  
  int p, t;
  double M, E, f;
  struct reb_simulation* r = reb_create_simulation();
  progress_t *progress = progress_new(body[0]->nt, 50);
  
  if (!settings.quiet) {
    progress->fmt = "[:bar] :percent :elapsed";
    progress_on(progress, PROGRESS_EVENT_PROGRESS, on_progress);
    printf("Computing orbits with REBOUND...\n");
  }
  
	// Set the timestep
	r->dt = settings.dt;

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
  reb_add(r, primary);

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
	  
	  // Create the particle in REBOUND
	  struct reb_particle planet = reb_tools_orbit_to_particle(r->G, primary, body[p]->m,
	                               body[p]->a, body[p]->ecc, body[p]->inc, 
	                               body[p]->Omega, body[p]->w, f);
		reb_add(r, planet);
	 
	}
	
	// Move to center of mass frame
	reb_move_to_com(r);
	
	// Store the initial positions
	for (p = 0; p < np; p++) {
	  body[p]->x[0] = r->particles[p].x;
	  body[p]->y[0] = r->particles[p].y;
	  body[p]->z[0] = r->particles[p].z;
	}
	
	// Integrate!
	for (t = 1; t < body[0]->nt; t++) {
    
    // Take one step
	  reb_integrate(r, body[0]->time[t] - body[0]->time[0]);
	  reb_integrator_synchronize(r);
	  
	  // Update body positions
	  for (p = 0; p < np; p++) {
	    body[p]->x[t] = r->particles[p].x;
	    body[p]->y[t] = r->particles[p].y;
	    body[p]->z[t] = r->particles[p].z;
	  }
	  
	  // Display the progress
	  if (!settings.quiet) {
	    if (t % 1000 == 0) 
	      progress_tick(progress, 1000);
    }

  }
  
  if (!settings.quiet)
    printf("\n");
  
  progress_free(progress);
  return 0;
  
}