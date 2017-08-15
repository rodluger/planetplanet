/**
@file orbit.c
@brief Orbital evolution routines.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "rebound.h"
#include "progress.h"

/**
Update the progress bar. Used internally.

*/
void on_progress(progress_data_t *data) {  
  progress_write(data->holder);
} 

/**
The arithmetic modulus, x mod y.

@param x
@param y
@return x mod y
*/
double modulus(double x, double y) {
  return x - y * floor(x / y);
} 
 
/**
Returns the sign of x.

@param x
@return sign(x)
*/ 
double sgn(double x) {
  return (x > 0) - (x < 0);
}

/**
The true anomaly as a function of the eccentric anomaly \a E
and the eccentricity \a ecc.

@param E The eccentric anomaly in radians
@param ecc The eccentricity
@return The true anomaly in radians

*/
double TrueAnomaly(double E, double ecc) {
  if (ecc == 0.) return E;
  else return 2. * atan2(pow(1. + ecc, 0.5) * sin(E / 2.), 
                         pow(1. - ecc, 0.5) * cos(E / 2.));
}

/**
Computes the eccentric anomaly from the mean anomaly and the eccentricity
using a fast iterative scheme. Adapted from Rory Barnes, based on equations 
in Murray & Dermott.

@param dMeanA The mean anomaly in radians
@param dEcc The eccentricity
@param tol The solver tolerance
@param maxiter The maximum number of iterations
@return The eccentric anomaly in radians

*/
double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter) {

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

/**  
Computes the eccentric anomaly from the mean anomaly and the eccentricity.
This is a simpler version of the Kepler solver, borrowed from
https://github.com/lkreidberg/batman/blob/master/c_src/_rsky.c

@param M The mean anomaly in radians
@param e The eccentricity
@param tol The solver tolerance
@param maxiter The maximum number of iterations
@return The eccentric anomaly in radians
*/
double EccentricAnomaly(double M, double e, double tol, int maxiter) {
  
  double E = M, eps = tol;                                                            // Kreidberg: eps = 1.0e-7;
  int iter;
  
  if (e == 0.) return M;                                                              // The trivial circular case
  
  for (iter = 1; iter <= maxiter; iter++) {
    E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
    if (fabs(E - e*sin(E) - M) <= eps) break;
  }
  
  if (iter >= maxiter) 
    return -1.;                                                                       // Solver didn't converge
  else
	  return E;
	
}

/**
Compute the instantaneous x, y, and z positions of all 
planets with a simple multi-Keplerian solver.

@param np The number of particles (bodies) in the system
@param body An array of BODY pointers to each of the bodies in the system
@param settings An instance of the SETTINGS class
@return The error code

*/
int Kepler(int np, BODY **body, SETTINGS settings){
  
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
      M = 2. * PI / body[p]->per * modulus(body[p]->time[t] - body[p]->tperi0, body[p]->per);                  

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

/**
Called at each step of the N-Body simulation. Currently does nothing!

*/
void heartbeat(struct reb_simulation* r){
  // Nothing!
}

/**
Compute the instantaneous x, y, and z positions of all 
planets using the REBOUND N-Body code.

@param np The number of particles (bodies) in the system
@param body An array of BODY pointers to each of the bodies in the system
@param settings An instance of the SETTINGS class
@return The error code

*/
int NBody(int np, BODY **body, SETTINGS settings) {
  /*
  
  */
  
  int p, t, i;
  double M, E, f, d;
  double cwf, swf, co, so, ci, si;
  int last_t;
  int moons = 0;
  struct reb_simulation* r = reb_create_simulation();
  progress_t *progress = progress_new(body[0]->nt, 50);
  
  if (!settings.quiet) {
    progress->fmt = "[:bar] :percent :elapsed";
    progress_on(progress, PROGRESS_EVENT_PROGRESS, on_progress);
    printf("Computing orbits with REBOUND...\n");
  }
  
	// Set the timestep
	r->dt = settings.timestep;

	// G in REARTH^3 / MEARTH / day^2
	r->G = GEARTH;
	
	// Settings for WHFAST: 11th order symplectic corrector
	r->ri_whfast.safe_mode 	= 0;
	r->ri_whfast.corrector 	= 11;		
	r->integrator = settings.integrator;
	r->heartbeat = heartbeat;
	r->exact_finish_time = 1;

  // Initialize the star
  struct reb_particle primary = {0};
  primary.m = body[0]->m;
  reb_add(r, primary);

	// Initialize the planets
	for (p = 1; p < np; p++) {
	  
	  // Get the true anomaly at the first timestep
    M = 2. * PI / body[p]->per * modulus(body[p]->time[0] - body[p]->tperi0, body[p]->per);                  
    if (settings.kepsolver == MDFAST)
      E = EccentricAnomalyFast(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
    else
      E = EccentricAnomaly(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
    if (E == -1) return ERR_KEPLER;
    f = TrueAnomaly(E, body[p]->ecc);  
	  
	  // Create the particle in REBOUND
    struct reb_particle planet = reb_tools_orbit_to_particle(r->G, r->particles[body[p]->host], body[p]->m,
                                 body[p]->a, body[p]->ecc, body[p]->inc, 
                                 body[p]->Omega, body[p]->w, f);
    reb_add(r, planet);
    
    // Is this a moon?
    if (body[p]->host != 0) moons++;
     
	}
	
	// Move to center of mass frame
	reb_move_to_com(r);
	
	// Store the initial positions
	last_t = 0;
	for (p = 0; p < np; p++) {
	  body[p]->x[0] = r->particles[p].x;
	  body[p]->y[0] = r->particles[p].y;
	  body[p]->z[0] = r->particles[p].z;
	}
	
	// Integrate!
	for (t = 1; t < body[0]->nt; t++) {
    
    // Do we need to synchronize this time step?
    if ((body[0]->time[t] - body[0]->time[last_t] >= settings.timestep) || (t == body[0]->nt - 1) || moons) {
    
      // Yes: take one step
      reb_integrate(r, body[0]->time[t] - body[0]->time[0]);
      reb_integrator_synchronize(r);
    
      // Update body positions
      for (p = 0; p < np; p++) {
        body[p]->x[t] = r->particles[p].x;
        body[p]->y[t] = r->particles[p].y;
        body[p]->z[t] = r->particles[p].z;
        
        // Go back and compute the ones we skipped
        // with a Keplerian solver using the current
        // osculating elements
        // NOTE: I am unable to assign a planet as the primary here! Leads to terrible orbital errors. Need to investigate why.
        struct reb_orbit orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], primary);

        for (i = last_t + 1; i < t; i++) {
          
          // Is this the star? If so, assume it hasn't moved
          if (p == 0) {
            
            // Sky coordinates
            body[p]->x[i] = body[p]->x[t];
            body[p]->y[i] = body[p]->y[t];
            body[p]->z[i] = body[p]->z[t];
            
          } else {
          
            // Mean anomaly
            M = modulus(orbit.M + orbit.n * (body[0]->time[i] - body[0]->time[t]), 2 * PI);
          
            // Eccentric anomaly
            if (settings.kepsolver == MDFAST)
              E = EccentricAnomalyFast(M, orbit.e, settings.keptol, settings.maxkepiter);
            else
              E = EccentricAnomaly(M, orbit.e, settings.keptol, settings.maxkepiter);
            if (E == -1) return ERR_KEPLER;
          
            // True anomaly
            f = TrueAnomaly(E, orbit.e);
          
            // Radial distance                                        
            d = orbit.a * (1. - orbit.e * orbit.e) / (1. + orbit.e * cos(f));
          
            // Murray and Dermott, p. 51
            cwf = cos(orbit.omega + f);
            swf = sin(orbit.omega + f);
            co = cos(orbit.Omega); 
            so = sin(orbit.Omega);
            ci = cos(orbit.inc);
            si = sin(orbit.inc);
          
            // Sky coordinates
            body[p]->x[i] = d * (co * cwf - so * swf * ci);
            body[p]->y[i] = d * (so * cwf + co * swf * ci);
            body[p]->z[i] = d * swf * si;
          
          }
          
        }
        
      }

	    // Reset
	    last_t = t;
	    
	  }
	  
	  // Display the progress
	  if (!settings.quiet) {
	    if (t % 1000 == 0) 
	      progress_tick(progress, 1000);
    }

  }
  
  
  // BUG: The orbital elements for moons are screwed up.
  

  // Update the orbital elements of all the bodies
  for (p = 0; p < np; p++) {
    
    struct reb_orbit orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], r->particles[body[p]->host]);
    body[p]->a = orbit.a;
    body[p]->ecc = orbit.e;
    body[p]->inc = orbit.inc;
    body[p]->w = orbit.omega;
    body[p]->Omega = orbit.Omega;
    body[p]->tperi0 = body[p]->time[body[p]->nt - 1] - body[p]->per * orbit.M / (2 * PI);
    body[p]->per = orbit.P;
    
  } 

    
  if (!settings.quiet)
    printf("\n");
  
  progress_free(progress);
  return 0;
  
}