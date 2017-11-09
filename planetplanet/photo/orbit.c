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

@param x A double
@param y A double
@return x mod y
*/
double modulus(double x, double y) {
  return x - y * floor(x / y);
}

/**
Returns the sign of x.

@param x A double
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
  double cwf, swf, co, so, ci, si, ecw, esw, amp;
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

      // Sky velocity. Based on differentiating the
      // expression on p. 51 of Murray and Dermott
      ecw = body[p]->ecc * cos(body[p]->w);
      esw = body[p]->ecc * sin(body[p]->w);
      amp = (2. * PI / body[p]->per) * body[p]->a / sqrt(1 - body[p]->ecc * body[p]->ecc);
      body[p]->vx[t] = -amp * (co * (esw + swf) - ci * so * (ecw + cwf));
      body[p]->vy[t] = amp * (ci * co * (ecw + cwf) - so * (esw + swf));
      body[p]->vz[t] = amp * si * (ecw + cwf);

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
@param halt_on_occultation Halt when a certain number of occultations of a given body have occurred?
@param occulted The index of the body being occulted (only used if `halt_on_occultation` = 1)
@param noccultors Size of occultor array (only used if `halt_on_occultation` = 1)
@param occultors The indices of the occultors of `occulted` (only used if `halt_on_occultation` = 1)
@param noccultations The number of occultations to look for (only used if `halt_on_occultation` = 1)
@param occultation_times The times of each of the occultations, computed by this function (only used if `halt_on_occultation` = 1)
@param occultation_inds The indices of the occultors corresponding to each of the occultations identified by this function (only used if `halt_on_occultation` = 1)
@param occultation_durs The duration of each of the occultations identified by this function (only used if `halt_on_occultation` = 1)
@return The error code

*/
int NBody(int np, BODY **body, SETTINGS settings, int halt_on_occultation, int occulted, int noccultors, int occultors[noccultors], int noccultations, double occultation_times[noccultations], int occultation_inds[noccultations], double occultation_durs[noccultations]) {

  int p, t, i, o, icenter;
  double M, E, f, d;
  double dx, dy, dz;
  double cwf, swf, co, so, ci, si, ecw, esw, amp;
  int last_t;
  int moons = 0;
  int cartesian = 0;
  int ocount = 0;
  struct reb_simulation* r = reb_create_simulation();
  progress_t *progress = progress_new(body[0]->nt, 50);
  int istart[noccultors];
  if (halt_on_occultation) {
    for (i = 0; i < noccultors; i++)
      istart[i] = -1;
  }

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
  r->ri_whfast.safe_mode = 0;
  r->ri_whfast.corrector = 11;
  r->integrator = settings.integrator;
  r->heartbeat = heartbeat;
  r->exact_finish_time = 1;

  // Initialize the particles
  for (p = 0; p < np; p++) {

      if (!(body[p]->cartesian)) {

        // We're doing orbital elements.
        // Get the true anomaly at the first timestep
        M = 2. * PI / body[p]->per * modulus(body[p]->time[0] - body[p]->tperi0, body[p]->per);
        if (settings.kepsolver == MDFAST)
          E = EccentricAnomalyFast(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
        else
          E = EccentricAnomaly(M, body[p]->ecc, settings.keptol, settings.maxkepiter);
        if (E == -1) return ERR_KEPLER;
        f = TrueAnomaly(E, body[p]->ecc);

        // Instantiate it w/ respect to either the center of mass or its host body
        struct reb_particle refbody;
        if (body[p]->host == -1) refbody = reb_get_com(r);
        else refbody = r->particles[body[p]->host];
        struct reb_particle particle = reb_tools_orbit_to_particle(r->G, refbody, body[p]->m,
                                       body[p]->a, body[p]->ecc, body[p]->inc,
                                       body[p]->Omega, body[p]->w, f);
        reb_add(r, particle);

      } else {

        // We're doing Cartesian coordinates.
        struct reb_particle particle = {0};
        particle.x = body[p]->x[0];
        particle.y = body[p]->y[0];
        particle.z = body[p]->z[0];
        particle.vx = body[p]->vx[0];
        particle.vy = body[p]->vy[0];
        particle.vz = body[p]->vz[0];
        particle.m = body[p]->m;
        reb_add(r, particle);
        cartesian = 1;

      }

      // Is this a moon?
      if (body[p]->host > 0) moons++;

  }

  // Move to center of mass frame
  reb_move_to_com(r);

  // Store the initial positions
  last_t = 0;
  for (p = 0; p < np; p++) {
      body[p]->x[0] = r->particles[p].x;
      body[p]->y[0] = r->particles[p].y;
      body[p]->z[0] = r->particles[p].z;
      body[p]->vx[0] = r->particles[p].vx;
      body[p]->vy[0] = r->particles[p].vy;
      body[p]->vz[0] = r->particles[p].vz;
  }

  // Integrate!
  for (t = 1; t < body[0]->nt; t++) {

    // Do we need to synchronize this time step?
    if ((body[0]->time[t] - body[0]->time[last_t] >= settings.timestep) || (t == body[0]->nt - 1) || moons || cartesian || (settings.nstars > 1)) {

      // Yes: take one step
      reb_integrate(r, body[0]->time[t] - body[0]->time[0]);
      reb_integrator_synchronize(r);

      // Update body positions
      for (p = 0; p < np; p++) {
        body[p]->x[t] = r->particles[p].x;
        body[p]->y[t] = r->particles[p].y;
        body[p]->z[t] = r->particles[p].z;
        body[p]->vx[t] = r->particles[p].vx;
        body[p]->vy[t] = r->particles[p].vy;
        body[p]->vz[t] = r->particles[p].vz;

        // Go back and compute the ones we skipped
        // with a Keplerian solver using the current
        // osculating elements.
        struct reb_particle refbody;
        if (body[p]->host == -1) refbody = reb_get_com(r);
        else refbody = r->particles[body[p]->host];
        struct reb_orbit orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], refbody);

        // TODO: Currently the Keplerian updates below are disabled
        // when there are exomoons, binary star systems, or any sort
        // of nested problem. In principle this could work, but I need
        // to add the position and velocity of the host to the position
        // and velocity vectors of the secondary body after computing
        // its evolution along the osculating orbit.
        for (i = last_t + 1; i < t; i++) {

          // Is this the star? If so, assume it hasn't moved
          if (p == 0) {

            // Sky coordinates
            body[p]->x[i] = body[p]->x[t];
            body[p]->y[i] = body[p]->y[t];
            body[p]->z[i] = body[p]->z[t];
            body[p]->vx[i] = body[p]->vx[t];
            body[p]->vy[i] = body[p]->vy[t];
            body[p]->vz[i] = body[p]->vz[t];

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

            // Sky velocity. Based on differentiating the
            // expression on p. 51 of Murray and Dermott
            ecw = orbit.e * cos(orbit.omega);
            esw = orbit.e * sin(orbit.omega);
            amp = orbit.n * orbit.a / sqrt(1 - orbit.e * orbit.e);
            body[p]->vx[i] = -amp * (co * (esw + swf) - ci * so * (ecw + cwf));
            body[p]->vy[i] = amp * (ci * co * (ecw + cwf) - so * (esw + swf));
            body[p]->vz[i] = amp * si * (ecw + cwf);

          }

        }

      }

      // SPECIAL: Are we checking for occultations?
      if (halt_on_occultation) {

        // Loop over the last set of Kepler steps
        for (i = last_t + 1; i <= t; i++) {

          // Loop over all possible occultors
          for (o = 0; o < noccultors; o++) {

            // Skip self
            if (occultors[o] == occulted) continue;

            // Compute the body-body separation between the two bodies
            dx = body[occultors[o]]->x[i] - body[occulted]->x[i];
            dy = body[occultors[o]]->y[i] - body[occulted]->y[i];
            dz = body[occultors[o]]->z[i] - body[occulted]->z[i];
            d = sqrt(dx * dx + dy * dy);

            // Is an occultation occurring?
            if ((d <= body[occulted]->r + body[occultors[o]]->r) && (dz < 0)){

              // YES. Is this the first index of the occultation?
              if (istart[o] < 0) {

                // YES. Let's keep track of it
                istart[o] = i;

              }

            } else if (istart[o] > 0) {

              // NO. But one just ended!
              // The index of the center of the occultation is...
              icenter = (i + istart[o]) / 2;

              // Log the time of the occultation and the occultor index
              occultation_times[ocount] = body[0]->time[icenter];
              occultation_inds[ocount] = occultors[o];
              occultation_durs[ocount++] = body[0]->time[i] - body[0]->time[istart[o]];

              // Check if we've already found `noccultations` occultations
              if (ocount == noccultations) {

                // Return!
                if (!settings.quiet)
                  printf("\nFound %d occultations. Returning...\n", noccultations);
                progress_free(progress);
                return 0;

              }

              // Reset the occultation flag for this occultor
              istart[o] = -1;

            }

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

  // TODO: BUG: The orbital elements for moons
  // and multiple star systems are screwed up. This
  // does not affect the photodynamical model! It is
  // only a problem if the user is interested in
  // the updated orbital parameters (period, semi-major
  // axis, etc) of the moon after running the N-body code.

  // Move to center of mass frame
  reb_move_to_com(r);

  for (p = 0; p < np; p++) {
    if (body[p]->host == -1) {
      struct reb_particle cm = reb_get_com(r);
      struct reb_orbit orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], cm);
      body[p]->a = orbit.a;
      body[p]->ecc = orbit.e;
      body[p]->inc = orbit.inc;
      body[p]->w = orbit.omega;
      body[p]->Omega = orbit.Omega;
      body[p]->tperi0 = body[p]->time[body[p]->nt - 1] - body[p]->per * orbit.M / (2 * PI);
      body[p]->per = orbit.P;
    } else {
      struct reb_orbit orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], r->particles[body[p]->host]);
      body[p]->a = orbit.a;
      body[p]->ecc = orbit.e;
      body[p]->inc = orbit.inc;
      body[p]->w = orbit.omega;
      body[p]->Omega = orbit.Omega;
      body[p]->tperi0 = body[p]->time[body[p]->nt - 1] - body[p]->per * orbit.M / (2 * PI);
      body[p]->per = orbit.P;
    }
  }

  // Log
  if (!settings.quiet) {
    printf("\n");
    if (halt_on_occultation) {
      if (ocount > 0)
        printf("Requested %d occultations, but only %d found.\n", noccultations, ocount);
      else
        printf("No occultations found over the specified time interval.\n");
    }
  }

  // Return!
  progress_free(progress);
  if (halt_on_occultation)
    return ERR_TOO_FEW_OCCS;
  else
    return 0;

}
