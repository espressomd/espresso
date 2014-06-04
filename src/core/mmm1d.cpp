/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file mmm1d.cpp  MMM1D algorithm for long range coulomb interaction.
 *
 *  For more information about MMM1D, see \ref mmm1d.hpp "mmm1d.h".
 */

#include "utils.hpp"
#include "mmm1d.hpp"
#include "polynom.hpp"
#include "specfunc.hpp"
#include "communication.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "tuning.hpp"
#include "interaction_data.hpp"
#include "mmm-common.hpp"
#include "errorhandling.hpp"

#ifdef ELECTROSTATICS

/** How many trial calculations */
#define TEST_INTEGRATIONS 1000

/** Largest numerically stable cutoff for Bessel function. Don't
    change without improving the formulas. */
#define MAXIMAL_B_CUT 30

/** minimal radius for the far formula in multiples of box_l[2] */
#define MIN_RAD 0.01

/** if you define this, the Besselfunctions are calculated up
    to machine precision, otherwise 10^-14, which should be
    definitely enough for daily life. */
#undef BESSEL_MACHINE_PREC

#ifndef BESSEL_MACHINE_PREC
#define K0 LPK0
#define K1 LPK1
#endif

/** inverse box dimensions and other constants */
/*@{*/
static double uz, L2, uz2, prefuz2, prefL3_i;
/*@}*/

MMM1D_struct mmm1d_params = { 0.05, 1e-5 };
/** From which distance a certain bessel cutoff is valid. Can't be part of the
    params since these get broadcasted. */
static double *bessel_radii;
  

static double far_error(int P, double minrad)
{
  // this uses an upper bound to all force components and the potential
  double rhores = 2*M_PI*uz*minrad;
  double pref = 4*uz*dmax(1, 2*M_PI*uz);

  return pref*K1(rhores*P)*exp(rhores)/rhores*(P - 1 + 1/rhores);
}

static double determine_minrad(double maxPWerror, int P)
{
  // bisection to search for where the error is maxPWerror
  double rgranularity = MIN_RAD*box_l[2];
  double rmin = rgranularity;
  double rmax = dmin(box_l[0], box_l[1]);
  double errmin = far_error(P, rmin);
  double errmax = far_error(P, rmax);
  if (errmin < maxPWerror) {
    // we can do almost all radii with this P
    return rmin;
  }
  if (errmax > maxPWerror) {
    // make sure that this switching radius cannot be reached
    return 2*dmax(box_l[0], box_l[1]);
  }

  while (rmax - rmin > rgranularity) {
    double c = 0.5*(rmin + rmax);
    double errc = far_error(P, c);
    if (errc > maxPWerror) {
      rmin = c;
    } else {
      rmax = c;
    }
  }
  return 0.5*(rmin + rmax);
}

static void determine_bessel_radii(double maxPWerror, int maxP)
{
  bessel_radii = (double *)realloc(bessel_radii, sizeof(double)*maxP);
  for (int P = 1; P <= maxP; ++P) {
    bessel_radii[P-1] = determine_minrad(maxPWerror, P);
    //printf("cutoff %d %f\n", P, bessel_radii[P-1]);
  }
}

static void prepare_polygamma_series(double maxPWerror, double maxrad2)
{
  /* polygamma, determine order */
  int n;
  double err;
  double rhomax2nm2, rhomax2 = uz2*maxrad2;
  /* rhomax2 < 1, so rhomax2m2 falls monotonously */

  n = 1;
  rhomax2nm2 = 1.0;
  do {
    create_mod_psi_up_to(n+1);

    /* |uz*z| <= 0.5 */
    err = 2*n*fabs(mod_psi_even(n, 0.5))*rhomax2nm2;
    rhomax2nm2 *= rhomax2;
    n++;
    // fprintf(stderr, "%f\n", err);
  }
  while (err > 0.1*maxPWerror);
}

int MMM1D_set_params(double switch_rad, double maxPWerror)
{
  mmm1d_params.far_switch_radius_2 = (switch_rad > 0) ? SQR(switch_rad) : -1;
  mmm1d_params.maxPWerror = maxPWerror;
  coulomb.method = COULOMB_MMM1D;

  mpi_bcast_coulomb_params();

  return 0;
}

int MMM1D_sanity_checks()
{
  char *errtxt;
  if (PERIODIC(0) || PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{022 MMM1D requires periodicity 0 0 1} ");
    return 1;
  }

  if (cell_structure.type != CELL_STRUCTURE_NSQUARE) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{023 MMM1D requires n-square cellsystem} ");
    return 1;
  }
  return 0;
}

void MMM1D_init()
{
  if (MMM1D_sanity_checks()) return;

  if (mmm1d_params.far_switch_radius_2 >= SQR(box_l[2]))
    mmm1d_params.far_switch_radius_2 = 0.8*SQR(box_l[2]);

  uz  = 1/box_l[2];
  L2  = box_l[2]*box_l[2];
  uz2 = uz*uz;
  prefuz2 = coulomb.prefactor*uz2;
  prefL3_i = prefuz2*uz;

  determine_bessel_radii(mmm1d_params.maxPWerror, MAXIMAL_B_CUT);
  prepare_polygamma_series(mmm1d_params.maxPWerror, mmm1d_params.far_switch_radius_2);
}

void add_mmm1d_coulomb_pair_force(double chpref, double d[3], double r2, double r, double force[3])
{
  int dim;
  double F[3];
  double rxy2, rxy2_d, z_d;
  double pref;
  double Fx, Fy, Fz;
  

  rxy2   = d[0]*d[0] + d[1]*d[1];
  rxy2_d = rxy2*uz2;
  z_d    = d[2]*uz;

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* near range formula */
    double sr, sz, r2nm1, rt, rt2, shift_z;
    int n;

    /* polygamma summation */
    sr = 0;
    sz = mod_psi_odd(0, z_d);

    r2nm1 = 1.0;
    for (n = 1; n < n_modPsi; n++) {
      double deriv = 2*n;
      double mpe   = mod_psi_even(n, z_d);
      double mpo   = mod_psi_odd(n, z_d);
      double r2n   = r2nm1*rxy2_d;

      sz +=         r2n*mpo;
      sr += deriv*r2nm1*mpe;

      if (fabs(deriv*r2nm1*mpe) < mmm1d_params.maxPWerror)
	break;

      r2nm1 = r2n;
    }
    // fprintf(stderr, "max_n %d\n", n);

    Fx = prefL3_i*sr*d[0];
    Fy = prefL3_i*sr*d[1];
    Fz = prefuz2*sz;

    /* real space parts */

    pref = coulomb.prefactor/(r2*r); 
    Fx += pref*d[0];
    Fy += pref*d[1];
    Fz += pref*d[2];

    shift_z = d[2] + box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    rt  = sqrt(rt2);
    pref = coulomb.prefactor/(rt2*rt); 
    Fx += pref*d[0];
    Fy += pref*d[1];
    Fz += pref*shift_z;

    shift_z = d[2] - box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    rt  = sqrt(rt2);
    pref = coulomb.prefactor/(rt2*rt); 
    Fx += pref*d[0];
    Fy += pref*d[1];
    Fz += pref*shift_z;

    F[0] = Fx;
    F[1] = Fy;
    F[2] = Fz;
  }
  else {
    /* far range formula */
    double rxy   = sqrt(rxy2);
    double rxy_d = rxy*uz;
    double sr = 0, sz = 0;
    int bp;

    for (bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp-1] < rxy)
        break;

      double fq = C_2PI*bp, k0, k1;
#ifdef BESSEL_MACHINE_PREC
      k0 = K0(fq*rxy_d);
      k1 = K1(fq*rxy_d);
#else
      LPK01(fq*rxy_d, &k0, &k1);
#endif
      sr += bp*k1*cos(fq*z_d);
      sz += bp*k0*sin(fq*z_d);
    }
    sr *= uz2*4*C_2PI;
    sz *= uz2*4*C_2PI;
    
    pref = coulomb.prefactor*(sr/rxy + 2*uz/rxy2);

    F[0] = pref*d[0];
    F[1] = pref*d[1];
    F[2] = coulomb.prefactor*sz;
  }

  for (dim = 0; dim < 3; dim++)
    force[dim] += chpref * F[dim];
}

double mmm1d_coulomb_pair_energy(Particle *p1, Particle *p2, double d[3], double r2, double r)
{
  double chpref = p1->p.q*p2->p.q;
  double rxy2, rxy2_d, z_d;
  double E;

  if (chpref == 0)
    return 0;

  rxy2   = d[0]*d[0] + d[1]*d[1];
  rxy2_d = rxy2*uz2;
  z_d    = d[2]*uz;

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* near range formula */
    double r2n, rt, shift_z;
    int n;

    E = -2*C_GAMMA;

    /* polygamma summation */
    r2n = 1.0;
    for (n = 0; n < n_modPsi; n++) {
      double add = mod_psi_even(n, z_d)*r2n;
      E -= add;
      
      if (fabs(add) < mmm1d_params.maxPWerror)
	break;
 
      r2n *= rxy2_d;
    }
    E *= coulomb.prefactor*uz;

    /* real space parts */

    E += coulomb.prefactor/r; 

    shift_z = d[2] + box_l[2];
    rt = sqrt(rxy2 + shift_z*shift_z);
    E += coulomb.prefactor/rt; 

    shift_z = d[2] - box_l[2];
    rt = sqrt(rxy2 + shift_z*shift_z);
    E += coulomb.prefactor/rt; 
  }
  else {
    /* far range formula */
    double rxy   = sqrt(rxy2);
    double rxy_d = rxy*uz;
    int bp;
    /* The first Bessel term will compensate a little bit the
       log term, so add them close together */
    E = -0.25*log(rxy2_d) + 0.5*(M_LN2 - C_GAMMA);
    for (bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp-1] < rxy)
        break;

      double fq = C_2PI*bp;
      E += K0(fq*rxy_d)*cos(fq*z_d);
    }
    E *= 4*coulomb.prefactor*uz;
  }

  return chpref*E;
}

int mmm1d_tune(char **log)
{
  char buffer[32 + 2*ES_DOUBLE_SPACE + ES_INTEGER_SPACE];
  double int_time, min_time=1e200, min_rad = -1;
  double maxrad = box_l[2]; /* N_psi = 2, theta=2/3 maximum for rho */
  double switch_radius;

  if (mmm1d_params.far_switch_radius_2 < 0) {
    /* determine besselcutoff and optimal switching radius. Should be around 0.33 */
    for (switch_radius = 0.2*maxrad;
	 switch_radius < 0.4*maxrad;
	 switch_radius += 0.025*maxrad) {
      if (switch_radius <= bessel_radii[MAXIMAL_B_CUT - 1]) {
        // this switching radius is too small for our Bessel series
        continue;
      }

      mmm1d_params.far_switch_radius_2 = SQR(switch_radius);

      coulomb.method = COULOMB_MMM1D;
      
      /* initialize mmm1d temporary structures */
      mpi_bcast_coulomb_params();

      /* perform force calculation test */
      int_time = time_force_calc(TEST_INTEGRATIONS);

      /* exit on errors */
      if (int_time < 0)
	return ES_ERROR;

      sprintf(buffer, "r= %f t= %f ms\n", switch_radius, int_time);
      *log = strcat_alloc(*log, buffer);

      if (int_time < min_time) {
	min_time = int_time;
	min_rad = switch_radius;
      }
      /* stop if all hope is vain... */
      else if (int_time > 2*min_time)
	break;
    }
    switch_radius = min_rad;
    mmm1d_params.far_switch_radius_2 = SQR(switch_radius);
  }
  else {
    if (mmm1d_params.far_switch_radius_2 <= SQR(bessel_radii[MAXIMAL_B_CUT - 1])) {
      // this switching radius is too small for our Bessel series
      *log = strcat_alloc(*log, "could not find reasonable bessel cutoff");
      return ES_ERROR;
    }
  }

  coulomb.method = COULOMB_MMM1D;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

#endif
