// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#include <mpi.h>
#include <tcl.h>
#include "mmm1d.h"
#include "polynom.h"
#include "specfunc.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "tuning.h"
#include "utils.h"
#include "interaction_data.h"
#include "debug.h"

#ifdef ELECTROSTATICS

#define C_2PI     (2*M_PI)
#define C_GAMMA   0.57721566490153286060651209008
#define C_2LOG4PI -5.0620484939385815859557831885
#define C_2PISQR  C_2PI*C_2PI
/* precision of polygamma functions. More is unnecessary, the Bessel
   functions are not better anyways... */
#define EPS 1e-10

/* How many trial calculations */
#define TEST_INTEGRATIONS 1000

/* Largest reasonable cutoff for Bessel function */
#define MAXIMAL_B_CUT 30

/* Granularity of the radius scan in multiples of box_l[2] */
#define RAD_STEPPING 0.1

/** modified polygamma functions. See Arnold,Holm 2002 */
static Polynom *modPsi = NULL;
static int      n_modPsi = 0;
static double modPsi_curerror = 1e100;

static double L_i, L2, L2_i, prefL2_i, prefL3_i;

MMM1D_struct mmm1d_params = { 0.05, 5, 1e-5 };

MDINLINE double mod_psi_even(int n, double x)
{ return evaluateAsTaylorSeriesAt(&modPsi[2*n],x*x); }

MDINLINE double mod_psi_odd(int n, double x)
{ return x*evaluateAsTaylorSeriesAt(&modPsi[2*n+1], x*x); }

static void preparePolygammaEven(int n, double binom, Polynom *series)
{
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv = 2*n;
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    maxx = 0.25;
    alloc_doublelist(series, 1);
    series->e[0] = 2*(1 - C_GAMMA);
    for (order = 1;; order += 1) {
      x_order = 2*order;
      coeff = -2*hzeta(x_order + 1, 2);
      if (fabs(maxx*coeff)*(4.0/3.0) < EPS)
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = coeff;
      maxx *= 0.25;
    }
    series->n = order;
  }
  else {
    // even, n > 0
    maxx = 1;
    pref = 2;
    init_doublelist(series);
    for (order = 0;; order++) {
      // only even exponents of x
      x_order = 2*order;
      coeff = pref*hzeta(1 + deriv + x_order, 2);
      if ((fabs(maxx*coeff)*(4.0/3.0) < EPS) && (x_order > deriv))
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = -binom*coeff;
      maxx *= 0.25;
      pref *= (1.0 + deriv/(x_order + 1));
      pref *= (1.0 + deriv/(x_order + 2));
    }
    series->n = order;
  }
}

static void preparePolygammaOdd(int n, double binom, Polynom *series)
{
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv  = 2*n + 1;
  maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  pref = 2*deriv*(1 + deriv);
  init_doublelist(series);
  for (order = 0;; order++) {
    // only odd exponents of x
    x_order = 2*order + 1;
    coeff = pref*hzeta(1 + deriv + x_order, 2);
    if ((fabs(maxx*coeff)*(4.0/3.0) < EPS) && (x_order > deriv))
      break;
    realloc_doublelist(series, order + 1);
    series->e[order] = -binom*coeff;
    maxx *= 0.25;
    pref *= (1.0 + deriv/(x_order + 1));
    pref *= (1.0 + deriv/(x_order + 2));
  }
  series->n = order;
}

double determine_bessel_cutoff(double switch_rad, double maxPWerror, int maxP)
{
  /* this calculates an upper bound to all force components and the potential */

  // fprintf(stderr, "determ: %f %f %d\n", switch_rad, maxPWerror, maxP);

  double err;
  double Lz = box_l[2];
  double fac = 2*M_PI/Lz*switch_rad;
  double pref = 4/Lz*dmax(1, 2*M_PI/Lz);
  int P = (int)ceil(3*box_l[2]/2/M_PI/switch_rad);
  do {
    err = pref*exp(-fac*(P-1))*(1 + P*(exp(fac) - 1))/SQR(1 - exp(fac));
    // fprintf(stderr, "%d %e\n", P, err); */
    P++;
  } while (err > maxPWerror && P <= maxP);
  P--;
  return P;
}

int set_mmm1d_params(Tcl_Interp *interp, double switch_rad,
		     int bessel_cutoff, double maxPWerror)
{
  char buffer[32 + 2*TCL_DOUBLE_SPACE];
  double int_time, min_time=1e200, min_rad = -1;
  double maxrad = box_l[2];
  /* more doesn't work with the current implementation as N_psi = 2 fixed */

  if (bessel_cutoff < 0 && switch_rad < 0) {
    /* determine besselcutoff and optimal switching radius */
    for (switch_rad = RAD_STEPPING*box_l[2]; switch_rad < maxrad; switch_rad += RAD_STEPPING*box_l[2]) {
      mmm1d_params.bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
      /* no reasonable cutoff possible */
      if (mmm1d_params.bessel_cutoff == MAXIMAL_B_CUT)
	continue;
      mmm1d_params.far_switch_radius_2 = switch_rad*switch_rad;

      /* initialize mmm1d temporary structures */
      mpi_bcast_coulomb_params();
      /* perform force calculation test */
      int_time = time_force_calc(TEST_INTEGRATIONS);

      sprintf(buffer, "r= %f c= %d t= %f ms\n", switch_rad, mmm1d_params.bessel_cutoff, int_time);
      // fprintf(stderr, buffer);
      Tcl_AppendResult(interp, buffer, (char*)NULL);

      if (int_time < min_time) {
	min_time = int_time;
	min_rad = switch_rad;
      }
      /* stop if all hope is vain... */
      else if (int_time > 2*min_time)
	break;
    }
    if (min_rad < 0) {
      fprintf(stderr, "set_mmm1d_params: internal error");
      errexit();
    }
    switch_rad    = min_rad;
    bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
  }

  if (bessel_cutoff < 0) {
    /* determine besselcutoff to achieve at least the given pairwise error */
    bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
    if (bessel_cutoff == MAXIMAL_B_CUT) {
      Tcl_AppendResult(interp, "could not find reasonable bessel cutoff", (char *)NULL);
      return TCL_ERROR;
    }
  }

  if (switch_rad <= 0 || switch_rad > box_l[2]) {
    Tcl_AppendResult(interp, "switching radius is not between 0 and box_l[2]", (char *)NULL);
    return TCL_ERROR;
  }
  if (bessel_cutoff <=0) {
    Tcl_AppendResult(interp, "bessel cutoff too small", (char *)NULL);
    return TCL_ERROR;
  }

  mmm1d_params.far_switch_radius_2 = switch_rad*switch_rad;
  mmm1d_params.bessel_cutoff = bessel_cutoff;

  mmm1d_params.maxPWerror = maxPWerror;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

void MMM1D_recalcTables()
{
  /* polygamma, determine order */
  int n;
  double binom, err;
  double rho2m2max;

  /* our tables are even better */
  if (modPsi_curerror < mmm1d_params.maxPWerror)
    return;

  if (modPsi != NULL) {
    for (n = 0; n < 2*n_modPsi; n++)
      realloc_doublelist(&modPsi[n], 0);
    modPsi = realloc(modPsi, n_modPsi = 0);
  }

  n = 0;
  binom = 1.0;
  rho2m2max = 1.0;
  do {
    n_modPsi++;
    modPsi = realloc(modPsi, 2*n_modPsi*sizeof(Polynom));
    
    preparePolygammaEven(n, binom, &modPsi[2*n_modPsi - 2]);
    preparePolygammaOdd(n, binom, &modPsi[2*n_modPsi - 1]);
    
    err = fabs(2*mod_psi_even(n,rho2m2max));
    rho2m2max *= 0.5;
    binom *= (-0.5 - n)/(double)(n+1);
    n++;
  }
  while (err > 0.1*mmm1d_params.maxPWerror);

  modPsi_curerror = mmm1d_params.maxPWerror;
}

void MMM1D_init()
{
  if (PERIODIC(0) != 0 || PERIODIC(1) != 0 || PERIODIC(2) != 1) {
    fprintf(stderr, "ERROR: MMM1D requires periodicity 0 0 1\n");
    errexit();
  }

  if (cell_structure.type != CELL_STRUCTURE_NSQUARE) {
    fprintf(stderr, "ERROR: MMM1D requires n-square cellsystem\n");
    errexit();
  }

  /* precalculate some constants */
  L_i  = 1/box_l[2];
  L2   = box_l[2]*box_l[2];
  L2_i = L_i*L_i;
  prefL2_i = coulomb.prefactor*L2_i;
  prefL3_i = prefL2_i*L_i;

  MMM1D_recalcTables();
}

void add_mmm1d_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double r2, double r)
{
  int dim;
  double F[3];
  double chpref = p1->p.q*p2->p.q;
  double rxy2, rxy2_d, z_d;
  double pref;
  double Fx, Fy, Fz;
  
  if (chpref == 0)
    return;

  rxy2   = d[0]*d[0] + d[1]*d[1];
  rxy2_d = rxy2*L2_i;
  z_d    = d[2]*L_i;

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

    Fx = prefL3_i*sr*d[0];
    Fy = prefL3_i*sr*d[1];
    Fz = prefL2_i*sz;

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
    double rxy_d = rxy*L_i;
    double sr = 0, sz = 0;
    int bp;

    for (bp = 1; bp < mmm1d_params.bessel_cutoff; bp++) {
      double fq = C_2PI*bp;    
      sr += bp*K1(fq*rxy_d)*cos(fq*z_d);
      sz += bp*K0(fq*rxy_d)*sin(fq*z_d);
    }
    sr *= L2_i*4*C_2PI;
    sz *= L2_i*4*C_2PI;
    
    pref = coulomb.prefactor*(sr/rxy + 2*L_i/rxy2);

    F[0] = pref*d[0];
    F[1] = pref*d[1];
    F[2] = coulomb.prefactor*sz;
  }

  for (dim = 0; dim < 3; dim++)
    F[dim] *= chpref;
  for (dim = 0; dim < 3; dim++) {
    p1->f.f[dim] += F[dim];
    p2->f.f[dim] -= F[dim];
  }
}

double mmm1d_coulomb_pair_energy(Particle *p1, Particle *p2, double d[3], double r2, double r)
{
  double chpref = p1->p.q*p2->p.q;
  double rxy2, rxy2_d, z_d;
  double E;

  if (chpref == 0)
    return 0;

  rxy2   = d[0]*d[0] + d[1]*d[1];
  rxy2_d = rxy2*L2_i;
  z_d    = d[2]*L_i;

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* near range formula */
    double r2n, rt, shift_z;
    int n;

    E = -C_GAMMA*L_i;

    /* polygamma summation */
    r2n = 1.0;
    for (n = 0; n < n_modPsi; n++) {
      double add = mod_psi_even(n, z_d)*r2n;
      E -= add;
      
      if (fabs(add) < mmm1d_params.maxPWerror)
	break;
 
      r2n *= rxy2_d;
    }

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
    double rxy_d = rxy*L_i;
    E = 0;
    int bp;

    for (bp = 1; bp < mmm1d_params.bessel_cutoff; bp++) {
      double fq = C_2PI*bp;    
      E += K0(fq*rxy_d)*cos(fq*z_d);
    }
    E *= coulomb.prefactor*L_i*4;

    E += -L_i*log(0.25*rxy2_d) - C_GAMMA*L_i;
  }

  return chpref*E;
}

#endif
