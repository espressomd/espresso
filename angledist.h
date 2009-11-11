// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef ANGLEDIST_H
#define ANGLEDIST_H
/** \file angledist.h
 *  Routines to calculate the angle and distance dependent (from a constraint) energy or/and and force
 *  for a particle triple.
 *  \ref forces.c
*/

#ifdef BOND_ANGLEDIST

#include "utils.h"

/************************************************************/

/** set parameters for the angledist potential. The type of the angledist potential
    is chosen via config.h and cannot be changed at runtime.
**/
MDINLINE int angledist_set_params(int bond_type, double bend, double phimin, double distmin, double phimax, double distmax)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angledist.bend = bend;
  bonded_ia_params[bond_type].p.angledist.phimin = phimin;
  bonded_ia_params[bond_type].p.angledist.distmin = distmin;
  bonded_ia_params[bond_type].p.angledist.phimax = phimax;
  bonded_ia_params[bond_type].p.angledist.distmax = distmax;
#ifdef BOND_ANGLEDIST_COSINE
#error angledist not implemented for BOND_ANGLEDIST_COSINE
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
#error angledist not implemented for BOND_ANGLEDIST_COSSQUARE
#endif
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLEDIST;
  bonded_ia_params[bond_type].num = 2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the angle potential
MDINLINE int inter_parse_angledist(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phimin, distmin, phimax, distmax;

  if (argc != 6) {
    Tcl_AppendResult(interp, "angledist needs 5 parameters: "
		     "<bend> <phimin> <distmin> <phimax> <distmax>", (char *) NULL);
    printf ("argc=%d\n",argc);
    return (TCL_ERROR);
  }

  if (! ARG_IS_D(1, bend)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<bend> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(2, phimin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(3, distmin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(4, phimax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimax> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(5, distmax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmax> ", (char *) NULL);
    return TCL_ERROR;
  }


  CHECK_VALUE(angledist_set_params(bond_type, bend, phimin, distmin, phimax, distmax), "bond type must be nonnegative");
}

// Function to calculate wall distance and phi0(dist)
// Called by calc_angledist_force
MDINLINE double calc_angledist_param(Particle *p_mid, Particle *p_left, Particle *p_right, 
                       Bonded_ia_parameters *iaparams)
{
  double cosine=0.0, vec1[3], vec2[3], d1i=0.0, d2i=0.0, dist1=0.0, dist2=0.0, fac=0.0, phi0=0.0;
  double pwdist=0.0, pwdist0=0.0, pwdist1=0.0, normal, force[3], folded_pos[3], phimn=0.0, distmn=0.0, phimx=0.0, distmx=0.0, drange=0.0;
  Constraint_wall wall;
  int j, k;
  int img[3];

  cosine=0.0;
  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* vectors are normalised so cosine is just cos(angle_between_vec1_and_vec2) */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  fac    = iaparams->p.angledist.bend; /* spring constant from .tcl file */
  phimn  = iaparams->p.angledist.phimin;
  distmn = iaparams->p.angledist.distmin;
  phimx  = iaparams->p.angledist.phimax;
  distmx = iaparams->p.angledist.distmax;

  /* folds coordinates of p_mid into original box */
  memcpy(folded_pos, p_mid->r.p, 3*sizeof(double));
  memcpy(img, p_mid->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);

  /* Calculates distance between p_mid and constraint */
  for(k=0;k<n_constraints;k++) {
    for (j=0; j<3; j++) {
      force[j] = 0;
    }
    switch(constraints[k].type) {
    case CONSTRAINT_WAL: 

      /* dist is distance of wall from origin */
      wall=constraints[k].c.wal;

      /* check that constraint vector is normalised */
      normal=0.0;
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }

      /* pwdist is distance of wall from p_mid */
      pwdist0=-1.0 * constraints[0].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist0 += folded_pos[j] * constraints[0].c.wal.n[j];
      }
      pwdist1=-1.0 * constraints[1].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist1 += folded_pos[j] * constraints[1].c.wal.n[j];
      }
      if (pwdist0 <= pwdist1) {
        pwdist = pwdist0;
      }
      else {
        pwdist = pwdist1;
      }

      /*get phi0(z)*/
      if (pwdist <= distmn) {
        phi0 = phimn;
      }
      else if (pwdist >= distmx && pwdist <= box_l[2]-wall.d-distmx) {
          phi0 = phimx;
      }
      else {
        drange = (pwdist-distmn)*PI/(distmx-distmn);
        phi0 = ((cos(drange-PI)+1.0)*(phimx-phimn))*0.5+phimn;
      }
      break;
    }
  }
  return phi0;
}


MDINLINE int calc_angledist_force(Particle *p_mid, Particle *p_left, Particle *p_right,
                                  Bonded_ia_parameters *iaparams, double force1[3], double force2[3])
{
  double cosine=0.0, vec1[3], vec2[3], d1i=0.0, d2i=0.0, dist2, fac=0.0, f1=0.0, f2=0.0, phi0=0.0;
  int j;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  /* NOTE The angledist is ONLY implemented for the HARMONIC case */
  phi0=calc_angledist_param(p_mid, p_left, p_right, iaparams);

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi=0.0,sinphi=0.0;
    if ( cosine >  TINY_COS_VALUE) cosine =  TINY_COS_VALUE;
    if ( cosine < -TINY_COS_VALUE) cosine = -TINY_COS_VALUE;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi - phi0)/sinphi;
  }
#endif

#ifdef BOND_ANGLEDIST_COSINE
  #error angledist ONLY implemented for harmonic case!
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
  #error angledist ONLY implemented for harmonic case!
#endif

  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force1[j] = (f1-f2);
    force2[j] = -f1;
  }
  return 0;
}

/** Computes the three body angle interaction energy (see \ref #inter, \ref #analyze). 
    @param p_mid        Pointer to second/middle particle.
    @param p_left       Pointer to first particle.
    @param p_right      Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param _energy   return energy pointer.
    @return 0.
*/



MDINLINE int angledist_energy(Particle *p_mid, Particle *p_left, Particle *p_right, 
       			      Bonded_ia_parameters *iaparams, double *_energy)
{
  int j;
  double cosine=0.0, d1i, d2i, dist1, dist2;
  double phi0;
  double vec1[3], vec2[3];

  phi0=calc_angledist_param(p_mid, p_left, p_right, iaparams);
  /*  fprintf(stdout,"\nIn angledist_energy:  phi0=%f\n",phi0*180.0/PI);*/

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  phi0=calc_angledist_param(p_mid, p_left, p_right, iaparams);

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi;
    phi =  acos(-cosine);
    *_energy = 0.5*iaparams->p.angledist.bend*SQR(phi - phi0);
  }
#endif
#ifdef BOND_ANGLEDIST_COSINE
  #error angledist ONLY implemented for harmonic case!
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
  #error angledist ONLY implemented for harmonic case!
#endif
  return 0;
}

#endif /* BOND_ANGLEDIST */
#endif /* ANGLEDIST_H */
