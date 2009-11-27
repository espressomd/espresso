// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef ANGLE_H
#define ANGLE_H
/** \file angle.h
 *  Routines to calculate the angle energy or/and and force 
 *  for a particle triple.
 *  \ref forces.c
*/

#ifdef BOND_ANGLE

#include "utils.h"

/************************************************************/

/** set parameters for the angle potential. The type of the angle potential
    is chosen via config.h and cannot be changed at runtime.
*/
MDINLINE int angle_set_params(int bond_type, double bend, double phi0)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angle.bend = bend;
  bonded_ia_params[bond_type].p.angle.phi0 = phi0;
#ifdef BOND_ANGLE_COSINE
  bonded_ia_params[bond_type].p.angle.cos_phi0 = cos(phi0);
  bonded_ia_params[bond_type].p.angle.sin_phi0 = sin(phi0);
#endif
#ifdef BOND_ANGLE_COSSQUARE
  bonded_ia_params[bond_type].p.angle.cos_phi0 = cos(phi0);
#endif
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLE;
  bonded_ia_params[bond_type].num = 2;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the angle potential
MDINLINE int inter_parse_angle(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phi0;

  /* the optional parameter phi0 is due to backwards compatibility and is set to PI if not given */
  if (argc != 2 && argc != 3) {
    Tcl_AppendResult(interp, "angle needs 1 or 2 parameters: "
		     "<bend> [<phi0>]", (char *) NULL);
    return (TCL_ERROR);
  }

  if (! ARG_IS_D(1, bend)) {
    Tcl_AppendResult(interp, "angle needs a DOUBLE parameter: "
		     "<bend> ", (char *) NULL);
    return TCL_ERROR;
  }

  /* special treatment of the optional parameter phi0 */
  if (argc == 3) {
    if (! ARG_IS_D(2, phi0)) {
      Tcl_AppendResult(interp, "angle needs a DOUBLE parameter: "
		       "<phi0> ", (char *) NULL);
      return TCL_ERROR;
    }
  } else {
    phi0 = PI;
  }
  CHECK_VALUE(angle_set_params(bond_type, bend, phi0), "bond type must be nonnegative");
}

/** Computes the three body angle interaction force and adds this
    force to the particle forces (see \ref #inter). 
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
*/
MDINLINE int calc_angle_force(Particle *p_mid, Particle *p_left, Particle *p_right,
			      Bonded_ia_parameters *iaparams, double force1[3], double force2[3])
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j;

  cosine=0.0;
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
  fac    = iaparams->p.angle.bend;

#ifdef BOND_ANGLE_HARMONIC
  {
    double phi,sinphi;
    if ( cosine >  TINY_COS_VALUE) cosine = TINY_COS_VALUE;
    if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi - iaparams->p.angle.phi0)/sinphi;
  }
#endif
#ifdef BOND_ANGLE_COSINE
  if ( cosine >  TINY_COS_VALUE ) cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  fac *= iaparams->p.angle.sin_phi0 * (cosine/sqrt(1-SQR(cosine))) + iaparams->p.angle.cos_phi0;
#endif
#ifdef BOND_ANGLE_COSSQUARE
  fac *= iaparams->p.angle.cos_phi0 + cosine;
#endif
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force1[j] = (f1-f2);
    force2[j] = -f1;
  }
  return 0;
}

/* The force on each particle due to a three-body bonded potential
   is computed. */
MDINLINE void calc_angle_3body_forces(Particle *p_mid, Particle *p_left,
              Particle *p_right, Bonded_ia_parameters *iaparams,
              double force1[3], double force2[3], double force3[3]) {

  int j;
  double pot_dep;
  double cos_phi;
  double sin_phi;
  double vec31[3];
  double vec21[3];
  double vec12[3]; // espresso convention
  double vec21_sqr;
  double vec31_sqr;
  double vec21_magn;
  double vec31_magn;
  double fj[3];
  double fk[3];
  double K, sin_phi0, cos_phi0;

  get_mi_vector(vec12, p_mid->r.p, p_left->r.p);
  for(j = 0; j < 3; j++)
    vec21[j] = -vec12[j];

  get_mi_vector(vec31, p_right->r.p, p_mid->r.p);
  vec21_sqr = sqrlen(vec21);
  vec21_magn = sqrt(vec21_sqr);
  vec31_sqr = sqrlen(vec31);
  vec31_magn = sqrt(vec31_sqr);
  cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  sin_phi = sqrt(1.0 - SQR(cos_phi));

  /* uncomment this block if interested in the angle 
  if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
  if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE; 
  phi = acos(cos_phi);
  */

  K = iaparams->p.angle.bend;
  sin_phi0 = iaparams->p.angle.sin_phi0;
  cos_phi0 = iaparams->p.angle.cos_phi0;

  // potential dependent term (dU/dphi * 1 / sin(phi))
  pot_dep = K * (sin_phi * cos_phi0 - cos_phi * sin_phi0) / sin_phi;

  for(j = 0; j < 3; j++) {
    fj[j] = vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] = vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3)
  for(j = 0; j < 3; j++) {
    force1[j] = force1[j] - pot_dep * (fj[j] + fk[j]);
    force2[j] = force2[j] + pot_dep * fj[j];
    force3[j] = force3[j] + pot_dep * fk[j];
  }
}


/** Computes the three body angle interaction energy (see \ref #inter, \ref #analyze). 
    @param p_mid        Pointer to first particle.
    @param p_left        Pointer to second/middle particle.
    @param p_right        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param _energy   return energy pointer.
    @return 0.
*/
MDINLINE int angle_energy(Particle *p_mid, Particle *p_left, Particle *p_right,
			     Bonded_ia_parameters *iaparams, double *_energy)
{
  double cosine, vec1[3], vec2[3],  d1i, d2i, dist2;
  int j;

  cosine=0.0;
  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
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
  /* bond angle energy */
#ifdef BOND_ANGLE_HARMONIC
  {
    double phi;
    phi =  acos(-cosine);
    *_energy = 0.5*iaparams->p.angle.bend*SQR(phi - iaparams->p.angle.phi0);
  }
#endif
#ifdef BOND_ANGLE_COSINE
  *_energy = iaparams->p.angle.bend*(cosine*iaparams->p.angle.cos_phi0 - sqrt(1-SQR(cosine))*iaparams->p.angle.sin_phi0+1);
#endif
#ifdef BOND_ANGLE_COSSQUARE
  *_energy = 0.5*iaparams->p.angle.bend*SQR(cosine + iaparams->p.angle.cos_phi0);
#endif
  return 0;
}

#endif /* BOND_ANGLE */
#endif /* ANGLE_H */
