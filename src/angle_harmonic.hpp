/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef ANGLE_HARMONIC_H
#define ANGLE_HARMONIC_H
/** \file angle_harmonic.hpp
 *  Routines to calculate the angle energy or/and and force 
 *  for a particle triple.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"

#ifdef BOND_ANGLE
#include "grid.hpp"

/** set parameters for the angle potential.

    \todo The type of the angle potential
    is chosen via config.hpp and cannot be changed at runtime.
*/
int angle_harmonic_set_params(int bond_type, double bend, double phi0);

/************************************************************/

/** Computes the three body angle interaction force and adds this
    force to the particle forces (see \ref tclcommand_inter). 
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
*/
inline int calc_angle_harmonic_force(Particle *p_mid, Particle *p_left, Particle *p_right,
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
  fac    = iaparams->p.angle_harmonic.bend;

  {
    double phi,sinphi;
    if ( cosine >  TINY_COS_VALUE) cosine = TINY_COS_VALUE;
    if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi - iaparams->p.angle_harmonic.phi0)/sinphi;
  }
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
inline void calc_angle_harmonic_3body_forces(Particle *p_mid, Particle *p_left,
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
  double fac;

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
  {
    double K, phi, phi0;
    if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
    if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE;
    phi = acos(cos_phi);

    K = iaparams->p.angle_harmonic.bend;
    phi0 = iaparams->p.angle_harmonic.phi0;

    // potential dependent term [dU/dphi = K * (phi - phi0)]
    pot_dep = K * (phi - phi0);
  }

  fac = pot_dep / sin_phi;

  for(j = 0; j < 3; j++) {
    fj[j] = vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] = vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3)
  for(j = 0; j < 3; j++) {
    force1[j] = force1[j] - fac * (fj[j] + fk[j]);
    force2[j] = force2[j] + fac * fj[j];
    force3[j] = force3[j] + fac * fk[j];
  }
}


/** Computes the three body angle interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p_mid        Pointer to first particle.
    @param p_left        Pointer to second/middle particle.
    @param p_right        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
*/
inline int angle_harmonic_energy(Particle *p_mid, Particle *p_left, Particle *p_right,
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
  {
    double phi;
    phi =  acos(-cosine);
    *_energy = 0.5*iaparams->p.angle_harmonic.bend*SQR(phi - iaparams->p.angle_harmonic.phi0);
  }
  return 0;
}

#endif /* BOND_ANGLE */
#endif /* ANGLE_HARMONIC_H */
