/*
  Copyright (C) 2010-2018 The ESPResSo project
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
#ifndef ANGLE_COSSQUARE_H
#define ANGLE_COSSQUARE_H
/** \file
 *  Routines to calculate the angle energy or/and and force
 *  for a particle triple.
 *  \ref forces.cpp
 */

#include "bonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

#ifdef BOND_ANGLE
#include "angle_common.hpp"
#include "grid.hpp"

/** set parameters for the angle potential.

    \todo The type of the angle potential
    is chosen via config.hpp and cannot be changed at runtime.
*/
int angle_cossquare_set_params(int bond_type, double bend, double phi0);

/************************************************************/

/** Computes the three body angle interaction force and adds this
    force to the particle forces.
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction.
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
*/
inline int calc_angle_cossquare_force(Particle const *p_mid,
                                      Particle const *p_left,
                                      Particle const *p_right,
                                      Bonded_ia_parameters const *iaparams,
                                      double force1[3], double force2[3]) {
  double cosine, vec1[3], vec2[3], d1i, d2i, fac;

  /* vector from p_left to p_mid */
  calc_angle_vector(p_mid->r.p, p_left->r.p, vec1, d1i);
  /* vector from p_mid to p_right */
  calc_angle_vector(p_right->r.p, p_mid->r.p, vec2, d2i);
  /* scalar product of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  fac = iaparams->p.angle_cossquare.bend;

  fac *= iaparams->p.angle_cossquare.cos_phi0 + cosine;

  calc_angle_force(force1, force2, vec1, vec2, d1i, d2i, cosine, fac);

  return 0;
}

/* The force on each particle due to a three-body bonded potential
   is computed. */
inline void calc_angle_cossquare_3body_forces(Particle const *p_mid,
                                              Particle const *p_left,
                                              Particle const *p_right,
                                              Bonded_ia_parameters const *iaparams,
                                              double force1[3],
                                              double force2[3],
                                              double force3[3]) {
  double pot_dep;
  double cos_phi;
  double sin_phi;
  double vec21[3];
  double vec31[3];
  double vec21_sqr;
  double vec31_sqr;
  double vec21_magn;
  double vec31_magn;
  double fac;

  calc_angle_3body_vector(p_mid->r.p, p_left->r.p, p_right->r.p,
                          cos_phi, sin_phi, vec21, vec31,
                          vec21_sqr, vec31_sqr, vec21_magn, vec31_magn);

  /* uncomment this block if interested in the angle
  if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
  if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE;
  phi = acos(cos_phi);
  */
  {
    double K, cos_phi0;
    K = iaparams->p.angle_cossquare.bend;
    cos_phi0 = iaparams->p.angle_cossquare.cos_phi0;

    // potential dependent term [dU/dphi = K * (sin_phi * cos_phi0 - cos_phi *
    // sin_phi)]
    pot_dep = K * (sin_phi * cos_phi0 - cos_phi * sin_phi);
  }

  fac = pot_dep / sin_phi;

  calc_angle_3body_force(cos_phi, fac, vec21, vec31, vec21_sqr, vec31_sqr,
                         vec21_magn, vec31_magn, force1, force2, force3);
}

/** Computes the three body angle interaction energy.
    @param p_mid        Pointer to first particle.
    @param p_left        Pointer to second/middle particle.
    @param p_right        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction.
    @param _energy   return energy pointer.
    @return 0.
*/
inline int angle_cossquare_energy(Particle const *p_mid,
                                  Particle const *p_left,
                                  Particle const *p_right,
                                  Bonded_ia_parameters const *iaparams,
                                  double *_energy) {
  double cosine, vec1[3], vec2[3], d1i, d2i;

  /* vector from p_left to p_mid */
  calc_angle_vector(p_mid->r.p, p_left->r.p, vec1, d1i);
  /* vector from p_mid to p_right */
  calc_angle_vector(p_right->r.p, p_mid->r.p, vec2, d2i);
  /* scalar product of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;
  /* bond angle energy */
  *_energy = 0.5 * iaparams->p.angle_cossquare.bend *
             Utils::sqr(cosine + iaparams->p.angle_cossquare.cos_phi0);
  return 0;
}

#endif /* BOND_ANGLE */
#endif /* ANGLE_COSSQUARE_H */
