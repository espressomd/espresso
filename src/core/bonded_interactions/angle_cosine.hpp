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
#ifndef ANGLE_COSINE_H
#define ANGLE_COSINE_H
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
 *
 *  \todo The type of the angle potential
 *  is chosen via config.hpp and cannot be changed at runtime.
 */
int angle_cosine_set_params(int bond_type, double bend, double phi0);

/************************************************************/

/** Computes the three-body angle interaction force.
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] force1    Force on particle 1.
 *  @param[out] force2    Force on particle 2.
 *  @retval 0
 */
inline int calc_angle_cosine_force(Particle const *p_mid,
                                   Particle const *p_left,
                                   Particle const *p_right,
                                   Bonded_ia_parameters const *iaparams,
                                   double force1[3], double force2[3]) {

  auto forceFactor = [&iaparams](double &cosine) {
    auto fac = iaparams->p.angle_cosine.bend;
    if (cosine > TINY_COS_VALUE)
      cosine = TINY_COS_VALUE;
    if (cosine < -TINY_COS_VALUE)
      cosine = -TINY_COS_VALUE;
    fac *= iaparams->p.angle_cosine.sin_phi0 *
               (cosine / sqrt(1 - Utils::sqr(cosine))) +
           iaparams->p.angle_cosine.cos_phi0;
    return fac;
  };

  calc_angle_generic_force(p_mid, p_left, p_right, forceFactor, force1, force2);

  return 0;
}

/* The force on each particle due to a three-body bonded potential
   is computed. */
inline void calc_angle_cosine_3body_forces(Particle const *p_mid,
                                           Particle const *p_left,
                                           Particle const *p_right,
                                           Bonded_ia_parameters const *iaparams,
                                           Vector3d &force1, Vector3d &force2,
                                           Vector3d &force3) {

  auto forceFactor = [&iaparams](double &cos_phi, double &sin_phi) {
    /* uncomment this block if interested in the angle
    if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
    if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE;
    phi = acos(cos_phi);
    */
    auto K = iaparams->p.angle_cosine.bend;
    auto sin_phi0 = iaparams->p.angle_cosine.sin_phi0;
    auto cos_phi0 = iaparams->p.angle_cosine.cos_phi0;
    // potential dependent term [dU/dphi = K * sin(phi - phi0)]
    // trig identity: sin(a - b) = sin(a)cos(b) - cos(a)sin(b)
    auto pot_dep = K * (sin_phi * cos_phi0 - cos_phi * sin_phi0);
    auto fac = pot_dep / sin_phi;
    return fac;
  };

  calc_angle_generic_3body_forces(p_mid, p_left, p_right, forceFactor, force1,
                                  force2, force3);
}

/** Computes the three-body angle interaction energy.
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int angle_cosine_energy(Particle const *p_mid, Particle const *p_left,
                               Particle const *p_right,
                               Bonded_ia_parameters const *iaparams,
                               double *_energy) {
  /* vector from p_left to p_mid */
  auto vec1 = get_mi_vector(p_mid->r.p, p_left->r.p);
  double d1i = 1.0 / vec1.norm();
  vec1 *= d1i;
  /* vector from p_mid to p_right */
  auto vec2 = get_mi_vector(p_right->r.p, p_mid->r.p);
  double d2i = 1.0 / vec2.norm();
  vec2 *= d2i;
  /* scalar product of vec1 and vec2 */
  double cosine = scalar(vec1, vec2);
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;
  /* bond angle energy */

  *_energy =
      iaparams->p.angle_cosine.bend *
      (cosine * iaparams->p.angle_cosine.cos_phi0 -
       sqrt(1 - Utils::sqr(cosine)) * iaparams->p.angle_cosine.sin_phi0 + 1);

  return 0;
}

#endif /* BOND_ANGLE */
#endif /* ANGLE_COSINE_H */
