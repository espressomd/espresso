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

#ifdef BOND_ANGLE
#include "angle_common.hpp"
#include "grid.hpp"
#include <tuple>

/** set parameters for the angle potential.

    \todo The type of the angle potential
    is chosen via config.hpp and cannot be changed at runtime.
*/
int angle_cossquare_set_params(int bond_type, double bend, double phi0);

/************************************************************/

/** Compute the three-body angle interaction force.
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] force1    Force on particle 1.
 *  @param[out] force2    Force on particle 2.
 *  @retval 0
 */
inline int calc_angle_cossquare_force(Particle const *p_mid,
                                      Particle const *p_left,
                                      Particle const *p_right,
                                      Bonded_ia_parameters const *iaparams,
                                      double force1[3], double force2[3]) {

  auto forceFactor = [&iaparams](double const cos_phi) {
    auto const K = iaparams->p.angle_cossquare.bend;
    auto const cos_phi0 = iaparams->p.angle_cossquare.cos_phi0;
    return K * (cos_phi0 + cos_phi);
  };

  calc_angle_generic_force(p_mid->r.p, p_left->r.p, p_right->r.p, forceFactor,
                           force1, force2, false);

  return 0;
}

/* The force on each particle due to a three-body bonded potential
   is computed. */
inline void calc_angle_cossquare_3body_forces(
    Particle const *p_mid, Particle const *p_left, Particle const *p_right,
    Bonded_ia_parameters const *iaparams, Vector3d &force1, Vector3d &force2,
    Vector3d &force3) {

  auto forceFactor = [&iaparams](double const cos_phi, double const sin_phi) {
    auto const K = iaparams->p.angle_cossquare.bend;
    auto const cos_phi0 = iaparams->p.angle_cossquare.cos_phi0;
    // potential dependent term [dU/dphi = K * (sin_phi * cos_phi0 - cos_phi *
    // sin_phi)]
    return K * (sin_phi * cos_phi0 - cos_phi * sin_phi) / sin_phi;
  };

  std::tie(force1, force2, force3) = calc_angle_generic_3body_forces(
      p_mid->r.p, p_left->r.p, p_right->r.p, forceFactor, false);
}

/** Computes the three-body angle interaction energy.
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int angle_cossquare_energy(Particle const *p_mid, Particle const *p_left,
                                  Particle const *p_right,
                                  Bonded_ia_parameters const *iaparams,
                                  double *_energy) {
  auto const vectors =
      calc_vectors_and_cosine(p_mid->r.p, p_left->r.p, p_right->r.p, true);
  auto const cos_phi = std::get<4>(vectors);
  auto const cos_phi0 = iaparams->p.angle_cossquare.cos_phi0;
  auto const K = iaparams->p.angle_cossquare.bend;
  *_energy = 0.5 * K * Utils::sqr(cos_phi + cos_phi0);
  return 0;
}

#endif /* BOND_ANGLE */
#endif /* ANGLE_COSSQUARE_H */
