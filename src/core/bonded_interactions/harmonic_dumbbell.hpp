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
#ifndef _HARMONIC_DUMBBELL_HPP
#define _HARMONIC_DUMBBELL_HPP
/** \file
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force
 *  for a particle pair.
 *  \ref forces.cpp
 */

/************************************************************/

#include "config.hpp"

#ifdef ROTATION
#include "bonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <utils/math/sqr.hpp>

/** set the parameters for the harmonic potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int harmonic_dumbbell_set_params(int bond_type, double k1, double k2, double r,
                                 double r_cut);

/** Computes the harmonic dumbbell bond length force and update torque.
 *  @param[in,out]  p1        First particle, torque gets updated.
 *  @param[in]      iaparams  Bonded parameters for the pair interaction.
 *  @param[in]      dx        %Distance between the particles.
 *  @param[out]     force     Force.
 *  @retval 0
 */
inline int
calc_harmonic_dumbbell_pair_force(Particle *p1,
                                  Bonded_ia_parameters const *iaparams,
                                  Utils::Vector3d const &dx, double *force) {
  auto const dist = dx.norm();

  if ((iaparams->p.harmonic_dumbbell.r_cut > 0.0) &&
      (dist > iaparams->p.harmonic_dumbbell.r_cut))
    return 1;

  auto const dr = dist - iaparams->p.harmonic_dumbbell.r;
  auto const normalizer = (dist > ROUND_ERROR_PREC) ? 1. / dist : 0.0;
  auto const fac = -iaparams->p.harmonic_dumbbell.k1 * dr * normalizer;

  for (int i = 0; i < 3; i++)
    force[i] = fac * dx[i];

  auto const dhat = Utils::Vector3d{dx[0], dx[1], dx[2]} * normalizer;
  auto const da = vector_product(dhat, p1->r.calc_director());

  p1->f.torque += iaparams->p.harmonic_dumbbell.k2 * da;
  return 0;
}

/** Computes the harmonic dumbbell bond length energy.
 *  @param[in]  p1        First particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int harmonic_dumbbell_pair_energy(Particle const *p1,
                                         Bonded_ia_parameters const *iaparams,
                                         Utils::Vector3d const &dx,
                                         double *_energy) {
  auto const dist = dx.norm();

  if ((iaparams->p.harmonic_dumbbell.r_cut > 0.0) &&
      (dist > iaparams->p.harmonic_dumbbell.r_cut))
    return 1;

  double dhat[3];
  dhat[0] = dx[0] / dist;
  dhat[1] = dx[1] / dist;
  dhat[2] = dx[2] / dist;

  double da[3];
  const Utils::Vector3d director1 = p1->r.calc_director();
  da[0] = dhat[1] * director1[2] - dhat[2] * director1[1];
  da[1] = dhat[2] * director1[0] - dhat[0] * director1[2];
  da[2] = dhat[0] * director1[1] - dhat[1] * director1[0];

  double torque[3];
  torque[0] = iaparams->p.harmonic_dumbbell.k2 * da[0];
  torque[1] = iaparams->p.harmonic_dumbbell.k2 * da[1];
  torque[2] = iaparams->p.harmonic_dumbbell.k2 * da[2];

  double diff[3];
  diff[0] = dhat[0] - director1[0];
  diff[1] = dhat[1] - director1[1];
  diff[2] = dhat[2] - director1[2];

  *_energy =
      0.5 * iaparams->p.harmonic_dumbbell.k1 *
          Utils::sqr(dist - iaparams->p.harmonic.r) +
      0.5 * iaparams->p.harmonic_dumbbell.k2 *
          (torque[0] * diff[0] + torque[1] * diff[1] + torque[2] * diff[2]);
  return 0;
}

#endif // ROTATION

#endif
