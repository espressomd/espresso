/*
  Copyright (C) 2012-2018 The ESPResSo project

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
#ifndef _OBJECT_IN_FLUID_OUT_DIRECTION_H
#define _OBJECT_IN_FLUID_OUT_DIRECTION_H
/** \file
 * Routines to calculate the outward direction of the
 * membrane using a particle quadruple (one particle and its 3 strategically
 * placed neighbors)
 */

#include "config.hpp"

#ifdef MEMBRANE_COLLISION

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/math/triangle_functions.hpp>
using Utils::get_n_triangle;

// set out_direction parameters
int oif_out_direction_set_params(int bond_type);

/** Computes the outward direction of the membrane from one particle and its
 *  three neighbors
 *
 *  @param[out] p1            The central particle.
 *  @param[in]  p2 , p3 , p4  The neighboring particles.
 *
 *  Computes the normal of triangle p2p3p4. This triangle was initially
 *  oriented in such a way that its normal already points out of the object.
 *  Normalizes and stores the result as @ref ParticleProperties::out_direction
 *  "out_direction" in @p p1.
 *  @return 0
 */
inline int calc_out_direction(
    Particle *const p1, Particle const *const p2, Particle const *const p3,
    Particle const *const p4,
    Bonded_ia_parameters * /* iaparams */) // first-fold-then-the-same approach
{
  auto const fp2 = unfolded_position(*p2);
  auto const fp3 = fp2 + get_mi_vector(p3->r.p, fp2);
  auto const fp4 = fp2 + get_mi_vector(p4->r.p, fp2);

  auto const n = get_n_triangle(fp2, fp3, fp4);
  auto const dn = n.norm();

  if (std::abs(dn) < 0.001) {
    printf("out_direction.hpp, calc_out_direction: Length of outward vector is "
           "close to zero!\n");
  }

  p1->p.out_direction = n / dn;

  return 0;
}

#endif /* ifdef MEMBRANE_COLLISION */
#endif
