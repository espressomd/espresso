/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MEMBRANE_COLLISION_H
#define MEMBRANE_COLLISION_H

/** \file
 *  Routines to calculate the membrane collision force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef MEMBRANE_COLLISION

#include "Particle.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/Vector.hpp>

int membrane_collision_set_params(int part_type_a, int part_type_b, double a,
                                  double n, double cut, double offset);

/** Resultant force due to a sigmoid potential between two
 *  particles at interatomic separation r
 */
inline double sigmoid_force_r(double a, double n, double r) {
  return (a / (1 + exp(n * r)));
}

/** Calculate membrane-collision force between particle p1 and p2 */
inline Utils::Vector3d
membrane_collision_pair_force(Particle const &p1, Particle const &p2,
                              IA_parameters const &ia_params,
                              Utils::Vector3d const &d, double dist) {
  /************************
   *
   * Description of implementation:
   * We have two particles, each belongs to the membrane of a different immersed
   *object For both particles we have the position of the particle - p, and in
   *part.p.out_direction are the coordinates of the outward normal vector (with
   *respect to the immersed object).
   *
   * Algorithm:
   * 1) compute angle between outward normal vectors out1 and out2 to check if
   *they are almost collinear (in this case nothing happens, this is at the edge
   *of contact zone and forces should not be applied because they would quite
   *likely have wrong direction) 2) in other cases, repulsive forces are applied
   *in the direction out1-out2
   *
   *********************/
  Utils::Vector3d force{};
  if (dist < (ia_params.membrane.cut + ia_params.membrane.offset)) {

    auto const r_off = dist - ia_params.membrane.offset;
    // offset needs to be checked for the unphysical case when r_off should be
    // negative

    if (r_off > 0.0) {

      auto const out1 = p1.p.out_direction;
      auto const out2 = p2.p.out_direction;
      // check whether out_direction was set
      if (fabs(out1[0]) + fabs(out1[1]) + fabs(out1[2]) + fabs(out2[0]) +
              fabs(out2[1]) + fabs(out2[2]) <
          SMALL_OIF_MEMBRANE_CUTOFF) {
        throw std::runtime_error("out direction not set");
      }

      // this is the direction in which the repulsive forces will be applied and
      // its norm
      auto const dir = out1 - out2;
      auto const ndir = dir.norm();

      // for very small angles the force should not be applied - these happen at
      // the crossing of the boundary and would result in oscillation
      auto const product = out1 * out2;
      auto const angle = acos(product);

      if (fabs(angle) > SMALL_OIF_MEMBRANE_CUTOFF) {
        auto const fac =
            sigmoid_force_r(ia_params.membrane.a, ia_params.membrane.n, r_off) /
            dist;
        force = (-fac / ndir) * dir;
      }
    }
  }
  return force;
}

#endif /* ifdef MEMBRANE_COLLISION */
#endif
