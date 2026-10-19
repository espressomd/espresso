/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "integrators/steepest_descent.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

bool SteepestDescent::propagate(CellStructure &cell_structure) const {
  // Maximal force encountered on node
  auto f_max_local = -std::numeric_limits<double>::max();

  // Iteration over all local particles
  for (auto &p : cell_structure.local_particles()) {
    auto f = 0.0;

    auto const force = p.force();
    auto pos = p.pos();
    // For all Cartesian coordinates
    for (auto j = 0u; j < 3u; ++j) {
      // Skip, if coordinate is fixed
      if (!p.is_fixed_along(j)) {
        // Square of force on particle
        f += Utils::sqr(force[j]);

        // Positional increment, crop to maximum allowed by user
        auto const dp =
            std::clamp(gamma * force[j], -max_displacement, max_displacement);

        // Move particle
        pos[j] += dp;
      }
    }
#ifdef ESPRESSO_ROTATION
    {
      // Rotational increment
      auto const torque = Utils::Vector3d(p.torque());
      auto const dq = gamma * torque; // Vector parallel to torque
      auto const t = torque.norm2();

      // Normalize rotation axis and compute amount of rotation
      auto const l = dq.norm();
      if (l > 0.0) {
        auto const axis = dq / l;
        auto const angle = std::clamp(l, -max_displacement, max_displacement);

        // Rotate the particle around axis dq by amount l
        local_rotate_particle(p, axis, angle);
      }

      f_max_local = std::max(f_max_local, t);
    }
#endif
    // Note maximum force/torque encountered
    f_max_local = std::max(f_max_local, f);
  }

  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);

  // Synchronize maximum force/torque encountered
  auto const f_max_global = boost::mpi::all_reduce(
      comm_cart, f_max_local, boost::mpi::maximum<double>());

  return std::sqrt(f_max_global) < f_max;
}

SteepestDescent::SteepestDescent(double f_max, double gamma,
                                 double max_displacement)
    : f_max{f_max}, gamma{gamma}, max_displacement{max_displacement} {
  if (f_max < 0.0) {
    throw std::runtime_error("The maximal force must be positive.");
  }
  if (gamma < 0.0) {
    throw std::runtime_error("The dampening constant must be positive.");
  }
  if (max_displacement < 0.0) {
    throw std::runtime_error("The maximal displacement must be positive.");
  }
}
