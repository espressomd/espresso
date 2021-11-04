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

#include "integrators/steepest_descent.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/algorithm/clamp.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

/** Currently active steepest descent instance */
static SteepestDescentParameters params{};

bool steepest_descent_step(const ParticleRange &particles) {
  // Maximal force encountered on node
  auto f_max = -std::numeric_limits<double>::max();

  // Iteration over all local particles
  for (auto &p : particles) {
    auto f = 0.0;

    // For all Cartesian coordinates
    for (int j = 0; j < 3; j++) {
      // Skip, if coordinate is fixed
      if (!(p.p.ext_flag & COORD_FIXED(j)))
        // Skip positional increments of virtual particles
        if (!p.p.is_virtual) {
          // Square of force on particle
          f += Utils::sqr(p.f.f[j]);

          // Positional increment, crop to maximum allowed by user
          auto const dp = boost::algorithm::clamp(params.gamma * p.f.f[j],
                                                  -params.max_displacement,
                                                  params.max_displacement);

          // Move particle
          p.r.p[j] += dp;
        }
    }
#ifdef ROTATION
    {
      // Rotational increment
      auto const dq = params.gamma * p.f.torque; // Vector parallel to torque
      auto const t = p.f.torque.norm2();

      // Normalize rotation axis and compute amount of rotation
      auto const l = dq.norm();
      if (l > 0.0) {
        auto const axis = dq / l;
        auto const angle = boost::algorithm::clamp(l, -params.max_displacement,
                                                   params.max_displacement);

        // Rotate the particle around axis dq by amount l
        local_rotate_particle(p, axis, angle);
      }

      f_max = std::max(f_max, t);
    }
#endif
    // Note maximum force/torque encountered
    f_max = std::max(f_max, f);
  }

  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);

  // Synchronize maximum force/torque encountered
  namespace mpi = boost::mpi;
  auto const f_max_global =
      mpi::all_reduce(comm_cart, f_max, mpi::maximum<double>());

  return sqrt(f_max_global) < params.f_max;
}

void mpi_bcast_steepest_descent_local() {
  boost::mpi::broadcast(comm_cart, params, 0);
}

REGISTER_CALLBACK(mpi_bcast_steepest_descent_local)

/** Broadcast steepest descent parameters */
void mpi_bcast_steepest_descent() {
  mpi_call_all(mpi_bcast_steepest_descent_local);
}

void steepest_descent_init(const double f_max, const double gamma,
                           const double max_displacement) {
  if (f_max < 0.0) {
    throw std::runtime_error("The maximal force must be positive.");
  }
  if (gamma < 0.0) {
    throw std::runtime_error("The dampening constant must be positive.");
  }
  if (max_displacement < 0.0) {
    throw std::runtime_error("The maximal displacement must be positive.");
  }

  params = SteepestDescentParameters{f_max, gamma, max_displacement};

  mpi_bcast_steepest_descent();
}
