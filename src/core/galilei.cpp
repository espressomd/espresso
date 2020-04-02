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
/** \file
 *
 */

#include "galilei.hpp"
#include "cells.hpp"
#include "event.hpp"
#include "grid.hpp"

#include <boost/range/numeric.hpp>

/* Stop the particle motion by setting the
   velocity of each particle to zero */
void local_kill_particle_motion(int omega, const ParticleRange &particles) {
  for (auto &p : particles) {
    if (omega) {
      p.m = {};
    } else {
      p.m.v = {};
    }
  }
}

/* Set all the forces acting on the particles
   to zero */
void local_kill_particle_forces(int torque, const ParticleRange &particles) {
  for (auto &p : particles) {
    if (torque) {
      p.f = {};
    } else {
      p.f.f = {};
    }
  }
}

/* Calculate the CMS of the system */
std::pair<Utils::Vector3d, double> local_system_CMS() {
  return boost::accumulate(
      cell_structure.local_particles(), std::pair<Utils::Vector3d, double>{},
      [](auto sum, const Particle &p) {
        if (not p.p.is_virtual) {
          return std::pair<Utils::Vector3d, double>{
              sum.first +
                  p.p.mass * unfolded_position(p.r.p, p.l.i, box_geo.length()),
              sum.second + p.p.mass};
        }
        return std::pair<Utils::Vector3d, double>{sum.first, sum.second};
      });
}

std::pair<Utils::Vector3d, double> local_system_CMS_velocity() {
  return boost::accumulate(
      cell_structure.local_particles(), std::pair<Utils::Vector3d, double>{},
      [](auto sum, const Particle &p) {
        if (not p.p.is_virtual) {
          return std::pair<Utils::Vector3d, double>{
              sum.first + p.p.mass * p.m.v, sum.second + p.p.mass};
        }
        return std::pair<Utils::Vector3d, double>{sum.first, sum.second};
      });
}

/* Remove the CMS velocity */
void local_galilei_transform(const Utils::Vector3d &cmsvel) {
  for (auto &p : cell_structure.local_particles()) {
    p.m.v -= cmsvel;
  }
}
