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

#include "galilei.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"

#include <boost/range/numeric.hpp>
#include <boost/serialization/utility.hpp>

#include <utils/Vector.hpp>

#include <algorithm>
#include <utility>

/** Stop particle motion by setting the velocity of each particle to zero. */
void local_kill_particle_motion(int omega, const ParticleRange &particles) {
  for (auto &p : particles) {
    if (omega) {
      p.m = {};
    } else {
      p.v() = {};
    }
  }
}

/** Set all the forces acting on the particles to zero */
void local_kill_particle_forces(int torque, const ParticleRange &particles) {
  for (auto &p : particles) {
    if (torque) {
      p.f = {};
    } else {
      p.force() = {};
    }
  }
}

/** Calculate the CMS of the system */
std::pair<Utils::Vector3d, double> local_system_CMS() {
  return boost::accumulate(
      cell_structure.local_particles(), std::pair<Utils::Vector3d, double>{},
      [](auto sum, const Particle &p) {
        if (not p.is_virtual()) {
          return std::pair<Utils::Vector3d, double>{
              sum.first + p.mass() * unfolded_position(p.pos(), p.image_box(),
                                                       box_geo.length()),
              sum.second + p.mass()};
        }
        return std::pair<Utils::Vector3d, double>{sum.first, sum.second};
      });
}

/** Calculate the CMS velocity of the system */
std::pair<Utils::Vector3d, double> local_system_CMS_velocity() {
  return boost::accumulate(
      cell_structure.local_particles(), std::pair<Utils::Vector3d, double>{},
      [](auto sum, const Particle &p) {
        if (not p.is_virtual()) {
          return std::pair<Utils::Vector3d, double>{
              sum.first + p.mass() * p.v(), sum.second + p.mass()};
        }
        return std::pair<Utils::Vector3d, double>{sum.first, sum.second};
      });
}

/** Remove the CMS velocity */
void local_galilei_transform(const Utils::Vector3d &cmsvel) {
  for (auto &p : cell_structure.local_particles()) {
    p.v() -= cmsvel;
  }
}

void mpi_kill_particle_motion_local(int rotation) {
  local_kill_particle_motion(rotation, cell_structure.local_particles());
  on_particle_change();
}

REGISTER_CALLBACK(mpi_kill_particle_motion_local)

void mpi_kill_particle_motion(int rotation) {
  mpi_call_all(mpi_kill_particle_motion_local, rotation);
}

void mpi_kill_particle_forces_local(int torque) {
  local_kill_particle_forces(torque, cell_structure.local_particles());
  on_particle_change();
}

REGISTER_CALLBACK(mpi_kill_particle_forces_local)

void mpi_kill_particle_forces(int torque) {
  mpi_call_all(mpi_kill_particle_forces_local, torque);
}

struct pair_sum {
  template <class T, class U>
  auto operator()(std::pair<T, U> l, std::pair<T, U> r) const {
    return std::pair<T, U>{l.first + r.first, l.second + r.second};
  }
};

Utils::Vector3d mpi_system_CMS() {
  auto const data =
      mpi_call(Communication::Result::reduction, pair_sum{}, local_system_CMS);
  return data.first / data.second;
}

REGISTER_CALLBACK_REDUCTION(local_system_CMS_velocity, pair_sum{})

Utils::Vector3d mpi_system_CMS_velocity() {
  auto const data = mpi_call(Communication::Result::reduction, pair_sum{},
                             local_system_CMS_velocity);
  return data.first / data.second;
}

REGISTER_CALLBACK_REDUCTION(local_system_CMS, pair_sum{})

void mpi_galilei_transform_local(Utils::Vector3d const &cmsvel) {
  local_galilei_transform(cmsvel);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_galilei_transform_local)

void mpi_galilei_transform() {
  auto const cmsvel = mpi_system_CMS_velocity();
  mpi_call_all(mpi_galilei_transform_local, cmsvel);
}
