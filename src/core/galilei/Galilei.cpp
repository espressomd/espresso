/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "galilei/Galilei.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "system/System.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>

#include <utils/Vector.hpp>

#include <tuple>

void Galilei::kill_particle_motion(System::System &system, bool omega) const {
#ifndef ROTATION
  std::ignore = omega;
#endif
  for (auto &p : system.cell_structure->local_particles()) {
    p.v() = {};
#ifdef ROTATION
    if (omega) {
      p.omega() = {};
    }
#endif // ROTATION
  }
  system.on_particle_change();
}

void Galilei::kill_particle_forces(System::System &system, bool torque) const {
#ifndef ROTATION
  std::ignore = torque;
#endif
  for (auto &p : system.cell_structure->local_particles()) {
    p.force() = {};
#ifdef ROTATION
    if (torque) {
      p.torque() = {};
    }
#endif // ROTATION
  }
  system.on_particle_change();
}

Utils::Vector3d
Galilei::calc_system_cms_position(System::System const &system) const {
  auto const &box_geo = *system.box_geo;
  auto total_mass = 0.;
  auto cms_pos = Utils::Vector3d{};
  for (auto const &p : system.cell_structure->local_particles()) {
    if (not p.is_virtual()) {
      total_mass += p.mass();
      cms_pos += p.mass() * box_geo.unfolded_position(p.pos(), p.image_box());
    }
  }
  total_mass = boost::mpi::all_reduce(comm_cart, total_mass, std::plus<>());
  cms_pos = boost::mpi::all_reduce(comm_cart, cms_pos, std::plus<>());
  cms_pos /= total_mass;
  return cms_pos;
}

Utils::Vector3d
Galilei::calc_system_cms_velocity(System::System const &system) const {
  auto total_mass = 0.;
  auto cms_vel = Utils::Vector3d{};
  for (auto const &p : system.cell_structure->local_particles()) {
    if (not p.is_virtual()) {
      total_mass += p.mass();
      cms_vel += p.mass() * p.v();
    }
  }
  total_mass = boost::mpi::all_reduce(comm_cart, total_mass, std::plus<>());
  cms_vel = boost::mpi::all_reduce(comm_cart, cms_vel, std::plus<>());
  cms_vel /= total_mass;
  return cms_vel;
}

void Galilei::galilei_transform(System::System &system) const {
  auto const cms_vel = calc_system_cms_velocity(system);
  for (auto &p : system.cell_structure->local_particles()) {
    p.v() -= cms_vel;
  }
  system.on_particle_change();
}
