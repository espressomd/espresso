/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "VirtualSitesRelative.hpp"
#include "forces_inline.hpp"
#include <utils/math/quaternion.hpp>
#include <utils/math/sqr.hpp>

#ifdef VIRTUAL_SITES_RELATIVE

#include "cells.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "rotation.hpp"

namespace {
Utils::Vector4d updated_quaternion(
    Particle const *p_real,
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  return multiply_quaternions(p_real->r.quat, vs_rel.quat);
}

Utils::Vector3d
updated_pos(Particle const *p_real,
            const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  // Calculate the quaternion defining the orientation of the vector connecting
  // the virtual site and the real particle
  // This is obtained, by multiplying the quaternion representing the director
  // of the real particle with the quaternion of the virtual particle, which
  // specifies the relative orientation.
  auto const director =
      Utils::convert_quaternion_to_director(
          Utils::multiply_quaternions(p_real->r.quat, vs_rel.rel_orientation))
          .normalize();

  return p_real->r.p + director * vs_rel.distance;
}

Utils::Vector3d
updated_vel(const Particle *p_real,
            const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  auto const d =
      Utils::convert_quaternion_to_director(
          Utils::multiply_quaternions(p_real->r.quat, vs_rel.rel_orientation))
          .normalize() *
      vs_rel.distance;

  // Get omega of real particle in space-fixed frame
  auto const omega_space_frame =
      convert_vector_body_to_space(*p_real, p_real->m.omega);
  // Obtain velocity from v=v_real particle + omega_real_particle \times
  // director
  return vector_product(omega_space_frame, d) + p_real->m.v;
}

Particle *get_real_particle(
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  auto p_real = get_local_particle_data(vs_rel.to_particle_id);
  if (!p_real) {
    throw std::runtime_error("No real particle associated with virtual site.");
  }

  return p_real;
}

} // namespace

void VirtualSitesRelative::update(bool recalc_positions) const {
  // Ghost update logic
  if (n_nodes > 1) {
    auto const data_parts =
        (recalc_positions ? GHOSTTRANS_POSITION : GHOSTTRANS_NONE) |
        (get_have_velocity() ? (GHOSTTRANS_POSITION | GHOSTTRANS_MOMENTUM)
                             : GHOSTTRANS_NONE);

    ghost_communicator(&cell_structure.exchange_ghosts_comm, data_parts);
  }
  for (auto &p : cell_structure.local_cells().particles()) {
    if (!p.p.is_virtual)
      continue;

    const Particle *p_real = get_real_particle(p.p.vs_relative);

    if (recalc_positions) {
      auto const new_pos = updated_pos(p_real, p.p.vs_relative);
      /* The shift has to respect periodic boundaries: if the reference
       * particles is not in the same image box, we potentially avoid to shift
       * to the other side of the box. */
      p.r.p += get_mi_vector(new_pos, p.r.p, box_geo);

      if ((p.r.p - p.l.p_old).norm2() > Utils::sqr(0.5 * skin))
        set_resort_particles(Cells::RESORT_LOCAL);
    }

    if (get_have_velocity())
      p.m.v = updated_vel(p_real, p.p.vs_relative);

    if (get_have_quaternion())
      p.r.quat = updated_quaternion(p_real, p.p.vs_relative);
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesRelative::back_transfer_forces_and_torques() const {
  ghost_communicator(&cell_structure.collect_ghost_force_comm,
                     GHOSTTRANS_FORCE);
  init_forces_ghosts(cell_structure.ghost_cells().particles());

  // Iterate over all the particles in the local cells
  for (auto &p : cell_structure.local_cells().particles()) {
    // We only care about virtual particles
    if (p.p.is_virtual) {
      // First obtain the real particle responsible for this virtual particle:
      Particle *p_real = get_real_particle(p.p.vs_relative);

      // The rules for transferring forces are:
      // F_realParticle +=F_virtualParticle
      // T_realParticle +=f_realParticle \times
      // (r_virtualParticle-r_realParticle)

      // Add forces and torques
      p_real->f.torque +=
          vector_product(get_mi_vector(p.r.p, p_real->r.p, box_geo), p.f.f) +
          p.f.torque;
      p_real->f.f += p.f.f;
    }
  }
}

// Rigid body contribution to scalar pressure and stress tensor
Utils::Matrix<double, 3, 3> VirtualSitesRelative::stress_tensor() const {
  Utils::Matrix<double, 3, 3> stress_tensor = {};

  for (auto &p : cell_structure.local_cells().particles()) {
    if (!p.p.is_virtual)
      continue;

    // First obtain the real particle responsible for this virtual particle:
    const Particle *p_real = get_real_particle(p.p.vs_relative);

    // Get distance vector pointing from real to virtual particle, respecting
    // periodic boundary i
    // conditions
    auto const d = get_mi_vector(p_real->r.p, p.r.p, box_geo);

    stress_tensor += tensor_product(d, p.f.f);
  }

  return stress_tensor;
}
#endif
