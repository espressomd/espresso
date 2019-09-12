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
#include <utils/math/sqr.hpp>

#ifdef VIRTUAL_SITES_RELATIVE

#include "cells.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "rotation.hpp"

void VirtualSitesRelative::update(bool recalc_positions) const {
  // Ghost update logic
  if (n_nodes > 0) {
    if (recalc_positions or get_have_velocity()) {
      ghost_communicator(&cell_structure.update_ghost_pos_comm);
    }
  }
  for (auto &p : local_cells.particles()) {
    if (!p.p.is_virtual)
      continue;

    if (recalc_positions)
      update_pos(p);

    if (get_have_velocity())
      update_vel(p);

    if (get_have_quaternion())
      update_virtual_particle_quaternion(p);
  }
}

void VirtualSitesRelative::update_virtual_particle_quaternion(
    Particle &p) const {
  const Particle *p_real = local_particles[p.p.vs_relative.to_particle_id];
  if (!p_real) {
    throw std::runtime_error(
        "virtual_sites_relative.cpp - update_mol_pos_particle(): No real "
        "particle associated with virtual site.\n");
  }
  multiply_quaternions(p_real->r.quat, p.p.vs_relative.quat, p.r.quat);
#ifdef DIPOLES
  // When dipoles are enabled, update dipole moment
#endif
}

void VirtualSitesRelative::update_pos(Particle &p) const {
  // First obtain the real particle responsible for this virtual particle:
  // Find the 1st real particle in the topology for the virtual particle's
  // mol_id
  const Particle *p_real = local_particles[p.p.vs_relative.to_particle_id];
  // Check, if a real particle was found
  if (!p_real) {
    runtimeErrorMsg()
        << "virtual_sites_relative.cpp - update_mol_pos_particle(): No real "
           "particle associated with virtual site.\n";
    return;
  }

  // Calculate the quaternion defining the orientation of the vector connecting
  // the virtual site and the real particle
  // This is obtained, by multiplying the quaternion representing the director
  // of the real particle with the quaternion of the virtual particle, which
  // specifies the relative orientation.
  auto const director =
      convert_quat_to_director(
          multiply_quaternions(p_real->r.quat, p.p.vs_relative.rel_orientation))
          .normalize();

  auto const new_pos = p_real->r.p + director * p.p.vs_relative.distance;
  /* The shift has to respect periodic boundaries: if the reference particles
   * is not in the same image box, we potentially avoid to shift to the other
   * side of the box. */
  auto const shift = get_mi_vector(new_pos, p.r.p, box_geo);
  p.r.p += shift;

  if ((p.r.p - p.l.p_old).norm2() > Utils::sqr(0.5 * skin))
    set_resort_particles(Cells::RESORT_LOCAL);
}

void VirtualSitesRelative::update_vel(Particle &p) const {
  // First obtain the real particle responsible for this virtual particle:
  Particle *p_real = local_particles[p.p.vs_relative.to_particle_id];
  // Check, if a real particle was found
  if (!p_real) {
    runtimeErrorMsg()
        << "virtual_sites_relative.cpp - update_mol_pos_particle(): No real "
           "particle associated with virtual site.\n";
    return;
  }

  auto const d = get_mi_vector(p.r.p, p_real->r.p, box_geo);

  // Get omega of real particle in space-fixed frame
  Utils::Vector3d omega_space_frame =
      convert_vector_body_to_space(*p_real, p_real->m.omega);
  // Obtain velocity from v=v_real particle + omega_real_particle \times
  // director
  p.m.v = vector_product(omega_space_frame, d) + p_real->m.v;
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesRelative::back_transfer_forces_and_torques() const {
  ghost_communicator(&cell_structure.collect_ghost_force_comm);
  init_forces_ghosts(cell_structure.ghost_cells().particles());

  // Iterate over all the particles in the local cells
  for (auto &p : local_cells.particles()) {
    // We only care about virtual particles
    if (p.p.is_virtual) {
      // First obtain the real particle responsible for this virtual particle:
      Particle *p_real = local_particles[p.p.vs_relative.to_particle_id];

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
void VirtualSitesRelative::pressure_and_stress_tensor_contribution(
    double *pressure, double *stress_tensor) const {
  // Division by 3 volume is somewhere else. (pressure.cpp after all pressure
  // calculations) Iterate over all the particles in the local cells

  for (auto &p : local_cells.particles()) {
    if (!p.p.is_virtual)
      continue;

    update_pos(p);

    // First obtain the real particle responsible for this virtual particle:
    const Particle *p_real = local_particles[p.p.vs_relative.to_particle_id];

    // Get distance vector pointing from real to virtual particle, respecting
    // periodic boundary i
    // conditions
    auto const d = get_mi_vector(p_real->r.p, p.r.p, box_geo);

    // Stress tensor contribution
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        stress_tensor[k * 3 + l] += p.f.f[k] * d[l];

    // Pressure = 1/3 trace of stress tensor
    // but the 1/3 is applied somewhere else.
    *pressure += p.f.f * d;
  }
}

#endif
