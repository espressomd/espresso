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

#ifdef VIRTUAL_SITES_RELATIVE

#include "cells.hpp"
#include "forces_inline.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "rotation.hpp"

#include <utils/math/quaternion.hpp>
#include <utils/math/sqr.hpp>
#include <utils/math/tensor_product.hpp>

namespace {
/**
 * @brief Orientation of the virtual site.
 * @param p_ref Reference particle for the virtual site.
 * @param vs_rel Parameters for the virtual site.
 * @return Orientation quaternion of the virtual site.
 */
Utils::Vector4d
orientation(Particle const *p_ref,
            const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  return multiply_quaternions(p_ref->r.quat, vs_rel.quat);
}

/**
 * @brief Vector pointing from the real particle
 *        to the virtual site.
 *
 * @return Relative distance.
 */
Utils::Vector3d connection_vector(
    Particle const *p_ref,
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  // Calculate the quaternion defining the orientation of the vector connecting
  // the virtual site and the real particle
  // This is obtained, by multiplying the quaternion representing the director
  // of the real particle with the quaternion of the virtual particle, which
  // specifies the relative orientation.
  auto const director =
      Utils::convert_quaternion_to_director(
          Utils::multiply_quaternions(p_ref->r.quat, vs_rel.rel_orientation))
          .normalize();

  return vs_rel.distance * director;
}

/**
 * @brief Position of the virtual site
 * @param p_ref Reference particle for the virtual site.
 * @param vs_rel Parameters for the virtual site.
 * @return Position of the virtual site.
 */
Utils::Vector3d
position(Particle const *p_ref,
         const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  return p_ref->r.p + connection_vector(p_ref, vs_rel);
}

/**
 * @brief Velocity of the virtual site
 * @param p_ref Reference particle for the virtual site.
 * @param vs_rel Parameters for the virtual site.
 * @return Velocity of the virtual site.
 */
Utils::Vector3d
velocity(const Particle *p_ref,
         const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  auto const d = connection_vector(p_ref, vs_rel);

  // Get omega of real particle in space-fixed frame
  auto const omega_space_frame =
      convert_vector_body_to_space(*p_ref, p_ref->m.omega);
  // Obtain velocity from v=v_real particle + omega_real_particle \times
  // director
  return vector_product(omega_space_frame, d) + p_ref->m.v;
}

/**
 * @brief Get reference particle.
 *
 * @param vs_rel Parameters to get the reference particle for.
 * @return Pointer to reference particle.
 */
Particle *get_reference_particle(
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  auto p_ref = cell_structure.get_local_particle(vs_rel.to_particle_id);
  if (!p_ref) {
    throw std::runtime_error("No real particle associated with virtual site.");
  }

  return p_ref;
}

/**
 * @brief Constraint force on the real particle.
 *
 *  Calculates the force exerted by the constraint on the
 *  reference particle.
 */
ParticleForce constraint_force(
    const ParticleForce &f, const Particle *p_ref,
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  return {f.f,
          vector_product(connection_vector(p_ref, vs_rel), f.f) + f.torque};
}

/**
 * @brief Constraint force to hold the particles at its prescribed position.
 *
 * @param f Force on the virtual site.
 * @param p_ref Reference particle.
 * @param vs_rel Parameters.
 * @return Constraint force.
 */
auto constraint_stress(
    const Utils::Vector3d &f, const Particle *p_ref,
    const ParticleProperties::VirtualSitesRelativeParameters &vs_rel) {
  /* The constraint force is minus the force on the particle, make it force
   * free. The counter force is translated by the connection vector to the
   * real particle, hence the virial stress is */
  return tensor_product(connection_vector(p_ref, vs_rel), -f);
}
} // namespace

void VirtualSitesRelative::update() const {
  cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
                               Cells::DATA_PART_MOMENTUM);

  for (auto &p : cell_structure.local_particles()) {
    if (!p.p.is_virtual)
      continue;

    const Particle *p_ref = get_reference_particle(p.p.vs_relative);

    auto const new_pos = position(p_ref, p.p.vs_relative);
    /* The shift has to respect periodic boundaries: if the reference
     * particles is not in the same image box, we potentially avoid to shift
     * to the other side of the box. */
    p.r.p += get_mi_vector(new_pos, p.r.p, box_geo);

    p.m.v = velocity(p_ref, p.p.vs_relative);

    if (get_have_quaternion())
      p.r.quat = orientation(p_ref, p.p.vs_relative);

    if ((p.r.p - p.l.p_old).norm2() > Utils::sqr(0.5 * skin))
      cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  } // namespace
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesRelative::back_transfer_forces_and_torques() const {
  cell_structure.ghosts_reduce_forces();

  init_forces_ghosts(cell_structure.ghost_particles());

  // Iterate over all the particles in the local cells
  for (auto &p : cell_structure.local_particles()) {
    // We only care about virtual particles
    if (p.p.is_virtual) {
      // First obtain the real particle responsible for this virtual particle:
      Particle *p_ref = get_reference_particle(p.p.vs_relative);

      // Add forces and torques
      p_ref->f += constraint_force(p.f, p_ref, p.p.vs_relative);
    }
  }
}

// Rigid body contribution to scalar pressure and stress tensor
Utils::Matrix<double, 3, 3> VirtualSitesRelative::stress_tensor() const {
  Utils::Matrix<double, 3, 3> stress_tensor = {};

  for (auto &p : cell_structure.local_particles()) {
    if (!p.p.is_virtual)
      continue;

    // First obtain the real particle responsible for this virtual particle:
    const Particle *p_ref = get_reference_particle(p.p.vs_relative);

    stress_tensor += constraint_stress(p.f.f, p_ref, p.p.vs_relative);
  }

  return stress_tensor;
}
#endif
