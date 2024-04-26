/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"
#ifdef VIRTUAL_SITES_RELATIVE

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/quaternion.hpp>

#include <cassert>

/**
 * @brief Vector pointing from the real particle to the virtual site.
 *
 * @return Relative distance.
 */
static auto connection_vector(Particle const &p_ref, Particle const &p_vs) {
  // Calculate the quaternion defining the orientation of the vector connecting
  // the virtual site and the real particle
  // This is obtained, by multiplying the quaternion representing the director
  // of the real particle with the quaternion of the virtual particle, which
  // specifies the relative orientation.
  auto const director = Utils::convert_quaternion_to_director(
                            p_ref.quat() * p_vs.vs_relative().rel_orientation)
                            .normalize();

  return p_vs.vs_relative().distance * director;
}

/**
 * @brief Velocity of the virtual site
 * @param p_ref Reference particle for the virtual site.
 * @param p_vs Virtual site.
 * @return Velocity of the virtual site.
 */
static Utils::Vector3d velocity(Particle const &p_ref, Particle const &p_vs) {
  auto const d = connection_vector(p_ref, p_vs);

  // Get omega of real particle in space-fixed frame
  auto const omega_space_frame =
      convert_vector_body_to_space(p_ref, p_ref.omega());
  // Obtain velocity from v = v_real particle + omega_real_particle * director
  return vector_product(omega_space_frame, d) + p_ref.v();
}

/**
 * @brief Get real particle tracked by a virtual site.
 *
 * @param cell_structure Cell structure.
 * @param p Virtual site.
 * @return Pointer to real particle.
 */
static Particle *get_reference_particle(CellStructure &cell_structure,
                                        Particle const &p) {
  auto const &vs_rel = p.vs_relative();
  if (vs_rel.to_particle_id == -1) {
    runtimeErrorMsg() << "Particle with id " << p.id()
                      << " is a dangling virtual site";
    return nullptr;
  }
  auto p_ref_ptr = cell_structure.get_local_particle(vs_rel.to_particle_id);
  if (!p_ref_ptr) {
    runtimeErrorMsg() << "No real particle with id " << vs_rel.to_particle_id
                      << " for virtual site with id " << p.id();
  }
  return p_ref_ptr;
}

/**
 * @brief Constraint force to hold the particle at its prescribed position.
 *
 * @param p_ref Reference particle.
 * @param p_vs Virtual site.
 * @return Constraint force.
 */
static auto constraint_stress(Particle const &p_ref, Particle const &p_vs) {
  /* The constraint force is minus the force on the particle, make it force
   * free. The counter force is translated by the connection vector to the
   * real particle, hence the virial stress is */
  return tensor_product(-p_vs.force(), connection_vector(p_ref, p_vs));
}

static bool is_vs_relative_trans(Particle const &p) {
  return p.propagation() & PropagationMode::TRANS_VS_RELATIVE;
}
static bool is_vs_relative_rot(Particle const &p) {
  return p.propagation() & PropagationMode::ROT_VS_RELATIVE;
}

static bool is_vs_relative(Particle const &p) {
  return (is_vs_relative_trans(p) or is_vs_relative_rot(p));
}

void vs_relative_update_particles(CellStructure &cell_structure,
                                  BoxGeometry const &box_geo) {
  cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
                               Cells::DATA_PART_MOMENTUM);

  auto const particles = cell_structure.local_particles();
  for (auto &p : particles) {
    if (!is_vs_relative(p)) {
      continue;
    }

    auto const *p_ref_ptr = get_reference_particle(cell_structure, p);
    if (!p_ref_ptr)
      continue;

    auto const &p_ref = *p_ref_ptr;

    // position update
    if (is_vs_relative_trans(p)) {
      p.image_box() = p_ref.image_box();
      p.pos() = p_ref.pos() + connection_vector(p_ref, p);
      p.v() = velocity(p_ref, p);

      if (box_geo.type() == BoxType::LEES_EDWARDS) {
        auto push = LeesEdwards::Push(box_geo);
        push(p, -1); // includes a position fold
      } else {
        box_geo.fold_position(p.pos(), p.image_box());
      }
    }

    // Orientation update
    if (is_vs_relative_rot(p)) {
      p.quat() = p_ref.quat() * p.vs_relative().quat;
    }
  }

  if (cell_structure.check_resort_required()) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void vs_relative_back_transfer_forces_and_torques(
    CellStructure &cell_structure) {
  cell_structure.ghosts_reduce_forces();

  init_forces_ghosts(cell_structure.ghost_particles());

  // Iterate over all the particles in the local cells
  for (auto &p : cell_structure.local_particles()) {
    if (!is_vs_relative(p))
      continue;

    auto *p_ref_ptr = get_reference_particle(cell_structure, p);
    assert(p_ref_ptr != nullptr);

    auto &p_ref = *p_ref_ptr;
    if (is_vs_relative_trans(p)) {
      p_ref.force() += p.force();
      p_ref.torque() += vector_product(connection_vector(p_ref, p), p.force());
    }

    if (is_vs_relative_rot(p)) {
      p_ref.torque() += p.torque();
    }
  }
}

// Rigid body contribution to scalar pressure and pressure tensor
Utils::Matrix<double, 3, 3>
vs_relative_pressure_tensor(CellStructure const &cell_structure) {
  Utils::Matrix<double, 3, 3> pressure_tensor = {};

  for (auto const &p : cell_structure.local_particles()) {
    if (is_vs_relative_trans(p)) {
      if (auto const *p_ref_ptr = cell_structure.get_local_particle(
              p.vs_relative().to_particle_id)) {
        pressure_tensor += constraint_stress(*p_ref_ptr, p);
      }
    }
  }

  return pressure_tensor;
}
#endif // VIRTUAL_SITES_RELATIVE
