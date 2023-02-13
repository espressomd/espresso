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

#include "VirtualSitesCenterOfMass.hpp"

// #ifdef VIRTUAL_SITES_CENTER_OF_MASS

#include "Particle.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/quaternion.hpp>


void VirtualSitesCenterOfMass::update() const {
  struct ComInfo {
    double total_mass = 0.0;
    Utils::Vector3d weighted_position_sum = {0., 0., 0.};
  };

  std::map<int, ComInfo> com_by_mol_id;

  // com_by_mol_id initialization
  for (const auto &[mol_id, vs_id]: vitual_site_id_for_mol_id) {
    com_by_mol_id[mol_id] = ComInfo;
  }

  auto const particles = cell_structure.local_particles();

  for (auto &p : particles) {
    if (com_by_mol_id.find(p.mol_id) != com_by_mol_id.end()) {
      com_by_mol_id[p.mol_id].total_mass += p.mass();
      com_by_mol_id[p.mol_id].weighted_position_sum += p.mass() * p.pos();
      }
    }

  double com;
  int vs_id;
  auto const *vs_ptr;

  for (auto &[mol_id, com_info]: com_by_mol_id){
    com = com_info.weighted_position_sum / com_info.total_mass;
    vs_id =  vitual_site_id_for_mol_id[mol_id];
    vs_ptr = cell_structure.get_local_particle(vs_id);
    if (vs_ptr == nullptr){
      continue;
    } else {
      vs_ptr->pos() = com;
      vs_ptr->mass() = com_info.total_mass;

      auto folded_pos = vs_ptr->pos();
      auto image_box = Utils::Vector3i{};
      Utils::Vector3i n_shifts{};
      fold_position(folded_pos, image_box, box_geo);

      vs_ptr->pos() = folded_pos;
      
    }
  }  
}



  // cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
  //                              Cells::DATA_PART_MOMENTUM);
  
  // for (auto &p : particles) {
  //   auto const *p_ref_ptr = get_reference_particle(p);
  //   if (!p_ref_ptr)
  //     continue;

  //   auto const &p_ref = *p_ref_ptr;
  //   auto new_pos = p_ref.pos() + connection_vector(p_ref, p);
  //   /* The shift has to respect periodic boundaries: if the reference
  //    * particles is not in the same image box, we potentially avoid shifting
  //    * to the other side of the box. */
  //   auto shift = box_geo.get_mi_vector(new_pos, p.pos());
  //   p.pos() += shift;
  //   Utils::Vector3i image_shift{};
  //   fold_position(shift, image_shift, box_geo);
  //   p.image_box() = p_ref.image_box() - image_shift;

  //   p.v() = velocity(p_ref, p);

  //   if (box_geo.type() == BoxType::LEES_EDWARDS) {
  //     auto const &lebc = box_geo.lees_edwards_bc();
  //     auto const shear_dir = lebc.shear_direction;
  //     auto const shear_normal = lebc.shear_plane_normal;
  //     auto const le_vel = lebc.shear_velocity;
  //     Utils::Vector3i n_shifts{};
  //     fold_position(new_pos, n_shifts, box_geo);
  //     p.v()[shear_dir] -= n_shifts[shear_normal] * le_vel;
  //   }

  //   if (have_quaternions())
  //     p.quat() = p_ref.quat() * p.vs_relative().quat;
  // }

  // if (cell_structure.check_resort_required(particles, skin)) {
  //   cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  // }




// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesCenterOfMass::back_transfer_forces_and_torques() const {
  cell_structure.ghosts_reduce_forces();

  init_forces_ghosts(cell_structure.ghost_particles());

  // Iterate over all the particles in the local cells
  for (auto &p : cell_structure.local_particles()) {
    auto *p_ref_ptr = get_reference_particle(p);
    if (!p_ref_ptr)
      continue;

    // Add forces and torques
    auto &p_ref = *p_ref_ptr;
    p_ref.force() += p.force();
    p_ref.torque() +=
        vector_product(connection_vector(p_ref, p), p.force()) + p.torque();
  }
}

// Rigid body contribution to scalar pressure and pressure tensor
Utils::Matrix<double, 3, 3> VirtualSitesCenterOfMass::pressure_tensor() const {
  Utils::Matrix<double, 3, 3> pressure_tensor = {};

  for (auto const &p : cell_structure.local_particles()) {
    auto const *p_ref_ptr = get_reference_particle(p);
    if (p_ref_ptr) {
      pressure_tensor += constraint_stress(*p_ref_ptr, p);
    }
  }

  return pressure_tensor;
}
// #endif // VIRTUAL_SITES_CENTER_OF_MASS
