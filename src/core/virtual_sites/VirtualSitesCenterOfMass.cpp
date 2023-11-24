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
#include "communication.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/quaternion.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <functional>
#include <unordered_map>

void VirtualSitesCenterOfMass::update() {

  // com_by_mol_id initialization
  for (const auto &[mol_id, vs_id] : vitual_site_id_for_mol_id) {
    com_by_mol_id.emplace(std::make_pair(mol_id, std::make_shared<ComInfo>()));
  }

  auto const particles = cell_structure.local_particles();

  // Update com_by_mol_id
  for (const auto &p : particles) {
    if (com_by_mol_id.find(p.mol_id()) != com_by_mol_id.end()) {
      com_by_mol_id[p.mol_id()]->total_mass += p.mass();
      com_by_mol_id[p.mol_id()]->weighted_position_sum +=
          p.mass() * p.pos(); // Are these the unforlded positions?
    }
  }

  // Reduction operation
  for (auto &kv : com_by_mol_id) {
    auto const tot_mass =
        boost::mpi::all_reduce(comm_cart, kv.second->total_mass, std::plus());
    kv.second->total_mass = tot_mass;
    auto const weighted_position_sum = boost::mpi::all_reduce(
        comm_cart, kv.second->weighted_position_sum, std::plus());
    kv.second->weighted_position_sum = weighted_position_sum;
  }

  Utils::Vector3d com;
  int vs_id;

  for (auto &[mol_id, com_info] : com_by_mol_id) {
    com = com_info->weighted_position_sum /
          com_info->total_mass; // Is '/' definend in Utils::Vector3d?
    vs_id = vitual_site_id_for_mol_id[mol_id];
    auto vs_ptr = cell_structure.get_local_particle(
        vs_id); // When has the VS been defined?
    if (vs_ptr == nullptr) {
      continue;
    } else {
      vs_ptr->pos() = com;
      vs_ptr->mass() = com_info->total_mass;

      auto folded_pos = vs_ptr->pos();
      auto image_box = Utils::Vector3i{};
      Utils::Vector3i n_shifts{};
      fold_position(folded_pos, image_box, box_geo);

      vs_ptr->pos() = folded_pos;
    }
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles

void VirtualSitesCenterOfMass::back_transfer_forces() {

  // cell_structure.ghosts_reduce_forces();
  // init_forces_ghosts(cell_structure.ghost_particles());

  //int p_mol_id;
  //int vs_id;
  for (auto &p : cell_structure.local_particles()) {
    //p_mol_id = p.mol_id();
    //vs_id = vitual_site_id_for_mol_id[p.mol_id()];
    auto vs_ptr = cell_structure.get_local_particle(vitual_site_id_for_mol_id[p.mol_id()]);
    p.force() +=
        (p.mass() / com_by_mol_id[p.mol_id()]->total_mass) * vs_ptr->force();
  }
}

// #endif // VIRTUAL_SITES_CENTER_OF_MASS
