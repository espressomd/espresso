/*
 * Copyright (C) 2010-2025 The ESPResSo project
 * Copyright (C) 2010,2011 Rudolf Weeber
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
// #ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

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

#include <boost/mpi/collectives/all_reduce.hpp>

#include <functional>
#include <unordered_map>

/**  @brief Store information about the center of mass  of the molecule*/
struct ComInfo {
double total_mass = 0.0;
Utils::Vector3d weighted_position = {0., 0., 0.};
};

static bool is_vs_com(Particle const &p) {
  return p.propagation() & PropagationMode::TRANS_VS_CENTER_OF_MASS;
}

void vs_com_update_particles(CellStructure &cell_structure,
                             BoxGeometry const &box_geo) {

    cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
                               Cells::DATA_PART_MOMENTUM);

    // Store virtual site center of mass particles
    // (mold_id: vs_com_id)
    std::unordered_map<int, int> vitual_site_id_for_mol_id;
    // Store com information for each molecule id
    // (mold_id: com_info)
    std::unordered_map<int, std::shared_ptr<ComInfo>> m_com_by_mol_id;

    cell_structure.for_each_local_particle([&](Particle &p) {

    if (is_vs_com(p)) { // get vs_com particle
        vitual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
    }
    else if (!p.is_virtual()) { // update m_com_by_mol_id
        if (m_com_by_mol_id.find(p.mol_id()) == m_com_by_mol_id.end()) {
            m_com_by_mol_id[p.mol_id()] = 
                std::make_shared<ComInfo>(); // m_com_by_mol_id initialization
        }
        m_com_by_mol_id[p.mol_id()]->total_mass += p.mass();
        auto const pos_unfolded =
            box_geo.unfolded_position(p.pos(), p.image_box());
        m_com_by_mol_id[p.mol_id()]->weighted_position += p.mass() * pos_unfolded;
    }
    });

    for (const auto &[mol_id, com_info] : m_com_by_mol_id) {
        auto com = com_info->weighted_position / com_info->total_mass;
        auto const vs_id = vitual_site_id_for_mol_id[mol_id];
        auto const vs_ptr = cell_structure.get_local_particle(vs_id);
        if (vs_ptr == nullptr) {
        continue;
        } else {
        auto folded_pos = com;
        auto image_box = Utils::Vector3i{};

        box_geo.fold_position(folded_pos, image_box);

        vs_ptr->image_box() = image_box;
        vs_ptr->mass() = com_info->total_mass;
        vs_ptr->pos() = folded_pos;
        }
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void vs_com_back_transfer_forces_and_torques(
    CellStructure &cell_structure) {

    cell_structure.ghosts_reduce_forces();
    init_forces_ghosts(cell_structure);
    // Store virtual site center of mass particles
    // (mold_id: vs_com_id)
    std::unordered_map<int, int> vitual_site_id_for_mol_id;
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) { // get vs_com particle
            vitual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
        }
    });

    // Iterate over all the particles in the local cells
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) return; // Check if particle is a virtual site center of mass

    auto const vs_id = vitual_site_id_for_mol_id.at(p.mol_id());
    auto vs_ptr = cell_structure.get_local_particle(vs_id);
    p.force() += (p.mass() / vs_ptr->mass()) * vs_ptr->force();

  });
}

// #endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS