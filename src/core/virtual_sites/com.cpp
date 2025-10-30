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
#include "communication.hpp"
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
#include <boost/mpi/collectives/all_gather.hpp>

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

    cell_structure.ghosts_update(Cells::DATA_PART_POSITION);

    // Store virtual site center of mass particles
    // (mol_id: vs_com_id)
    std::unordered_map<int, int> virtual_site_id_for_mol_id;
    // Store com information for each molecule id
    // (mol_id: com_info)
    std::unordered_map<int, std::shared_ptr<ComInfo>> m_com_by_mol_id;

    cell_structure.for_each_local_particle([&](Particle &p) {

        // std::cout << "Particle ID: " << p.id()
        //           << " Pos: (" << p.pos()[0] << ", "
        //           << p.pos()[1] << ", "
        //           << p.pos()[2] << ")\n";

    if (is_vs_com(p)) { // get vs_com particle
        virtual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
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

    // --------------------------------------------------------------------------------------------------

    auto const rank = comm_cart.rank();

    std::cout << "Rank: " << rank 
              << " Number of unique mol_ids: " << m_com_by_mol_id.size()
              << std::endl;

    // log the rank and each entry of m_com_by_mol_id
    for (const auto &[mol_id, com_info] : m_com_by_mol_id) {
        std::cout << "Rank: " << rank
                  <<  " Mol ID: " << mol_id 
                  << " Total Mass: " << com_info->total_mass 
                  << " Weighted Position: (" 
                  << com_info->weighted_position[0] << ", "
                  << com_info->weighted_position[1] << ", "
                  << com_info->weighted_position[2] << ")"
                  << std::endl;
    }

    // get a list of all molids that need to be communicated
    std::vector<int> local_mol_ids;
    for (const auto &[mol_id, _] : m_com_by_mol_id) {
        local_mol_ids.emplace_back(mol_id);
    }
    std::sort(local_mol_ids.begin(), local_mol_ids.end());
    local_mol_ids.erase(std::unique(local_mol_ids.begin(), local_mol_ids.end()), local_mol_ids.end());

    // Communicate the list of molids to all processes
    std::vector<std::vector<int>> all_mol_ids;
    boost::mpi::all_gather(comm_cart, local_mol_ids, all_mol_ids);
    // flatten the list of lists into a single list
    std::vector<int> flattened_mol_ids{};
    for (const auto &mol_id_list : all_mol_ids) {
        flattened_mol_ids.insert(flattened_mol_ids.end(), mol_id_list.begin(), mol_id_list.end());
    }
    std::sort(flattened_mol_ids.begin(), flattened_mol_ids.end());
    flattened_mol_ids.erase(std::unique(flattened_mol_ids.begin(), flattened_mol_ids.end()), flattened_mol_ids.end());

    // MPI-Allreduce for total mass and weighted position
    for (const auto mol_id : flattened_mol_ids) {
        double local_total_mass = 0.0;
        Utils::Vector3d local_weighted_position = {0., 0., 0.};

        if (m_com_by_mol_id.find(mol_id) != m_com_by_mol_id.end()) {
            local_total_mass = m_com_by_mol_id[mol_id]->total_mass;
            local_weighted_position = m_com_by_mol_id[mol_id]->weighted_position;
        }

        double global_total_mass = 0.0;
        Utils::Vector3d global_weighted_position = {0., 0., 0.};

        boost::mpi::all_reduce(comm_cart, local_total_mass, global_total_mass, std::plus());
        boost::mpi::all_reduce(comm_cart, local_weighted_position, global_weighted_position, std::plus());

        if (m_com_by_mol_id.find(mol_id) != m_com_by_mol_id.end()) {
            m_com_by_mol_id[mol_id]->total_mass = global_total_mass;
            m_com_by_mol_id[mol_id]->weighted_position = global_weighted_position;
        }
    }

    std::cout << "Rank: " << rank 
              << " After Allreduce Number of unique mol_ids: " << m_com_by_mol_id.size()
              << std::endl;
    
    for (const auto &[mol_id, com_info] : m_com_by_mol_id) {
        std::cout << "Rank: " << rank
                  <<  " Mol ID: " << mol_id 
                  << " Total Mass: " << com_info->total_mass 
                  << " Weighted Position: (" 
                  << com_info->weighted_position[0] << ", "
                  << com_info->weighted_position[1] << ", "
                  << com_info->weighted_position[2] << ")"
                  << std::endl;
    }

    // --------------------------------------------------------------------------------------------------

    // for (auto &[mol_id, com_info] : m_com_by_mol_id) {
    //     auto const tot_mass =
    //         boost::mpi::all_reduce(comm_cart, com_info->total_mass, std::plus());
    //     com_info->total_mass = tot_mass;
    //     auto const weighted_position = boost::mpi::all_reduce(
    //         comm_cart, com_info->weighted_position, std::plus());
    //     com_info->weighted_position = weighted_position;
    // }


    // auto const particles = cell_structure.local_particles();
    // for (auto &p : particles) {
    //     if (is_vs_com(p)) { // get vs_com particle
    //         virtual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
    //     }
    //     else if (!p.is_virtual()) { // update m_com_by_mol_id
    //         if (m_com_by_mol_id.find(p.mol_id()) == m_com_by_mol_id.end()) {
    //             m_com_by_mol_id[p.mol_id()] = 
    //                 std::make_shared<ComInfo>(); // m_com_by_mol_id initialization
    //         }
    //         m_com_by_mol_id[p.mol_id()]->total_mass += p.mass();
    //         auto const pos_unfolded =
    //             box_geo.unfolded_position(p.pos(), p.image_box());
    //         m_com_by_mol_id[p.mol_id()]->weighted_position += p.mass() * pos_unfolded;
    //     }
    // }

    // // Reduction operation
    // for (auto &kv : m_com_by_mol_id) {
    //     auto const tot_mass =
    //         boost::mpi::all_reduce(comm_cart, kv.second->total_mass, std::plus());
    //     kv.second->total_mass = tot_mass;
    //     auto const weighted_position = boost::mpi::all_reduce(
    //         comm_cart, kv.second->weighted_position, std::plus());
    //     kv.second->weighted_position = weighted_position;
    // }


    // Iterate over all the molecules and set the position and mass of the
    // associated virtual site com particle
    for (const auto &[mol_id, com_info] : m_com_by_mol_id) {
        if (virtual_site_id_for_mol_id.find(mol_id) == 
        virtual_site_id_for_mol_id.end()) {
            continue; // No virtual site for this molecule id
        }
        auto com = com_info->weighted_position / com_info->total_mass;
        auto const vs_id = virtual_site_id_for_mol_id[mol_id];
        auto const vs_ptr = cell_structure.get_local_particle(vs_id);
        if (vs_ptr == nullptr) {
        continue;
        }

        auto folded_pos = com;
        auto image_box = Utils::Vector3i{};

        box_geo.fold_position(folded_pos, image_box);

        vs_ptr->image_box() = image_box;
        vs_ptr->mass() = com_info->total_mass;
        vs_ptr->pos() = folded_pos;

        std::cout << "Rank: " << rank
                  << " Mol ID: " << mol_id 
                  << " COM Pos: (" 
                  << com[0] << ", "
                  << com[1] << ", "
                  << com[2] << ")"
                  << " Folded Pos: (" 
                  << folded_pos[0] << ", "
                  << folded_pos[1] << ", "
                  << folded_pos[2] << ")"
                  << " Total Mass: " << com_info->total_mass 
                  << std::endl;
    }

    if (cell_structure.check_resort_required()) {
        cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
    }
    cell_structure.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES);

}



// struct COM_Output {
//     Utils::Vector3d position;
//     double mass;

//     COM_Output() : position{0., 0., 0.}, mass{0.} {}
//     COM_Output(Utils::Vector3d const &pos, double m) : position{pos}, mass{m} {}
//     // COM_Output(COM_Output const &other) = default;
//     // COM_Output &operator=(COM_Output const &other) = default;
//     // COM_Output(COM_Output &&other) noexcept = default;
//     // COM_Output &operator=(COM_Output &&other) noexcept = default;
//     // ~COM_Output() = default;
// };

// COM_Output center_of_mass(CellStructure &cell_structure,
//                           BoxGeometry const &box_geo, 
//                           int mol_id) {

//   Utils::Vector3d local_com{};
//   double local_mass = 0.;

//   for (auto const &p : cell_structure.local_particles()) {
//     if ((p.mol_id() == mol_id or mol_id == -1) and not p.is_virtual()) {
//       local_com += box_geo.unfolded_position(p.pos(), p.image_box()) * p.mass();
//       local_mass += p.mass();
//     }
//   }
//   Utils::Vector3d com{};
//   double mass = 1.; // placeholder value to avoid division by zero
//   std::cout << "Rank: " << comm_cart.rank() 
//             << "Performing reduction for mol_id: " << mol_id
//             << std::endl;
//   boost::mpi::reduce(comm_cart, local_com, com, std::plus<>(), 0);
//   boost::mpi::reduce(comm_cart, local_mass, mass, std::plus<>(), 0);
//   return COM_Output{com / mass, mass};
// }

// void vs_com_update_particles(CellStructure &cell_structure,
//                              BoxGeometry const &box_geo) {

//     cell_structure.ghosts_update(Cells::DATA_PART_POSITION |
//                                Cells::DATA_PART_MOMENTUM);

//     // Store virtual site center of mass particles
//     // (mol_id: vs_com_id)
//     std::unordered_map<int, int> virtual_site_id_for_mol_id;
//     // Store com information for each molecule id
//     // (mol_id: com_info)
//     std::unordered_map<int, std::shared_ptr<ComInfo>> m_com_by_mol_id;

//     cell_structure.for_each_local_particle([&](Particle &p) {

//     if (is_vs_com(p)) { // get vs_com particle
//         virtual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
//     }
//     else if (!p.is_virtual()) { // update m_com_by_mol_id
//         if (m_com_by_mol_id.find(p.mol_id()) == m_com_by_mol_id.end()) {
//             m_com_by_mol_id[p.mol_id()] = 
//                 std::make_shared<ComInfo>(); // m_com_by_mol_id initialization
//         }
//     }
//     });

//     for (const auto &[mol_id, _] : m_com_by_mol_id) {
//         auto com_info = center_of_mass(cell_structure, box_geo, mol_id);
//         m_com_by_mol_id[mol_id]->total_mass = com_info.mass;
//         m_com_by_mol_id[mol_id]->weighted_position = com_info.position * com_info.mass;
//     }

//     // --------------------------------------------------------------------------------------------------

//     // Iterate over all the molecules and set the position and mass of the
//     // associated virtual site com particle
//     for (const auto &[mol_id, com_info] : m_com_by_mol_id) {
//         if (virtual_site_id_for_mol_id.find(mol_id) == 
//         virtual_site_id_for_mol_id.end()) {
//             continue; // No virtual site for this molecule id
//         }
//         auto com = com_info->weighted_position / com_info->total_mass;
//         auto const vs_id = virtual_site_id_for_mol_id[mol_id];
//         auto const vs_ptr = cell_structure.get_local_particle(vs_id);
//         if (vs_ptr == nullptr) {
//         continue;
//         } else {
//         auto folded_pos = com;
//         auto image_box = Utils::Vector3i{};

//         box_geo.fold_position(folded_pos, image_box);

//         vs_ptr->image_box() = image_box;
//         vs_ptr->mass() = com_info->total_mass;
//         vs_ptr->pos() = folded_pos;
//         }
//   }
// }

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void vs_com_back_transfer_forces_and_torques(
    CellStructure &cell_structure) {

    cell_structure.ghosts_reduce_forces();
    init_forces_ghosts(cell_structure);
    // Store virtual site center of mass particles
    // (mold_id: vs_com_id)
    std::unordered_map<int, int> virtual_site_id_for_mol_id;
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) { // get vs_com particle
            virtual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
        }
    });

    // Iterate over all the particles in the local cells
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) return; // Check if particle is a virtual site center of mass
    
        if (virtual_site_id_for_mol_id.find(p.mol_id()) == 
            virtual_site_id_for_mol_id.end()) {
            return; // No virtual site for this molecule id
        }
        auto const vs_id = virtual_site_id_for_mol_id.at(p.mol_id());
        auto vs_ptr = cell_structure.get_local_particle(vs_id);
        p.force() += (p.mass() / vs_ptr->mass()) * vs_ptr->force();
  });

  cell_structure.ghosts_update(Cells::DATA_PART_FORCE);
}

// #endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS