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
#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

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
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/shared_ptr.hpp>

#include <functional>
#include <unordered_map>

/**
 * @brief Stores center of mass information for a molecule.
 *
 * Holds the total mass and the mass-weighted position sum for a molecule,
 * used in center of mass calculations and MPI communication.
 *
 * Members:
 * - total_mass: The sum of the masses of all constituent particles.
 * - weighted_position: The sum of (mass * position) for all constituent particles.
 *
 * Supports Boost.Serialization for MPI communication.
 */
struct ComInfo {
double total_mass = 0.0;
Utils::Vector3d weighted_position = {0., 0., 0.};

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & total_mass;
        ar & weighted_position;
    }
};

static bool is_vs_com(Particle const &p) {
  return p.propagation() & PropagationMode::TRANS_VS_CENTER_OF_MASS;
}

/**
 * @brief Synchronize an unordered_map across all MPI ranks.
 *
 * Gathers all key-value pairs from all ranks, broadcasts the full set to every rank,
 * and reconstructs the map so all processes have the same data.
 *
 * @tparam T1 Key type (must be serializable by Boost.MPI).
 * @tparam T2 Value type (must be serializable by Boost.MPI).
 * @param map_data The map to synchronize (input/output).
 * @param comm_cart The MPI communicator.
 *
 * Complexity: O(N) per rank, where N is the total number of key-value pairs across all ranks.
 *
 * All ranks must call this collectively. Requires Boost.MPI and compatible serialization for T1, T2.
 */
template <typename T1, typename T2>
void communicate_map(std::unordered_map<T1, T2> &map_data, boost::mpi::communicator const &comm_cart) {
    // Gather all keys and values from all processes
    std::vector<T1> keys;
    std::vector<T2> values;
    for (const auto &[key, value] : map_data) {
        keys.push_back(key);
        values.push_back(value);
    }

    Utils::Mpi::gather_buffer(keys, comm_cart);
    Utils::Mpi::gather_buffer(values, comm_cart);
    boost::mpi::broadcast(comm_cart, keys, 0);
    boost::mpi::broadcast(comm_cart, values, 0);

    // Clear the original map and reconstruct it with gathered data
    map_data.clear();
    for (size_t i = 0; i < keys.size(); ++i) {
        map_data[keys[i]] = values[i];
    }   
}

void vs_com_update_particles(CellStructure &cell_structure,
                             BoxGeometry const &box_geo) {

    // Update ghost positions before computing centers of mass
    cell_structure.ghosts_update(Cells::DATA_PART_POSITION);

    // Store virtual site center of mass particles
    // (mol_id: vs_com_id)
    std::unordered_map<int, int> virtual_site_id_for_mol_id;
    // Store com information for each molecule id
    // (mol_id: com_info)
    std::unordered_map<int, std::shared_ptr<ComInfo>> m_com_by_mol_id;

    cell_structure.for_each_local_particle([&](Particle &p) {

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

    // Reduction of m_com_by_mol_id across all processes
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

    // communicate m_com_by_mol_id across all processes
    communicate_map(m_com_by_mol_id, comm_cart);

    // communicate the virtual_site_id_for_mol_id across all processes
    communicate_map(virtual_site_id_for_mol_id, comm_cart);

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

    }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void vs_com_back_transfer_forces_and_torques(
    CellStructure &cell_structure) {

    cell_structure.ghosts_reduce_forces();
    init_forces_ghosts(cell_structure);

    // Store forces for virtual site com particles
    // (vs_com_id: force)
    std::unordered_map<int, Utils::Vector3d> force_for_vs_id;

    // Store virtual site center of mass particles
    // (mold_id: vs_com_id)
    std::unordered_map<int, int> virtual_site_id_for_mol_id;
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) { // get vs_com particle
            virtual_site_id_for_mol_id[p.vs_com().to_molecule_id] = p.id();
            force_for_vs_id[p.id()] = p.force();
        }
    });
    
    // gather and broadcast virtual_site_id_for_mol_id across all processes
    std::vector<std::vector<int>> tmp;
    for (const auto &[mol_id, vs_id] : virtual_site_id_for_mol_id) {
        tmp.emplace_back(std::vector<int>{mol_id, vs_id});
    }
    Utils::Mpi::gather_buffer(tmp, comm_cart);
    boost::mpi::broadcast(comm_cart, tmp, 0);
    std::unordered_map<int, int> all_tmp;
    for (const auto &inner_vec : tmp) {
        all_tmp[inner_vec[0]] = inner_vec[1];
    }
    std::unordered_map<int, int>(all_tmp.begin(), all_tmp.end()).swap(virtual_site_id_for_mol_id);

    // communicate force_for_vs_id, namely the force acting on the virtual site com particles to all processes 
    communicate_map(force_for_vs_id, comm_cart);

    // Iterate over all the particles in the local cells
    cell_structure.for_each_local_particle([&](Particle &p) {
        if (is_vs_com(p)) return; // Check if particle is a virtual site center of mass
    
        if (virtual_site_id_for_mol_id.find(p.mol_id()) == 
            virtual_site_id_for_mol_id.end()) {
            return; // No virtual site for this molecule id
        }
        auto const vs_id = virtual_site_id_for_mol_id.at(p.mol_id());
        auto vs_ptr = cell_structure.get_local_particle(vs_id);
        p.force() += (p.mass() / vs_ptr->mass()) * force_for_vs_id.at(vs_id);

  });
}

#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS