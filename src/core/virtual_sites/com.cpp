/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

/**
 * @file
 * This file contains routine to handle virtual sites at the center of mass of
 * a bunch of other particles (say, a molecule). Forces acting on the center of
 * mass are distributed back onto the constituent particles.
 * The position/velocity/mass of the virtual site at center of mass
 * is calculated from the positions/velocities/masses of many particles.
 *
 * Virtual sites are like particles, but they will not be integrated.
 * Step performed for virtual sites:
 * - update virtual sites
 * - calculate forces
 * - distribute forces
 * - move non-virtual particles
 * - update virtual sites
 */

#include <config/config.hpp>

#ifdef ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/for_each_particle.hpp"
#include "cells.hpp"
#include "communication.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Center of mass information for a molecule.
 *
 * Hold the total mass and the mass-weighted position sum for a molecule,
 * used in center of mass calculations and MPI communication.
 */
struct ComInfo {
  /** @brief Sum of the masses of all constituent particles. */
  double total_mass = 0.;
  /** @brief Sum of (mass * position) for all constituent particles. */
  Utils::Vector3d weighted_position = {0., 0., 0.};

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, unsigned int const version) {
    ar & total_mass;
    ar & weighted_position;
  }
};

using ParticleId = decltype(ParticleProperties::identity);
using MoleculeId = decltype(ParticleProperties::mol_id);

static bool is_vs_com(Particle const &p) {
  return p.propagation() & PropagationMode::TRANS_VS_CENTER_OF_MASS;
}

std::optional<int> get_pid_for_vs_com(CellStructure &cell_structure,
                                      int mol_id) {
  auto constexpr parallel_execution_policy = false;
  std::vector<int> pids;
  cell_structure.for_each_local_particle(
      [&](Particle const &p) {
        if (is_vs_com(p) and p.mol_id() == mol_id) {
          pids.emplace_back(p.id());
        }
      },
      parallel_execution_policy);
  Utils::Mpi::gather_buffer(pids, comm_cart);
  std::optional<int> result = std::nullopt;
  if (not pids.empty()) {
    assert(pids.size() == 1ul);
    result = pids.front();
  }
  return result;
}

/**
 * @brief Synchronize a map across all MPI ranks.
 *
 * Gathers all key-value pairs from all ranks, broadcasts the full set to every
 * rank, and reconstructs the map so all processes have the same data.
 */
template <typename K, typename V>
void gather_buffer_map(std::unordered_map<K, V> &map_data) {
  std::vector<std::pair<K, V>> flat_map{map_data.begin(), map_data.end()};
  Utils::Mpi::gather_buffer(flat_map, comm_cart);
  boost::mpi::broadcast(comm_cart, flat_map, 0);
  map_data = {flat_map.begin(), flat_map.end()};
}

void vs_com_update_particles(CellStructure &cell_structure,
                             BoxGeometry const &box_geo) {

  auto constexpr parallel_execution_policy = false;
  // Update ghost positions before computing centers of mass
  cell_structure.ghosts_update(Cells::DATA_PART_POSITION);

  // Store virtual site center of mass particles
  // (mol_id: vs_com_id)
  std::unordered_map<ParticleId, MoleculeId> virtual_site_id_for_mol_id;
  // Store com information for each molecule id
  // (mol_id: com_info)
  std::unordered_map<MoleculeId, ComInfo> m_com_by_mol_id;

  cell_structure.for_each_local_particle(
      [&](Particle const &p) {
        if (is_vs_com(p)) {
          assert(not virtual_site_id_for_mol_id.contains(p.mol_id()));
          virtual_site_id_for_mol_id[p.mol_id()] = p.id();
        } else if (not p.is_virtual()) {
          if (not m_com_by_mol_id.contains(p.mol_id())) {
            m_com_by_mol_id[p.mol_id()] = ComInfo{};
          }
          auto const pos_unfolded =
              box_geo.unfolded_position(p.pos(), p.image_box());
          m_com_by_mol_id[p.mol_id()].total_mass += p.mass();
          m_com_by_mol_id[p.mol_id()].weighted_position +=
              p.mass() * pos_unfolded;
        }
      },
      parallel_execution_policy);

  // Reduction of m_com_by_mol_id across all processes
  // get a list of all molids that need to be communicated
  std::unordered_set<MoleculeId> local_mol_ids;
  for (auto const &mol_id : m_com_by_mol_id | std::views::keys) {
    local_mol_ids.insert(mol_id);
  }

  // Communicate the list of molids to all processes
  std::vector<std::unordered_set<MoleculeId>> global_mol_ids{};
  boost::mpi::all_gather(comm_cart, local_mol_ids, global_mol_ids);
  std::unordered_set<MoleculeId> unique_mol_ids{};
  for (auto const &mol_id_set : global_mol_ids) {
    for (auto const &mol_id : mol_id_set) {
      unique_mol_ids.insert(mol_id);
    }
  }
  std::vector<MoleculeId> flattened_mol_ids{unique_mol_ids.begin(),
                                            unique_mol_ids.end()};
  std::ranges::sort(flattened_mol_ids);

  // MPI-Allreduce for total mass and weighted position
  for (auto const mol_id : flattened_mol_ids) {
    double local_total_mass = 0.;
    Utils::Vector3d local_weighted_position = {0., 0., 0.};

    if (m_com_by_mol_id.contains(mol_id)) {
      local_total_mass = m_com_by_mol_id[mol_id].total_mass;
      local_weighted_position = m_com_by_mol_id[mol_id].weighted_position;
    }

    auto const total_mass =
        boost::mpi::all_reduce(comm_cart, local_total_mass, std::plus{});
    auto const weighted_position =
        boost::mpi::all_reduce(comm_cart, local_weighted_position, std::plus{});

    if (m_com_by_mol_id.contains(mol_id)) {
      m_com_by_mol_id[mol_id].total_mass = total_mass;
      m_com_by_mol_id[mol_id].weighted_position = weighted_position;
    }
  }

  // communicate m_com_by_mol_id across all processes
  gather_buffer_map(m_com_by_mol_id);

  // communicate the virtual_site_id_for_mol_id across all processes
  gather_buffer_map(virtual_site_id_for_mol_id);

  // Iterate over all the molecules and set the position and mass of the
  // associated virtual site com particle
  for (auto const &[mol_id, com_info] : m_com_by_mol_id) {
    if (not virtual_site_id_for_mol_id.contains(mol_id)) {
      continue;
    }
    auto const vs_id = virtual_site_id_for_mol_id[mol_id];
    // Mutable optional view -- the fields are written through below.
    if (auto vs_ptr = cell_structure.get_local_particle(vs_id)) {
      auto folded_pos = com_info.weighted_position / com_info.total_mass;
      auto image_box = Utils::Vector3i{};

      box_geo.fold_position(folded_pos, image_box);

      vs_ptr->image_box() = image_box;
      vs_ptr->mass() = com_info.total_mass;
      vs_ptr->pos() = folded_pos;
    }
  }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void vs_com_back_transfer_forces_and_torques(CellStructure &cell_structure) {

  auto constexpr parallel_execution_policy = false;
  cell_structure.ghosts_reduce_forces();
  cell_structure.ghosts_reset_forces();

  // Store forces for virtual site com particles
  // (vs_com_id: force)
  std::unordered_map<ParticleId, Utils::Vector3d> force_for_vs_id;
  // (vs_com_id: mass)
  std::unordered_map<ParticleId, double> mass_for_vs_id;

  // Store virtual site center of mass particles
  // (mol_id: vs_com_id)
  std::unordered_map<MoleculeId, ParticleId> virtual_site_id_for_mol_id;
  cell_structure.for_each_local_particle(
      [&](Particle const &p) {
        if (is_vs_com(p)) { // get vs_com particle
          assert(not virtual_site_id_for_mol_id.contains(p.mol_id()));
          virtual_site_id_for_mol_id[p.mol_id()] = p.id();
          force_for_vs_id[p.id()] = p.force();
          mass_for_vs_id[p.id()] = p.mass();
        }
      },
      parallel_execution_policy);

  // communicate virtual_site_id_for_mol_id across all processes
  gather_buffer_map(virtual_site_id_for_mol_id);
  // communicate force_for_vs_id, namely the force acting on the virtual site
  // com particles to all processes
  gather_buffer_map(force_for_vs_id);
  // communicate mass_for_vs_id, namely the mass of the virtual site com
  // particles to all processes
  gather_buffer_map(mass_for_vs_id);

  // Iterate over all the particles in the local cells
  cell_structure.for_each_local_particle(
      [&](Particle &p) {
        if (is_vs_com(p) or
            not virtual_site_id_for_mol_id.contains(p.mol_id())) {
          return;
        }
        auto const vs_id = virtual_site_id_for_mol_id.at(p.mol_id());
        p.force() +=
            (p.mass() / mass_for_vs_id.at(vs_id)) * force_for_vs_id.at(vs_id);
      },
      parallel_execution_policy);
}

#endif // ESPRESSO_VIRTUAL_SITES_CENTER_OF_MASS
