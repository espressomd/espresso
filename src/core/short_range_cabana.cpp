/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#include "cell_system/CellStructure.hpp"

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#ifdef SHARED_MEMORY_PARALLELISM

#include <Cabana_Core.hpp>
#include "cabana_data.hpp"
#include "custom_verlet_list.hpp"
#include <cassert>
#include <unordered_set>
#include <utility>



template <class SliceDouble3, class SliceInt>
inline void write_particle(Particle const &p, std::unordered_map<int, int> const &id_to_index, SliceDouble3 &s_position, SliceDouble3 &s_force, SliceInt &s_id, SliceInt &s_type) {
  auto const pos = p.pos();
  auto const id = id_to_index.at(p.id());
  s_position(id, 0) = pos[0];
  s_position(id, 1) = pos[1];
  s_position(id, 2) = pos[2];
  s_id(id) = p.id();
  s_type(id) = p.type();
  s_force(id, 0) = 0.0;
  s_force(id, 1) = 0.0;
  s_force(id, 2) = 0.0;
}

template <class BondKernel,
          class VerletCriterion = detail::True>
void cabana_short_range(BondKernel bond_kernel,
                      CellStructure &cell_structure, double pair_cutoff,
                      double bond_cutoff, BoxGeometry &box_geo,
                      InteractionsNonBonded &nonbonded_ias,
                      ParticleRange particles, ParticleRange ghost_particles,
                      VerletCriterion const &verlet_criterion = {}) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif

  #ifdef CALIPER
  CALI_MARK_BEGIN("Espresso - Bond Kernel");
  #endif

  assert(cell_structure.get_resort_particles() == Cells::RESORT_NONE);

  if (bond_cutoff >= 0.) {
    cell_structure.bond_loop(bond_kernel);
  }

  #ifdef CALIPER
    CALI_MARK_END("Espresso - Bond Kernel");
  #endif

  // Cabana short range loop
  if (pair_cutoff > 0.) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ===================================================
    // Setup Cabana Variables
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Setup");
#endif
    // Dont know where to do this better
    using data_types = Cabana::MemberTypes<double[3], double[3], int, int, int>;
    using memory_space = Kokkos::SharedSpace;
    using execution_space = Kokkos::DefaultExecutionSpace;

    using ListAlgorithm = Cabana::HalfNeighborTag;
    using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

    const int vector_length = 8;
#ifdef CALIPER
    CALI_MARK_END("Cabana - Setup");
#endif

    // ===================================================
    // Count unique particles and create Index map
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Index map");
#endif
    std::unordered_map<int, int> id_to_index;
    int index = 0;

    bool const rebuild = cell_structure.get_rebuild_verlet_list();
    //bool const rebuild = true;

    CabanaData saved_data;

    // Load saved data if we do not have to rebuild
    if (!rebuild) {
      saved_data = cell_structure.get_cabana_data();
    }

    // If we have to rebuild, we need to count the particles and create a new map
    if (rebuild) {
      
      for (auto const& p : particles) {
        id_to_index[p.id()] = index;
        index++;
      }

      for (auto const& p : ghost_particles) {
        if (not id_to_index.contains(p.id())) {
          id_to_index[p.id()] = index;
          index++;
        }
      }
    } else {
      // If we do not rebuild we can use the saved map
      id_to_index = saved_data.get_id_to_index();
      index = id_to_index.size();
    }

    const int number_of_unique_particles = index;
#ifdef CALIPER
    CALI_MARK_END("Cabana - Index map");
#endif

    // ===================================================
    // Create and fill particle storage
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Fill particle storage");
#endif
    Cabana::AoSoA<data_types, memory_space, vector_length> particle_storage("particles", number_of_unique_particles);
    auto slice_position = Cabana::slice<0>(particle_storage);
    auto slice_force = Cabana::slice<1>(particle_storage);
    auto slice_id = Cabana::slice<2>(particle_storage);
    auto slice_type = Cabana::slice<3>(particle_storage);
    
    for (auto const& p : particles) {
      write_particle(p, id_to_index, slice_position, slice_force, slice_id, slice_type);
    }

    for (auto const& p : ghost_particles) {
      // if the ghost is not in the previous map, but mpi moved it to this rank?
      // it will not have neighbors because we did not rebuild the verlet list.
      if (not id_to_index.contains(p.id())) {
        continue;
      }
      write_particle(p, id_to_index, slice_position, slice_force, slice_id, slice_type);
    }
#ifdef CALIPER
    CALI_MARK_END("Cabana - Fill particle storage");
#endif

    // ===================================================
    // Get Verlet Pairs and Fill list
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Verlet List");
#endif
    ListType verlet_list;
    
    // Rebuild verlet list if needed
    if (rebuild) {

      verlet_list = ListType(slice_position, 0, slice_position.size(), 64);
      
      auto kernel = [&](Particle const &p1, Particle const &p2) {
        verlet_list.addNeighbor(id_to_index.at(p1.id()), id_to_index.at(p2.id()));
      };

      cell_structure.cabana_verlet_list_loop(kernel, verlet_criterion);
    } else {
      // Else use the saved verlet list
      verlet_list = saved_data.get_verlet_list();
    }

    // Save data for next iteration if we just rebuilt
    if (rebuild) {
      CabanaData new_data(verlet_list, id_to_index);
      cell_structure.set_cabana_data(std::make_unique<CabanaData>(new_data));
    }

    // fill customverletlist with pairs
    auto first_neighbor_kernel = KOKKOS_LAMBDA(const int i, const int j) {

        Utils::Vector3d const pi = {slice_position(i, 0), slice_position(i, 1), slice_position(i, 2)};
        Utils::Vector3d const pj = {slice_position(j, 0), slice_position(j, 1), slice_position(j, 2)};

        Utils::Vector3d const dist_vec = box_geo.get_mi_vector(pi, pj);

        auto const dist = dist_vec.norm();

        //if (dist > pair_cutoff) {
        //  return;
        //}

        IA_parameters const& ia_param = nonbonded_ias.get_ia_param(slice_type(i), slice_type(j));

        const ParticleForce kokkos_force = calc_central_radial_force(ia_param, dist_vec, dist);

        //if (kokkos_force.f[0] == 0.0 && kokkos_force.f[1] == 0.0 && kokkos_force.f[2] == 0.0) {
        //  return;
        //}

        Kokkos::atomic_add(&slice_force(i, 0), kokkos_force.f[0]);
        Kokkos::atomic_add(&slice_force(i, 1), kokkos_force.f[1]);
        Kokkos::atomic_add(&slice_force(i, 2), kokkos_force.f[2]);
        
        Kokkos::atomic_add(&slice_force(j, 0), -kokkos_force.f[0]);
        Kokkos::atomic_add(&slice_force(j, 1), -kokkos_force.f[1]);
        Kokkos::atomic_add(&slice_force(j, 2), -kokkos_force.f[2]);  
    };
#ifdef CALIPER
    CALI_MARK_END("Cabana - Verlet List");
#endif

    // ===================================================
    // Execute Kernel
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Execute Kernel");
#endif
    Kokkos::RangePolicy<execution_space> policy(0, particle_storage.size());


    // TODO: Add option to switch "SerialOpTag" Between "TeamOpTag"
    // Feels like TeamOpTag is faster, atleast for large particle numbers
    Cabana::neighbor_parallel_for(policy, first_neighbor_kernel, verlet_list,
                                    Cabana::FirstNeighborsTag(),
                                    Cabana::TeamOpTag(), "verlet_list");

    Kokkos::fence();

#ifdef CALIPER
    CALI_MARK_END("Cabana - Execute Kernel");
#endif

    // ===================================================
    // Add forces to particles
    // ===================================================
#ifdef CALIPER
    CALI_MARK_BEGIN("Cabana - Particle Forces");
#endif
    for (auto & p : particles) {
        auto const id = id_to_index.at(p.id());
        Utils::Vector3d f_vec{slice_force(id,0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
    }

    std::unordered_set<int> processed_ids;

    for (auto & p : ghost_particles) {
        int const pid = p.id();
        // Check if the particle has already been processed
        if (processed_ids.find(pid) != processed_ids.end()) {
          continue;
        }

        // Check if the ghost particle is in the map, i.e. was used during force calculation
        if (id_to_index.find(pid) == id_to_index.end()) {
          continue;
        }

        auto const id = id_to_index.at(pid);

        // Only add forces to ghost particles that are not as normal particles in the map,
        // as they have already been added to the force calculation
        if (id < particles.size()) {
          continue;
        }

        processed_ids.insert(pid);

        Utils::Vector3d f_vec{slice_force(id, 0), slice_force(id, 1), slice_force(id, 2)};
        
        ParticleForce f(f_vec);
        p.force_and_torque() += f;
    }
#ifdef CALIPER
    CALI_MARK_END("Cabana - Particle Forces");
#endif

  }
}

#endif