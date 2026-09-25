/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "ExclusionRadius.hpp"

#include "cells.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

static auto get_real_particle(boost::mpi::communicator const &comm, int p_id) {
  assert(p_id >= 0);
  auto &system = System::get_system();
  auto ptr = system.cell_structure->get_local_particle(p_id);
  if (ptr != nullptr and ptr->is_ghost()) {
    ptr = nullptr;
  }
  assert(boost::mpi::all_reduce(comm, static_cast<int>(ptr != nullptr),
                                std::plus<>()) == 1);
  return ptr;
}

void ExclusionRadius::set_exclusion_range(double range) {
  if (range < 0.) {
    throw std::domain_error("Invalid value for exclusion range");
  }
  exclusion_range = range;
  recalc_derived_parameters();
}

void ExclusionRadius::set_exclusion_radius_per_type(map_type const &map) {
  for (auto const &[type, exclusion_radius] : map) {
    if (exclusion_radius < 0.) {
      throw std::domain_error("Invalid exclusion radius for type " +
                              std::to_string(type) + ": radius " +
                              std::to_string(exclusion_radius));
    }
  }
  exclusion_radius_per_type = map;
  recalc_derived_parameters();
}

void ExclusionRadius::recalc_derived_parameters() {
  m_max_exclusion_range = exclusion_range;
  for (auto const &item : exclusion_radius_per_type) {
    auto const radius = item.second;
    m_max_exclusion_range = std::max(m_max_exclusion_range, 2. * radius);
  }
}

/**
 * Check if an inserted particle is too close to neighboring particles.
 */
bool ExclusionRadius::check_exclusion_range(int p_id, int p_type) {
  /* Check the exclusion radius of the inserted particle */
  if (exclusion_radius_per_type.contains(p_type)) {
    if (exclusion_radius_per_type[p_type] == 0.) {
      return false;
    }
  }

  auto p1_ptr = get_real_particle(m_comm, p_id);

  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto &cell_structure = *system.cell_structure;

  /* Test whether particle @p p2 lies inside the exclusion range of the
   * inserted particle located at @p p1_pos. */
  auto const is_inside_exclusion_range = [&](Utils::Vector3d const &p1_pos,
                                             Particle const &p2) {
    double excluded_distance;
    if (not exclusion_radius_per_type.contains(p_type) or
        not exclusion_radius_per_type.contains(p2.type())) {
      excluded_distance = exclusion_range;
    } else if (exclusion_radius_per_type[p2.type()] == 0.) {
      return false;
    } else {
      excluded_distance = exclusion_radius_per_type[p_type] +
                          exclusion_radius_per_type[p2.type()];
    }
    auto const d_min = box_geo.get_mi_vector(p2.pos(), p1_pos).norm();
    return d_min < excluded_distance;
  };

  if (neighbor_search_order_n) {
    /* Exhaustive O(N) search. The candidate id list is global (every rank
     * holds the full list), but each real particle is owned by exactly one
     * rank. We broadcast the position of the inserted particle to every rank,
     * let each rank distance-test the candidates it owns locally (real,
     * non-ghost copies, to avoid double-counting ghosts), and OR-reduce the
     * partial results. This keeps the check correct even when a candidate
     * sits within exclusion_range of the inserted particle but beyond the
     * ghost layer of the inserted particle's domain, where get_local_particle
     * would return nullptr. */
    auto all_ids = get_particle_ids_parallel();
    /* remove the inserted particle id */
    std::erase(all_ids, p_id);

    /* broadcast the position of the inserted particle from its owning rank */
    auto p1_pos = (p1_ptr != nullptr) ? p1_ptr->pos() : Utils::Vector3d{};
    auto const owner_rank =
        boost::mpi::all_reduce(m_comm, (p1_ptr != nullptr) ? m_comm.rank() : -1,
                               boost::mpi::maximum<int>());
    boost::mpi::broadcast(m_comm, p1_pos, owner_rank);

    bool local_touched = false;
    for (auto const p2_id : all_ids) {
      auto const p2_ptr = cell_structure.get_local_particle(p2_id);
      /* only the owning rank (real, non-ghost copy) tests each candidate */
      if (p2_ptr != nullptr and not p2_ptr->is_ghost() and
          is_inside_exclusion_range(p1_pos, *p2_ptr)) {
        local_touched = true;
        break;
      }
    }
    return boost::mpi::all_reduce(m_comm, local_touched, std::logical_or<>());
  }

  /* Neighbor search via the cell structure: the neighbor ids returned by
   * get_short_range_neighbors are local/ghost on the inserted particle's rank
   * and within the interaction range, so resolving them with
   * get_local_particle is correct here. */
  std::vector<int> particle_ids;
  {
    auto &mutable_system = System::get_system();
    mutable_system.on_observable_calc();
    auto const local_ids =
        get_short_range_neighbors(mutable_system, p_id, m_max_exclusion_range);
    assert(p1_ptr == nullptr or !!local_ids);
    if (local_ids) {
      particle_ids = std::move(*local_ids);
    }
  }

  bool local_touched = false;
  if (p1_ptr != nullptr) {
    auto const &p1 = *p1_ptr;
    /* Check if the inserted particle within the exclusion radius of any other
     * particle */
    for (auto const p2_id : particle_ids) {
      if (auto const p2_ptr = cell_structure.get_local_particle(p2_id)) {
        if (is_inside_exclusion_range(p1.pos(), *p2_ptr)) {
          local_touched = true;
          break;
        }
      }
    }
    if (m_comm.rank() != 0) {
      m_comm.send(0, 1, local_touched);
    }
  } else if (m_comm.rank() == 0) {
    m_comm.recv(boost::mpi::any_source, 1, local_touched);
  }
  boost::mpi::broadcast(m_comm, local_touched, 0);
  return local_touched;
}

bool ExclusionRadius::check_exclusion_range(int pid) {
  auto const *p = get_real_particle(m_comm, pid);
  assert(boost::mpi::all_reduce(m_comm, static_cast<int>(p != nullptr),
                                std::plus<>()) == 1);
  int type_local = -1;
  if (m_comm.rank() == 0) {
    if (p) {
      type_local = p->type();
    } else {
      m_comm.recv(boost::mpi::any_source, 42, type_local);
    }
  } else if (p) {
    m_comm.send(0, 42, p->type());
  }
  boost::mpi::broadcast(m_comm, type_local, 0);
  return check_exclusion_range(pid, type_local);
}
