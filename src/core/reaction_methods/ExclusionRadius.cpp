/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "reaction_methods/ExclusionRadius.hpp"

#include "cells.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "particle_node.hpp"

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
  auto ptr = ::cell_structure.get_local_particle(p_id);
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
  if (exclusion_radius_per_type.count(p_type) != 0) {
    if (exclusion_radius_per_type[p_type] == 0.) {
      return false;
    }
  }

  auto p1_ptr = get_real_particle(m_comm, p_id);

  std::vector<int> particle_ids;
  if (neighbor_search_order_n) {
    auto all_ids = get_particle_ids_parallel();
    /* remove the inserted particle id */
    all_ids.erase(std::remove(all_ids.begin(), all_ids.end(), p_id),
                  all_ids.end());
    particle_ids = all_ids;
  } else {
    on_observable_calc();
    auto const local_ids =
        get_short_range_neighbors(p_id, m_max_exclusion_range);
    assert(p1_ptr == nullptr or !!local_ids);
    if (local_ids) {
      particle_ids = std::move(*local_ids);
    }
  }

  auto within_exclusion_range = false;
  if (p1_ptr != nullptr) {
    auto &p1 = *p1_ptr;

    /* Check if the inserted particle within any exclusion radius */
    for (auto const p2_id : particle_ids) {
      if (auto const p2_ptr = ::cell_structure.get_local_particle(p2_id)) {
        auto const &p2 = *p2_ptr;
        double excluded_distance;
        if (exclusion_radius_per_type.count(p_type) == 0 or
            exclusion_radius_per_type.count(p2.type()) == 0) {
          excluded_distance = exclusion_range;
        } else if (exclusion_radius_per_type[p2.type()] == 0.) {
          continue;
        } else {
          excluded_distance = exclusion_radius_per_type[p_type] +
                              exclusion_radius_per_type[p2.type()];
        }

        auto const d_min = ::box_geo.get_mi_vector(p2.pos(), p1.pos()).norm();

        if (d_min < excluded_distance) {
          within_exclusion_range = true;
          break;
        }
      }
    }
    if (m_comm.rank() != 0) {
      m_comm.send(0, 1, within_exclusion_range);
    }
  } else if (m_comm.rank() == 0) {
    m_comm.recv(boost::mpi::any_source, 1, within_exclusion_range);
  }
  boost::mpi::broadcast(m_comm, within_exclusion_range, 0);
  return within_exclusion_range;
}

bool ExclusionRadius::check_exclusion_range(int pid) {
  int type_local = 0;
  if (auto p = get_real_particle(m_comm, pid)) {
    type_local = p->type();
  }
  auto const type =
      boost::mpi::all_reduce(m_comm, type_local, std::plus<int>());
  return check_exclusion_range(pid, type);
}
