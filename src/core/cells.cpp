/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *
 *  This file contains functions for the cell system.
 *
 *  Implementation of cells.hpp.
 */

#include "cells.hpp"

#include "cell_system/Cell.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cell_system/HybridDecomposition.hpp"

#include "Particle.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/serialization/set.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

/** Type of cell structure in use */
CellStructure cell_structure{box_geo};

/**
 * @brief Get pairs of particles that are closer than a distance and fulfill a
 * filter criterion.
 *
 * It uses link_cell to get pairs out of the cellsystem
 * by a simple distance criterion and applies the filter on both particles.
 *
 * Pairs are sorted so that first.id < second.id
 */
template <class Filter>
std::vector<std::pair<int, int>> get_pairs_filtered(double const distance,
                                                    Filter filter) {
  std::vector<std::pair<int, int>> ret;
  auto const cutoff2 = Utils::sqr(distance);
  auto const pair_kernel = [cutoff2, &filter, &ret](Particle const &p1,
                                                    Particle const &p2,
                                                    Distance const &d) {
    if (d.dist2 < cutoff2 and filter(p1) and filter(p2))
      ret.emplace_back(p1.id(), p2.id());
  };

  cell_structure.non_bonded_loop(pair_kernel);

  /* Sort pairs */
  for (auto &pair : ret) {
    if (pair.first > pair.second)
      std::swap(pair.first, pair.second);
  }

  return ret;
}

namespace detail {
static auto get_max_neighbor_search_range() {
  return *boost::min_element(cell_structure.max_range());
}
static void search_distance_sanity_check_max_range(double const distance) {
  /* get_pairs_filtered() finds pairs via the non_bonded_loop. The maximum
   * finding range is therefore limited by the decomposition that is used.
   */
  auto const max_range = get_max_neighbor_search_range();
  if (distance > max_range) {
    throw std::domain_error("pair search distance " + std::to_string(distance) +
                            " bigger than the decomposition range " +
                            std::to_string(max_range));
  }
}
static void search_distance_sanity_check_cell_structure(double const distance) {
  if (cell_structure.decomposition_type() ==
      CellStructureType::CELL_STRUCTURE_HYBRID) {
    throw std::runtime_error("Cannot search for neighbors in the hybrid "
                             "decomposition cell system");
  }
}
static void search_neighbors_sanity_checks(double const distance) {
  search_distance_sanity_check_max_range(distance);
  search_distance_sanity_check_cell_structure(distance);
}
} // namespace detail

boost::optional<std::vector<int>>
get_short_range_neighbors(int const pid, double const distance) {
  detail::search_neighbors_sanity_checks(distance);
  std::vector<int> ret;
  auto const cutoff2 = Utils::sqr(distance);
  auto const kernel = [cutoff2, &ret](Particle const &, Particle const &p2,
                                      Utils::Vector3d const &vec) {
    if (vec.norm2() < cutoff2) {
      ret.emplace_back(p2.id());
    }
  };
  auto const p = ::cell_structure.get_local_particle(pid);
  if (p and not p->is_ghost()) {
    ::cell_structure.run_on_particle_short_range_neighbors(*p, kernel);
    return {ret};
  }
  return {};
}

/**
 * @brief Get pointers to all interacting neighbors of a central particle.
 */
static auto get_interacting_neighbors(Particle const &p) {
  auto const distance = *boost::min_element(::cell_structure.max_range());
  detail::search_neighbors_sanity_checks(distance);
  std::vector<Particle const *> ret;
  auto const cutoff2 = Utils::sqr(distance);
  auto const kernel = [cutoff2, &ret](Particle const &, Particle const &p2,
                                      Utils::Vector3d const &vec) {
    if (vec.norm2() < cutoff2) {
      ret.emplace_back(&p2);
    }
  };
  ::cell_structure.run_on_particle_short_range_neighbors(p, kernel);
  return ret;
}

std::vector<std::pair<int, int>> get_pairs(double const distance) {
  detail::search_neighbors_sanity_checks(distance);
  return get_pairs_filtered(distance, [](Particle const &) { return true; });
}

std::vector<std::pair<int, int>>
get_pairs_of_types(double const distance, std::vector<int> const &types) {
  detail::search_neighbors_sanity_checks(distance);
  return get_pairs_filtered(distance, [types](Particle const &p) {
    return std::any_of(types.begin(), types.end(),
                       // NOLINTNEXTLINE(bugprone-exception-escape)
                       [p](int const type) { return p.type() == type; });
  });
}

std::vector<PairInfo> non_bonded_loop_trace(int const rank) {
  std::vector<PairInfo> pairs;
  auto const pair_kernel = [&pairs, rank](Particle const &p1,
                                          Particle const &p2,
                                          Distance const &d) {
    pairs.emplace_back(p1.id(), p2.id(), p1.pos(), p2.pos(), d.vec21, rank);
  };
  cell_structure.non_bonded_loop(pair_kernel);
  return pairs;
}

std::vector<NeighborPIDs> get_neighbor_pids() {
  std::vector<NeighborPIDs> ret;
  auto kernel = [&ret](Particle const &p,
                       std::vector<Particle const *> const &neighbors) {
    std::vector<int> neighbor_pids;
    neighbor_pids.reserve(neighbors.size());
    for (auto const &neighbor : neighbors) {
      neighbor_pids.emplace_back(neighbor->id());
    }
    ret.emplace_back(p.id(), neighbor_pids);
  };
  for (auto const &p : ::cell_structure.local_particles()) {
    kernel(p, get_interacting_neighbors(p));
  }
  return ret;
}

void set_hybrid_decomposition(std::set<int> n_square_types,
                              double cutoff_regular) {
  cell_structure.set_hybrid_decomposition(comm_cart, cutoff_regular, box_geo,
                                          local_geo, n_square_types);
  on_cell_structure_change();
}

void set_regular_decomposition(bool without_ghost_force_reduction) {
  cell_structure.set_regular_decomposition(comm_cart, interaction_range(),
                                           box_geo, local_geo,
                                           without_ghost_force_reduction);
  on_cell_structure_change();
}

void cells_re_init(CellStructureType new_cs) {
  switch (new_cs) {
  case CellStructureType::CELL_STRUCTURE_REGULAR: {
    auto &current_regular_decomposition =
        dynamic_cast<RegularDecomposition const &>(
            std::as_const(cell_structure).decomposition());
    cell_structure.set_regular_decomposition(
        comm_cart, interaction_range(), box_geo, local_geo,
        current_regular_decomposition.get_without_ghost_force_reduction());
    break;
  }
  case CellStructureType::CELL_STRUCTURE_NSQUARE:
    cell_structure.set_atom_decomposition(comm_cart, box_geo, local_geo);
    break;
  case CellStructureType::CELL_STRUCTURE_HYBRID: {
    /* Get current HybridDecomposition to extract n_square_types */
    auto &current_hybrid_decomposition =
        dynamic_cast<HybridDecomposition const &>(
            std::as_const(cell_structure).decomposition());
    cell_structure.set_hybrid_decomposition(
        comm_cart, current_hybrid_decomposition.get_cutoff_regular(), box_geo,
        local_geo, current_hybrid_decomposition.get_n_square_types());
    break;
  }
  default:
    throw std::runtime_error("Unknown cell system type");
  }

  on_cell_structure_change();
}

void check_resort_particles() {
  auto const level = (cell_structure.check_resort_required(
                         cell_structure.local_particles(), skin))
                         ? Cells::RESORT_LOCAL
                         : Cells::RESORT_NONE;

  cell_structure.set_resort_particles(level);
}

void cells_update_ghosts(unsigned data_parts) {
  /* data parts that are only updated on resort */
  auto constexpr resort_only_parts =
      Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS;

  auto const global_resort =
      boost::mpi::all_reduce(comm_cart, cell_structure.get_resort_particles(),
                             std::bit_or<unsigned>());

  if (global_resort != Cells::RESORT_NONE) {
    int global = (global_resort & Cells::RESORT_GLOBAL)
                     ? CELL_GLOBAL_EXCHANGE
                     : CELL_NEIGHBOR_EXCHANGE;

    /* Resort cell system */
    cell_structure.resort_particles(global, box_geo);
    cell_structure.ghosts_count();
    cell_structure.ghosts_update(data_parts);

    /* Add the ghost particles to the index if we don't already
     * have them. */
    for (auto &p : cell_structure.ghost_particles()) {
      if (cell_structure.get_local_particle(p.id()) == nullptr) {
        cell_structure.update_particle_index(p.id(), &p);
      }
    }

    /* Particles are now sorted */
    cell_structure.clear_resort_particles();
  } else {
    /* Communication step: ghost information */
    cell_structure.ghosts_update(data_parts & ~resort_only_parts);
  }
}

Cell *find_current_cell(Particle const &p) {
  return cell_structure.find_current_cell(p);
}
