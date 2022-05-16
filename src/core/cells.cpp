/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "Particle.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_node.hpp"

#include "ParticleDecomposition.hpp"
#include "RegularDecomposition.hpp"

#include <utils/as_const.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <algorithm>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/serialization/set.hpp>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

/** Type of cell structure in use */
CellStructure cell_structure;

/**
 * @brief Get pairs of particles that are closer than a distance and fulfill a
 * filter criterion.
 *
 * It uses link_cell to get pairs out of the cellsystem
 * by a simple distance criterion and
 *
 * Pairs are sorted so that first.id < second.id
 */
template <class Filter>
std::vector<std::pair<int, int>> get_pairs_filtered(double const distance,
                                                    Filter filter) {
  std::vector<std::pair<int, int>> ret;
  on_observable_calc();
  auto const cutoff2 = distance * distance;
  auto pair_kernel = [&ret, &cutoff2, &filter](Particle const &p1,
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

namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive &ar, PairInfo &p, const unsigned int /* version */) {
  ar &p.id1;
  ar &p.id2;
  ar &p.pos1;
  ar &p.pos2;
  ar &p.vec21;
  ar &p.node;
}
} // namespace serialization
} // namespace boost

namespace detail {
static void search_distance_sanity_check(double const distance) {
  /* get_pairs_filtered() finds pairs via the non_bonded_loop. The maximum
   * finding range is therefore limited by the decomposition that is used.
   */
  auto const range = *boost::min_element(cell_structure.max_range());
  if (distance > range) {
    throw std::domain_error("pair search distance " + std::to_string(distance) +
                            " bigger than the decomposition range " +
                            std::to_string(range));
  }
}
} // namespace detail

std::vector<std::pair<int, int>> get_pairs(double const distance) {
  detail::search_distance_sanity_check(distance);
  auto pairs =
      get_pairs_filtered(distance, [](Particle const &) { return true; });
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

std::vector<std::pair<int, int>>
get_pairs_of_types(double const distance, std::vector<int> const &types) {
  detail::search_distance_sanity_check(distance);
  auto pairs = get_pairs_filtered(distance, [types](Particle const &p) {
    return std::any_of(types.begin(), types.end(),
                       [p](int const type) { return p.type() == type; });
  });
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

std::vector<PairInfo> non_bonded_loop_trace() {
  std::vector<PairInfo> pairs;
  auto pair_kernel = [&pairs](Particle const &p1, Particle const &p2,
                              Distance const &d) {
    pairs.emplace_back(p1.id(), p2.id(), p1.pos(), p2.pos(), d.vec21,
                       comm_cart.rank());
  };
  cell_structure.non_bonded_loop(pair_kernel);
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

void set_hybrid_decomposition(std::set<int> n_square_types,
                              double cutoff_regular) {
  cell_structure.set_hybrid_decomposition(comm_cart, cutoff_regular, box_geo,
                                          local_geo, n_square_types);
  on_cell_structure_change();
}

std::pair<std::size_t, std::size_t>
hybrid_parts_per_decomposition(HybridDecomposition const &hd) {
  Utils::Vector<std::size_t, 2> acc = {};
  boost::mpi::reduce(comm_cart, hd.parts_per_decomposition_local(), acc,
                     std::plus<>{}, 0);
  return std::make_pair(acc[0], acc[1]);
}

void cells_re_init(CellStructureType new_cs) {
  switch (new_cs) {
  case CellStructureType::CELL_STRUCTURE_REGULAR:
    cell_structure.set_regular_decomposition(comm_cart, interaction_range(),
                                             box_geo, local_geo);
    break;
  case CellStructureType::CELL_STRUCTURE_NSQUARE:
    cell_structure.set_atom_decomposition(comm_cart, box_geo, local_geo);
    break;
  case CellStructureType::CELL_STRUCTURE_HYBRID: {
    /* Get current HybridDecomposition to extract n_square_types */
    auto &current_hybrid_decomposition =
        dynamic_cast<const HybridDecomposition &>(
            Utils::as_const(cell_structure).decomposition());
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
    for (auto &part : cell_structure.ghost_particles()) {
      if (cell_structure.get_local_particle(part.id()) == nullptr) {
        cell_structure.update_particle_index(part.identity(), &part);
      }
    }

    /* Particles are now sorted */
    cell_structure.clear_resort_particles();
  } else {
    /* Communication step: ghost information */
    cell_structure.ghosts_update(data_parts & ~resort_only_parts);
  }
}

Cell *find_current_cell(const Particle &p) {
  return cell_structure.find_current_cell(p);
}
