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
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"

#include "DomainDecomposition.hpp"
#include "ParticleDecomposition.hpp"

#include <utils/as_const.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <algorithm>
#include <boost/range/algorithm/min_element.hpp>
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
      ret.emplace_back(p1.p.identity, p2.p.identity);
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

std::vector<PairInfo> non_bonded_loop_trace() {
  std::vector<PairInfo> ret;
  auto pair_kernel = [&ret](Particle const &p1, Particle const &p2,
                            Distance const &d) {
    ret.emplace_back(p1.p.identity, p2.p.identity, p1.r.p, p2.r.p, d.vec21,
                     comm_cart.rank());
  };

  cell_structure.non_bonded_loop(pair_kernel);

  return ret;
}

static auto mpi_get_pairs_local(double const distance) {
  auto pairs =
      get_pairs_filtered(distance, [](Particle const &) { return true; });
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

REGISTER_CALLBACK_MAIN_RANK(mpi_get_pairs_local)

static auto mpi_get_pairs_of_types_local(double const distance,
                                         std::vector<int> const &types) {
  auto pairs = get_pairs_filtered(distance, [types](Particle const &p) {
    return std::any_of(types.begin(), types.end(),
                       [p](int const type) { return p.p.type == type; });
  });
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

REGISTER_CALLBACK_MAIN_RANK(mpi_get_pairs_of_types_local)

namespace detail {
void search_distance_sanity_check(double const distance) {
  /* get_pairs_filtered() finds pairs via the non_bonded_loop. The maximum
   * finding range is therefore limited by the decomposition that is used.
   */
  auto range = *boost::min_element(cell_structure.max_range());
  if (distance > range) {
    runtimeErrorMsg() << "pair search distance " << distance
                      << " bigger than the decomposition range " << range;
  }
}
} // namespace detail

std::vector<std::pair<int, int>> mpi_get_pairs(double const distance) {
  detail::search_distance_sanity_check(distance);
  return mpi_call(::Communication::Result::main_rank, mpi_get_pairs_local,
                  distance);
}

std::vector<std::pair<int, int>>
mpi_get_pairs_of_types(double const distance, std::vector<int> const &types) {
  detail::search_distance_sanity_check(distance);
  return mpi_call(::Communication::Result::main_rank,
                  mpi_get_pairs_of_types_local, distance, types);
}

static void non_bonded_loop_trace_local() {
  auto pairs = non_bonded_loop_trace();
  Utils::Mpi::gather_buffer(pairs, comm_cart);
}

REGISTER_CALLBACK(non_bonded_loop_trace_local)

std::vector<PairInfo> mpi_non_bonded_loop_trace() {
  mpi_call(non_bonded_loop_trace_local);
  auto pairs = non_bonded_loop_trace();
  Utils::Mpi::gather_buffer(pairs, comm_cart);
  return pairs;
}

static void mpi_resort_particles_local(int global_flag) {
  cell_structure.resort_particles(global_flag);

  boost::mpi::gather(
      comm_cart, static_cast<int>(cell_structure.local_particles().size()), 0);
}

REGISTER_CALLBACK(mpi_resort_particles_local)

std::vector<int> mpi_resort_particles(int global_flag) {
  mpi_call(mpi_resort_particles_local, global_flag);
  cell_structure.resort_particles(global_flag);

  clear_particle_node();

  std::vector<int> n_parts;
  boost::mpi::gather(comm_cart,
                     static_cast<int>(cell_structure.local_particles().size()),
                     n_parts, 0);

  return n_parts;
}

void cells_re_init(int new_cs) {
  switch (new_cs) {
  case CELL_STRUCTURE_DOMDEC:
    cell_structure.set_domain_decomposition(comm_cart, interaction_range(),
                                            box_geo, local_geo);
    break;
  case CELL_STRUCTURE_NSQUARE:
    cell_structure.set_atom_decomposition(comm_cart, box_geo);
    break;
  default:
    throw std::runtime_error("Unknown cell system type");
  }

  on_cell_structure_change();
}

REGISTER_CALLBACK(cells_re_init)

void mpi_bcast_cell_structure(int cs) { mpi_call_all(cells_re_init, cs); }

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
    cell_structure.resort_particles(global);
    cell_structure.ghosts_count();
    cell_structure.ghosts_update(data_parts);

    /* Add the ghost particles to the index if we don't already
     * have them. */
    for (auto &part : cell_structure.ghost_particles()) {
      if (cell_structure.get_local_particle(part.p.identity) == nullptr) {
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

const DomainDecomposition *get_domain_decomposition() {
  return &dynamic_cast<const DomainDecomposition &>(
      Utils::as_const(cell_structure).decomposition());
}

void mpi_set_use_verlet_lists_local(bool use_verlet_lists) {
  cell_structure.use_verlet_list = use_verlet_lists;
}

REGISTER_CALLBACK(mpi_set_use_verlet_lists_local)

void mpi_set_use_verlet_lists(bool use_verlet_lists) {
  mpi_call_all(mpi_set_use_verlet_lists_local, use_verlet_lists);
}
