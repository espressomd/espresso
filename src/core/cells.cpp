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
#include "debug.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_loop.hpp"

#include "AtomDecomposition.hpp"
#include "DomainDecomposition.hpp"

#include <utils/NoOp.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <cstdio>

/** Type of cell structure in use */
CellStructure cell_structure;

/**
 * @brief Get pairs closer than distance from the cells.
 *
 * This is mostly for testing purposes and uses link_cell
 * to get pairs out of the cellsystem by a simple distance
 * criterion.
 *
 * Pairs are sorted so that first.id < second.id
 */
std::vector<std::pair<int, int>> get_pairs(double distance) {
  std::vector<std::pair<int, int>> ret;
  auto const cutoff2 = distance * distance;

  cells_update_ghosts(Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES);

  auto pair_kernel = [&ret, &cutoff2](Particle const &p1, Particle const &p2,
                                      Distance const &d) {
    if (d.dist2 < cutoff2)
      ret.emplace_back(p1.p.identity, p2.p.identity);
  };

  short_range_loop(Utils::NoOp{}, pair_kernel);

  /* Sort pairs */
  for (auto &pair : ret) {
    if (pair.first > pair.second)
      std::swap(pair.first, pair.second);
  }

  return ret;
}

void mpi_get_pairs_slave(int, int) {
  on_observable_calc();

  double distance;
  boost::mpi::broadcast(comm_cart, distance, 0);

  auto local_pairs = get_pairs(distance);

  Utils::Mpi::gather_buffer(local_pairs, comm_cart);
}

/**
 * @brief Collect pairs from all nodes.
 */
std::vector<std::pair<int, int>> mpi_get_pairs(double distance) {
  mpi_call(mpi_get_pairs_slave, 0, 0);
  on_observable_calc();

  boost::mpi::broadcast(comm_cart, distance, 0);

  auto pairs = get_pairs(distance);

  Utils::Mpi::gather_buffer(pairs, comm_cart);

  return pairs;
}

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

/************************************************************/

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

/*************************************************/

void cells_resort_particles(int global_flag) {
  n_verlet_updates++;

  cell_structure.resort_particles(global_flag);

#ifdef ADDITIONAL_CHECKS
  /* at the end of the day, everything should be consistent again */
  check_particle_consistency();
  check_particle_sorting();
#endif
}

/*************************************************/

void check_resort_particles() {
  const double skin2 = Utils::sqr(skin / 2.0);

  auto const level = (std::any_of(cell_structure.local_particles().begin(),
                                  cell_structure.local_particles().end(),
                                  [&skin2](Particle const &p) {
                                    return (p.r.p - p.l.p_old).norm2() > skin2;
                                  }))
                         ? Cells::RESORT_LOCAL
                         : Cells::RESORT_NONE;

  cell_structure.set_resort_particles(level);
}

/*************************************************/
void cells_update_ghosts(unsigned data_parts) {
  /* data parts that are only updated on resort */
  auto constexpr resort_only_parts =
      Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS;

  unsigned int localResort = cell_structure.get_resort_particles();
  auto const global_resort =
      boost::mpi::all_reduce(comm_cart, localResort, std::bit_or<unsigned>());

  if (global_resort != Cells::RESORT_NONE) {
    int global = (global_resort & Cells::RESORT_GLOBAL)
                     ? CELL_GLOBAL_EXCHANGE
                     : CELL_NEIGHBOR_EXCHANGE;

    /* Resort cell system */
    cells_resort_particles(global);

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
  assert(not cell_structure.get_resort_particles());

  if (p.l.ghost) {
    return nullptr;
  }

  return cell_structure.particle_to_cell(p);
}

const DomainDecomposition *get_domain_decomposition() {
  return &dynamic_cast<const DomainDecomposition &>(
      Utils::as_const(cell_structure).decomposition());
}

void cells_set_use_verlet_lists(bool use_verlet_lists) {
  cell_structure.use_verlet_list = use_verlet_lists;
}