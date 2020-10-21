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
#include "particle_data.hpp"

#include "DomainDecomposition.hpp"
#include "ParticleDecomposition.hpp"

#include <utils/math/sqr.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/sort.hpp>

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

  cell_structure.non_bonded_loop(pair_kernel);

  /* Sort pairs */
  for (auto &pair : ret) {
    if (pair.first > pair.second)
      std::swap(pair.first, pair.second);
  }

  return ret;
}

void mpi_get_pairs_local(int, int) {
  on_observable_calc();

  double distance;
  boost::mpi::broadcast(comm_cart, distance, 0);

  auto local_pairs = get_pairs(distance);

  Utils::Mpi::gather_buffer(local_pairs, comm_cart);
}

REGISTER_CALLBACK(mpi_get_pairs_local)

std::vector<std::pair<int, int>> mpi_get_pairs(double distance) {
  mpi_call(mpi_get_pairs_local, 0, 0);
  on_observable_calc();

  boost::mpi::broadcast(comm_cart, distance, 0);

  auto pairs = get_pairs(distance);

  Utils::Mpi::gather_buffer(pairs, comm_cart);

  return pairs;
}

void mpi_resort_particles_local(int global_flag, int) {
  cell_structure.resort_particles(global_flag);

  boost::mpi::gather(
      comm_cart, static_cast<int>(cell_structure.local_particles().size()), 0);
}

REGISTER_CALLBACK(mpi_resort_particles_local)

std::vector<int> mpi_resort_particles(int global_flag) {
  mpi_call(mpi_resort_particles_local, global_flag, 0);
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
