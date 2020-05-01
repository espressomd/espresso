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
#include "algorithm/link_cell.hpp"
#include "communication.hpp"
#include "debug.hpp"
#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "ghosts.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nsquare.hpp"

#include <utils/NoOp.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/sort.hpp>

#include <cstdio>

/** list of all cells. */
std::vector<Cell> cells;

/** Type of cell structure in use */
CellStructure cell_structure;

bool rebuild_verletlist = true;

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
                                      double dist2) {
    if (dist2 < cutoff2)
      ret.emplace_back(p1.p.identity, p2.p.identity);
  };

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    Algorithm::link_cell(
        boost::make_indirect_iterator(cell_structure.m_local_cells.begin()),
        boost::make_indirect_iterator(cell_structure.m_local_cells.end()),
        Utils::NoOp{}, pair_kernel, [](Particle const &p1, Particle const &p2) {
          return (p1.r.p - p2.r.p).norm2();
        });
    break;
  case CELL_STRUCTURE_NSQUARE:
    Algorithm::link_cell(
        boost::make_indirect_iterator(cell_structure.m_local_cells.begin()),
        boost::make_indirect_iterator(cell_structure.m_local_cells.end()),
        Utils::NoOp{}, pair_kernel, [](Particle const &p1, Particle const &p2) {
          return get_mi_vector(p1.r.p, p2.r.p, box_geo).norm2();
        });
    break;
  }

  /* Sort pairs */
  for (auto &pair : ret) {
    if (pair.first > pair.second)
      std::swap(pair.first, pair.second);
  }

  return ret;
}

void mpi_get_pairs_slave(int, int) {
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
  boost::mpi::broadcast(comm_cart, distance, 0);

  auto pairs = get_pairs(distance);

  Utils::Mpi::gather_buffer(pairs, comm_cart);

  return pairs;
}

/************************************************************/
/** \name Private Functions */
/************************************************************/
/*@{*/

/** Choose the topology init function of a certain cell system. */
void topology_init(int cs, double range) {
  /** broadcast the flag for using Verlet list */
  boost::mpi::broadcast(comm_cart, cell_structure.use_verlet_list, 0);

  switch (cs) {
  /* Default to DD */
  case CELL_STRUCTURE_NONEYET:
    topology_init(CELL_STRUCTURE_DOMDEC, range);
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type, range);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(comm_cart, range, box_geo, local_geo);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init(comm_cart);
    break;
  default:
    fprintf(stderr,
            "INTERNAL ERROR: attempting to sort the particles in an "
            "unknown way (%d)\n",
            cs);
    errexit();
  }
}

unsigned topology_check_resort(int cs, unsigned local_resort) {
  switch (cs) {
  case CELL_STRUCTURE_DOMDEC:
  case CELL_STRUCTURE_NSQUARE:
    return boost::mpi::all_reduce(comm_cart, local_resort,
                                  std::bit_or<unsigned>());
  default:
    return Cells::Resort::RESORT_GLOBAL;
  }
}

/** Go through ghost cells and remove the ghost entries from the
    local particle index. */
static void invalidate_ghosts() {
  for (auto const &p : cell_structure.ghost_particles()) {
    if (cell_structure.get_local_particle(p.identity()) == &p) {
      cell_structure.update_particle_index(p.identity(), nullptr);
    }
  }
}

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

/************************************************************/

void cells_re_init(int new_cs, double range) {
  invalidate_ghosts();

  auto local_parts = cell_structure.local_particles();
  std::vector<Particle> particles(local_parts.begin(), local_parts.end());

  cell_structure.m_local_cells.clear();
  cell_structure.m_ghost_cells.clear();
  cells.clear();

  topology_init(new_cs, range);
  cell_structure.min_range = range;

  for (auto &p : particles) {
    cell_structure.add_particle(std::move(p));
  }

  on_cell_structure_change();
}

/*************************************************/

int cells_get_n_particles() {
  return std::accumulate(
      cell_structure.m_local_cells.begin(), cell_structure.m_local_cells.end(),
      0, [](int n, const Cell *c) { return n + c->particles().size(); });
}

/*************************************************/

namespace {
/**
 * @brief Fold coordinates to box and reset the old position.
 */
void fold_and_reset(Particle &p) {
  fold_position(p.r.p, p.l.i, box_geo);

  p.l.p_old = p.r.p;
}

/**
 * @brief Sort and fold particles.
 *
 * This function folds the positions of all particles back into the
 * box and puts them back into the correct cells. Particles that do
 * not belong to this node are removed from the cell and returned.
 *
 * @param cs The cell system to be used.
 * @param cells Cells to iterate over.
 *
 * @returns List of Particles that do not belong on this node.
 */
ParticleList sort_and_fold_parts(const CellStructure &cs,
                                 Utils::Span<Cell *> cells,
                                 std::vector<Cell *> &modified_cells) {
  ParticleList displaced_parts;

  for (auto &c : cells) {
    for (auto it = c->particles().begin(); it != c->particles().end();) {
      fold_and_reset(*it);

      auto target_cell = cs.particle_to_cell(*it);

      /* Particle is in place */
      if (target_cell == c) {
        std::advance(it, 1);
        continue;
      }

      auto p = std::move(*it);
      it = c->particles().erase(it);
      modified_cells.push_back(c);

      /* Particle is not local */
      if (target_cell == nullptr) {
        displaced_parts.insert(std::move(p));
      }
      /* Particle belongs on this node but is in the wrong cell. */
      else if (target_cell != c) {
        target_cell->particles().insert(std::move(p));
        modified_cells.push_back(target_cell);
      }
    }
  }

  return displaced_parts;
}
} // namespace

void cells_resort_particles(int global_flag) {
  invalidate_ghosts();

  n_verlet_updates++;

  static std::vector<Cell *> modified_cells;
  modified_cells.clear();

  ParticleList displaced_parts = sort_and_fold_parts(
      cell_structure, cell_structure.local_cells(), modified_cells);

  for (auto const &p : displaced_parts) {
    cell_structure.update_particle_index(p.identity(), nullptr);
  }

  switch (cell_structure.type) {
  case CELL_STRUCTURE_NSQUARE:
    nsq_exchange_particles(global_flag, &displaced_parts, modified_cells);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_exchange_and_sort_particles(global_flag, &displaced_parts,
                                   modified_cells);
    break;
  }

  boost::sort(modified_cells);
  for (auto cell : modified_cells | boost::adaptors::uniqued) {
    cell_structure.update_particle_index(cell->particles());
  }

  if (not displaced_parts.empty()) {
    auto sort_cell = cell_structure.m_local_cells[0];

    for (auto &part : displaced_parts) {
      runtimeErrorMsg() << "Particle " << part.identity()
                        << " moved more than"
                           " one local box length in one timestep.";
      sort_cell->particles().insert(std::move(part));
    }

    cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
    cell_structure.update_particle_index(sort_cell->particles());
  } else {
#ifdef ADDITIONAL_CHECKS
    /* at the end of the day, everything should be consistent again */
    check_particle_consistency();
    check_particle_sorting();
#endif
  }

  rebuild_verletlist = true;

  displaced_parts.clear();

  on_resort_particles();
}

/*************************************************/

void cells_on_geometry_change(bool fast) {
  auto const range = interaction_range();
  cell_structure.min_range = range;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    dd_on_geometry_change(fast, range, box_geo, local_geo);
    break;
  case CELL_STRUCTURE_NSQUARE:
    break;
  }
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

  auto const global_resort = topology_check_resort(
      cell_structure.type, cell_structure.get_resort_particles());

  if (global_resort != Cells::RESORT_NONE) {
    int global = (global_resort & Cells::RESORT_GLOBAL)
                     ? CELL_GLOBAL_EXCHANGE
                     : CELL_NEIGHBOR_EXCHANGE;

    /* Resort cell system */
    cells_resort_particles(global);

    /* Communication step: number of ghosts and ghost information */
    ghost_communicator(&cell_structure.exchange_ghosts_comm,
                       GHOSTTRANS_PARTNUM);
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
