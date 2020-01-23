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
#include "layered.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nsquare.hpp"
#include "particle_data.hpp"

#include <utils/NoOp.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/iterator/indirect_iterator.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>

/** list of all cells. */
std::vector<Cell> cells;

/** Type of cell structure in use */
CellStructure cell_structure;

/** One of @ref Cells::Resort, announces the level of resort needed.
 */
unsigned resort_particles = Cells::RESORT_NONE;
bool rebuild_verletlist = true;

CellPList CellStructure::local_cells() {
  return {m_local_cells.data(), static_cast<int>(m_local_cells.size())};
}

CellPList CellStructure::ghost_cells() {
  return {m_ghost_cells.data(), static_cast<int>(m_ghost_cells.size())};
}

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

  cells_update_ghosts(GHOSTTRANS_POSITION | GHOSTTRANS_PROPRTS);

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
  case CELL_STRUCTURE_LAYERED:
    Algorithm::link_cell(
        boost::make_indirect_iterator(cell_structure.m_local_cells.begin()),
        boost::make_indirect_iterator(cell_structure.m_local_cells.end()),
        Utils::NoOp{}, pair_kernel, [](Particle const &p1, Particle const &p2) {
          auto vec21 = get_mi_vector(p1.r.p, p2.r.p, box_geo);
          vec21[2] = p1.r.p[2] - p2.r.p[2];

          return vec21.norm2();
        });
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

/** Choose the topology release function of a certain cell system. */
static void topology_release(int cs) {
  switch (cs) {
  case CELL_STRUCTURE_NONEYET:
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_release(cell_structure.type);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_release();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_release();
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_release();
    break;
  default:
    fprintf(stderr,
            "INTERNAL ERROR: attempting to sort the particles in an "
            "unknown way (%d)\n",
            cs);
    errexit();
  }
}

/** Choose the topology init function of a certain cell system. */
void topology_init(int cs, double range, CellPList local) {
  /** broadcast the flag for using Verlet list */
  boost::mpi::broadcast(comm_cart, cell_structure.use_verlet_list, 0);

  switch (cs) {
  /* Default to DD */
  case CELL_STRUCTURE_NONEYET:
    topology_init(CELL_STRUCTURE_DOMDEC, range, local);
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type, range, local);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(&local, node_grid, range);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init(&local);
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_init(&local, node_grid, range);
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
  case CELL_STRUCTURE_LAYERED:
    return boost::mpi::all_reduce(comm_cart, local_resort,
                                  std::bit_or<unsigned>());
  default:
    return true;
  }
}

/** Go through ghost cells and remove the ghost entries from \ref
    local_particles. */
static void invalidate_ghosts() {
  for (auto const &p : cell_structure.ghost_cells().particles()) {
    if (local_particles[p.identity()] == &p) {
      local_particles[p.identity()] = {};
    }
  }

  for (auto &c : cell_structure.m_ghost_cells) {
    c->n = 0;
  }
}

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

/************************************************************/

void cells_re_init(int new_cs, double range) {
  invalidate_ghosts();

  topology_release(cell_structure.type);
  /* MOVE old local_cell list to temporary buffer */
  std::vector<Cell *> old_local_cells;
  std::swap(old_local_cells, cell_structure.m_local_cells);

  /* MOVE old cells to temporary buffer */
  auto tmp_cells = std::move(cells);

  topology_init(
      new_cs, range,
      {old_local_cells.data(), static_cast<int>(old_local_cells.size())});
  cell_structure.min_range = range;

  clear_particle_node();

  for (auto &cell : tmp_cells) {
    cell.resize(0);
  }

  /* to enforce initialization of the ghost cells */
  resort_particles = Cells::RESORT_GLOBAL;

  on_cell_structure_change();
}

/************************************************************/

void realloc_cells(int size) {
  /* free all memory associated with cells to be deleted. */
  for (auto &c : cells) {
    c.resize(0);
  }
  /* resize the cell list */
  cells.resize(size);
}

/*************************************************/

void set_resort_particles(Cells::Resort level) {
  resort_particles |= level;
  assert(resort_particles & level);
}

unsigned const &get_resort_particles() { return resort_particles; }

/*************************************************/

int cells_get_n_particles() {
  return std::accumulate(cell_structure.m_local_cells.begin(),
                         cell_structure.m_local_cells.end(), 0,
                         [](int n, const Cell *c) { return n + c->n; });
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
} // namespace

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
ParticleList sort_and_fold_parts(const CellStructure &cs, CellPList cells) {
  ParticleList displaced_parts;

  for (auto &c : cells) {
    for (int i = 0; i < c->n; i++) {
      auto &p = c->part[i];

      fold_and_reset(p);

      auto target_cell = cs.particle_to_cell(p);

      if (target_cell == nullptr) {
        append_unindexed_particle(&displaced_parts,
                                  extract_indexed_particle(c, i));

        if (i < c->n) {
          i--;
        }
      } else if (target_cell != c) {
        move_indexed_particle(target_cell, c, i);

        if (i < c->n) {
          i--;
        }
      }
    }
  }

  return displaced_parts;
}

void cells_resort_particles(int global_flag) {

  invalidate_ghosts();

  clear_particle_node();
  n_verlet_updates++;

  ParticleList displaced_parts =
      sort_and_fold_parts(cell_structure, cell_structure.local_cells());

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED: {
    layered_exchange_and_sort_particles(global_flag, &displaced_parts);
  } break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_exchange_particles(global_flag, &displaced_parts);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_exchange_and_sort_particles(global_flag, &displaced_parts, node_grid);
    break;
  }

  if (0 != displaced_parts.n) {
    for (int i = 0; i < displaced_parts.n; i++) {
      auto &part = displaced_parts.part[i];
      runtimeErrorMsg() << "Particle " << part.identity()
                        << " moved more than"
                           " one local box length in one timestep.";
      resort_particles = Cells::RESORT_GLOBAL;
      append_indexed_particle(cell_structure.m_local_cells[0], std::move(part));
    }
  } else {
#ifdef ADDITIONAL_CHECKS
    /* at the end of the day, everything should be consistent again */
    check_particle_consistency();
    check_particle_sorting();
#endif
  }

  rebuild_verletlist = true;

  displaced_parts.clear();

  on_resort_particles(cell_structure.local_cells().particles());
}

/*************************************************/

void cells_on_geometry_change(int flags) {
  /* Consider skin only if there are actually interactions */
  auto const range = (max_cut > 0.) ? max_cut + skin : INACTIVE_CUTOFF;
  cell_structure.min_range = range;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    dd_on_geometry_change(flags, node_grid, range);
    break;
  case CELL_STRUCTURE_LAYERED:
    /* there is no fast version, always redo everything. */
    cells_re_init(CELL_STRUCTURE_LAYERED, range);
    break;
  case CELL_STRUCTURE_NSQUARE:
    break;
  }
}

/*************************************************/

void check_resort_particles() {
  const double skin2 = Utils::sqr(skin / 2.0);

  resort_particles |=
      (std::any_of(cell_structure.local_cells().particles().begin(),
                   cell_structure.local_cells().particles().end(),
                   [&skin2](Particle const &p) {
                     return (p.r.p - p.l.p_old).norm2() > skin2;
                   }))
          ? Cells::RESORT_LOCAL
          : Cells::RESORT_NONE;
}

/*************************************************/
void cells_update_ghosts(unsigned data_parts) {
  /* data parts that are only updated on resort */
  auto constexpr resort_only_parts = GHOSTTRANS_PROPRTS | GHOSTTRANS_BONDS;

  auto const global_resort =
      topology_check_resort(cell_structure.type, resort_particles);

  if (global_resort != Cells::RESORT_NONE) {
    int global = (global_resort & Cells::RESORT_GLOBAL)
                     ? CELL_GLOBAL_EXCHANGE
                     : CELL_NEIGHBOR_EXCHANGE;

    /* Resort cell system */
    cells_resort_particles(global);

    /* Communication step: number of ghosts and ghost information */
    ghost_communicator(&cell_structure.exchange_ghosts_comm,
                       GHOSTTRANS_PARTNUM);
    ghost_communicator(&cell_structure.exchange_ghosts_comm, data_parts);

    /* Particles are now sorted */
    resort_particles = Cells::RESORT_NONE;
  } else {
    /* Communication step: ghost information */
    ghost_communicator(&cell_structure.exchange_ghosts_comm,
                       data_parts & ~resort_only_parts);
  }
}

Cell *find_current_cell(const Particle &p) {
  assert(not resort_particles);

  if (p.l.ghost) {
    return nullptr;
  }

  return cell_structure.particle_to_cell(p);
}
