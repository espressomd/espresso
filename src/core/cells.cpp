/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *
 *  This file contains functions for the cell system.
 *
 *  Implementation of cells.hpp.
 */
#include "cells.hpp"
#include "algorithm/link_cell.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "ghosts.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "layered.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "nsquare.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "utils/NoOp.hpp"
#include "utils/mpi/gather_buffer.hpp"

#include <boost/iterator/indirect_iterator.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <boost/range/algorithm/transform.hpp>

/* Variables */

/** list of all cells. */
std::vector<Cell> cells;
/** list of pointers to all cells containing particles physically on the local
 * node. */
CellPList local_cells = {nullptr, 0, 0};
/** list of pointers to all cells containing ghosts. */
CellPList ghost_cells = {nullptr, 0, 0};

/** Type of cell structure in use */
CellStructure cell_structure = {
    CELL_STRUCTURE_NONEYET, true, {}, {}, {}, {}, nullptr, nullptr};

double max_range = 0.0;

/** On of Cells::Resort, announces the level of resort needed.
 */
unsigned resort_particles = Cells::RESORT_NONE;
int rebuild_verletlist = 1;

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

  cells_update_ghosts();

  auto pair_kernel = [&ret, &cutoff2](Particle const &p1, Particle const &p2,
                                      double dist2) {
    if (dist2 < cutoff2)
      ret.emplace_back(p1.p.identity, p2.p.identity);
  };

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    Algorithm::link_cell(boost::make_indirect_iterator(local_cells.begin()),
                         boost::make_indirect_iterator(local_cells.end()),
                         Utils::NoOp{}, pair_kernel,
                         [](Particle const &p1, Particle const &p2) {
                           return distance2(p1.r.p, p2.r.p);
                         });
    break;
  case CELL_STRUCTURE_NSQUARE:
    Algorithm::link_cell(boost::make_indirect_iterator(local_cells.begin()),
                         boost::make_indirect_iterator(local_cells.end()),
                         Utils::NoOp{}, pair_kernel,
                         [](Particle const &p1, Particle const &p2) {
                           double vec21[3];
                           get_mi_vector(vec21, p1.r.p, p2.r.p);
                           return sqrlen(vec21);
                         });
    break;
  case CELL_STRUCTURE_LAYERED:
    Algorithm::link_cell(boost::make_indirect_iterator(local_cells.begin()),
                         boost::make_indirect_iterator(local_cells.end()),
                         Utils::NoOp{}, pair_kernel,
                         [](Particle const &p1, Particle const &p2) {
                           double vec21[3];
                           get_mi_vector(vec21, p1.r.p, p2.r.p);
                           vec21[2] = p1.r.p[2] - p2.r.p[2];

                           return sqrlen(vec21);
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
void topology_init(int cs) {
  /** broadcast the flag for using Verlet list */
  boost::mpi::broadcast(comm_cart, cell_structure.use_verlet_list, 0);

  switch (cs) {
  case CELL_STRUCTURE_NONEYET:
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init();
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_init();
    break;
  default:
    fprintf(stderr,
            "INTERNAL ERROR: attempting to sort the particles in an "
            "unknown way (%d)\n",
            cs);
    errexit();
  }
}

/*@}*/

/************************************************************
 *            Exported Functions                            *
 ************************************************************/

/************************************************************/

void cells_re_init(int new_cs) {
  CellPList tmp_local;

  CELL_TRACE(fprintf(stderr, "%d: cells_re_init: convert type (%d->%d)\n",
                     this_node, cell_structure.type, new_cs));

  topology_release(cell_structure.type);
  /* MOVE old local_cell list to temporary buffer */
  memmove(&tmp_local, &local_cells, sizeof(CellPList));
  init_cellplist(&local_cells);

  /* MOVE old cells to temporary buffer */
  auto tmp_cells = std::move(cells);

  topology_init(new_cs);

  /* Add the particles to the new cs */
  cells_resort_particles(CELL_GLOBAL_EXCHANGE, tmp_local);

  /* finally deallocate the old cells */
  realloc_cellplist(&tmp_local, 0);

  for (auto &cell : tmp_cells) {
    cell.resize(0);
  }

  CELL_TRACE(fprintf(stderr, "%d: old cells deallocated\n", this_node));

  on_cell_structure_change();
}

/************************************************************/

void realloc_cells(int size) {
  CELL_TRACE(fprintf(stderr, "%d: realloc_cells %d\n", this_node, size));
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

void announce_resort_particles() {
  MPI_Allreduce(MPI_IN_PLACE, &resort_particles, 1, MPI_UNSIGNED, MPI_BOR,
                comm_cart);

  INTEG_TRACE(fprintf(stderr,
                      "%d: announce_resort_particles: resort_particles=%u\n",
                      this_node, resort_particles));
}

/*************************************************/

int cells_get_n_particles() {
  return std::accumulate(local_cells.begin(), local_cells.end(), 0,
                         [](int n, const Cell *c) { return n + c->n; });
}

/*************************************************/

namespace {
/**
 * @brief Fold coordinates to box and reset the old position.
 */
void fold_and_reset(Particle &p) {
  fold_position(p.r.p, p.l.i);

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

      auto target_cell = cs.position_to_cell(p.r.p);

      if (target_cell == nullptr) {
        append_particle(&displaced_parts, extract_particle(c, i));

        if (i < c->n) {
          i--;
        }
      } else if (target_cell != c) {
        move_particle(target_cell, c, i);

        if (i < c->n) {
          i--;
        }
      }
    }
  }

  return displaced_parts;
}

static ParticleDiff update_ghosts() {
    using Utils::make_span;

    auto part_diff = ParticleDiff{};

    for(auto &gc: ghost_cells) {
        auto parts = make_span(gc->part, gc->n);
        part_diff.removed.reserve(part_diff.removed.size() + parts.size());

        for(auto &gp: parts) {
            part_diff.removed.push_back(gp.identity());
            free_particle(&gp);
        }

        gc->n = 0;
    }

    ghost_communicator(&cell_structure.ghost_cells_comm);
    ghost_communicator(&cell_structure.exchange_ghosts_comm);

    for(auto &gc: ghost_cells) {
        auto parts = make_span(gc->part, gc->n);
        part_diff.added.reserve(part_diff.added.size() + parts.size());

        for(auto &gp: parts) {
            part_diff.added.push_back(&gp);
        }
    }

    return part_diff;
}

ParticleList topology_sort(int global, ParticleList& displaced_parts) {
  switch (cell_structure.type) {
    case CELL_STRUCTURE_LAYERED:
      return layered_exchange_and_sort_particles(global, &displaced_parts);
    case CELL_STRUCTURE_NSQUARE:
      return nsq_balance_particles();
    case CELL_STRUCTURE_DOMDEC:
      return dd_exchange_and_sort_particles(global, &displaced_parts);
  }
}

void cells_resort_particles(int global_flag, CellPList local_cells) {
  CELL_TRACE(fprintf(stderr, "%d: entering cells_resort_particles %d\n",
                     this_node, global_flag));
  auto part_diff = ParticleDiff{{}, {}};

  ParticleList displaced_parts =
      sort_and_fold_parts(cell_structure, local_cells);

  part_diff.removed.reserve(displaced_parts.n);
  for(auto const&p : Utils::make_span(displaced_parts.part, displaced_parts.n)) {
      part_diff.removed.emplace_back(p.identity());
  }

  ParticleList new_parts = topology_sort(global_flag, displaced_parts);

  std::vector<std::pair<Cell *, int>> updated_cells;
  updated_cells.reserve(new_parts.n);

  for(int p = 0; p < new_parts.n; p++) {
    auto cell = cell_structure.position_to_cell(new_parts.part[p].r.p);
    append_particle(cell, std::move(new_parts.part[p]));

    /* The new particle is now the last part in the cell */
    updated_cells.emplace_back(cell, cell->n - 1);
  }

  part_diff.added.resize(updated_cells.size());
  boost::transform(updated_cells, part_diff.added.begin(), [](auto &u) -> Particle * {
    return &(u.first->part[u.second]);
  });

  if (0 != displaced_parts.n) {
    for (int i = 0; i < displaced_parts.n; i++) {
      auto &part = displaced_parts.part[i];
      runtimeErrorMsg() << "Particle " << part.identity()
                        << " moved more than"
                           " one local box length in one timestep.";
      append_particle(local_cells.cell[0], std::move(part));
    }
    resort_particles = Cells::RESORT_GLOBAL;
  }

  /* Particles are now sorted */
  resort_particles = Cells::RESORT_NONE;

  on_resort_particles(part_diff);
}

void cells_resort_particles(int global_flag) {
  cells_resort_particles(global_flag, local_cells);
}

/*************************************************/

void cells_on_geometry_change(int flags) {
  if (max_cut > 0.0) {
    max_range = max_cut + skin;
  } else
    /* if no interactions yet, we also don't need a skin */
    max_range = 0.0;

  CELL_TRACE(fprintf(stderr, "%d: on_geometry_change with max range %f\n",
                     this_node, max_range));

  switch (cell_structure.type) {
  case CELL_STRUCTURE_DOMDEC:
    dd_on_geometry_change(flags);
    break;
  case CELL_STRUCTURE_LAYERED:
    /* there is no fast version, always redo everything. */
    cells_re_init(CELL_STRUCTURE_LAYERED);
    break;
  case CELL_STRUCTURE_NSQUARE:
    break;
  }
}

/*************************************************/

void check_resort_particles() {
  const double skin2 = Utils::sqr(skin / 2.0);

  resort_particles |= (std::any_of(local_cells.particles().begin(),
                                   local_cells.particles().end(),
                                   [&skin2](Particle const &p) {
                                     return distance2(p.r.p, p.l.p_old) > skin2;
                                   }))
                          ? Cells::RESORT_LOCAL
                          : Cells::RESORT_NONE;

  announce_resort_particles();
}

/*************************************************/
void cells_update_ghosts() {
  if (resort_particles) {
    int global = (resort_particles & Cells::RESORT_GLOBAL)
                     ? CELL_GLOBAL_EXCHANGE
                     : CELL_NEIGHBOR_EXCHANGE;

    /* Communication step: Exchange particles */
    cells_resort_particles(global);

    /* Communication step:  number of ghosts and ghost information */
    auto ghosts_change = update_ghosts();
    on_ghost_particles_change(ghosts_change);
  } else
    /* Communication step: ghost information */
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
}

Cell *find_current_cell(const Particle &p) {
  assert(not resort_particles);

  if (p.l.ghost) {
    return nullptr;
  }

  return cell_structure.position_to_cell(p.l.p_old);
}
