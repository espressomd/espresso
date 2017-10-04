/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file cells.cpp
 *
 *  This file contains functions for the cell system.
 *
 *  For more information on cells, see cells.hpp
 *   */
#include "cells.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "ghosts.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "layered.hpp"
#include "lees_edwards_domain_decomposition.hpp"
#include "nsquare.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "verlet.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>

/* Variables */

/** list of all cells. */
std::vector<Cell> cells;
/** list of pointers to all cells containing particles physically on the local
 * node. */
CellPList local_cells = {NULL, 0, 0};
/** list of pointers to all cells containing ghosts. */
CellPList ghost_cells = {NULL, 0, 0};

/** Type of cell structure in use */
CellStructure cell_structure = {CELL_STRUCTURE_NONEYET};

double max_range = 0.0;

int rebuild_verletlist = 0;

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Switch for choosing the topology release function of a certain
    cell system. */
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
    fprintf(stderr, "INTERNAL ERROR: attempting to sort the particles in an "
                    "unknown way (%d)\n",
            cs);
    errexit();
  }
}

/** Switch for choosing the topology init function of a certain
    cell system. */
static void topology_init(int cs, CellPList *local) {
  switch (cs) {
  case CELL_STRUCTURE_NONEYET:
    break;
  case CELL_STRUCTURE_CURRENT:
    topology_init(cell_structure.type, local);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_topology_init(local);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_topology_init(local);
    break;
  case CELL_STRUCTURE_LAYERED:
    layered_topology_init(local);
    break;
  default:
    fprintf(stderr, "INTERNAL ERROR: attempting to sort the particles in an "
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

  invalidate_ghosts();

  topology_release(cell_structure.type);
  /* MOVE old local_cell list to temporary buffer */
  memmove(&tmp_local, &local_cells, sizeof(CellPList));
  init_cellplist(&local_cells);

  /* MOVE old cells to temporary buffer */
  auto tmp_cells = std::move(cells);

  topology_init(new_cs, &tmp_local);

  particle_invalidate_part_node();

  /* finally deallocate the old cells */
  realloc_cellplist(&tmp_local, 0);

  for (auto &cell : tmp_cells) {
    realloc_particlelist(&cell, 0);
  }

  CELL_TRACE(fprintf(stderr, "%d: old cells deallocated\n", this_node));

  /* to enforce initialization of the ghost cells */
  resort_particles = 1;

  on_cell_structure_change();
}

/************************************************************/

void realloc_cells(int size) {
  CELL_TRACE(fprintf(stderr, "%d: realloc_cells %d\n", this_node, size));
  /* free all memory associated with cells to be deleted. */
  for (int i = size; i < cells.size(); i++) {
    realloc_particlelist(&cells[i], 0);
  }

  auto const old_size = cells.size();

  /* resize the cell list */
  cells.resize(size);

  /* initialize new cells */
  for (int i = old_size; i < size; i++) {
    init_particlelist(&cells[i]);
  }
}

/*************************************************/

void announce_resort_particles() {
  int sum;

  MPI_Allreduce(&resort_particles, &sum, 1, MPI_INT, MPI_SUM, comm_cart);
  resort_particles = (sum > 0) ? 1 : 0;

  INTEG_TRACE(fprintf(stderr,
                      "%d: announce_resort_particles: resort_particles=%d\n",
                      this_node, resort_particles));
}

/*************************************************/

int cells_get_n_particles() {
  using std::distance;
  return distance(local_cells.particles().begin(),
                  local_cells.particles().end());
}

/*************************************************/

void cells_resort_particles(int global_flag) {
  CELL_TRACE(fprintf(stderr, "%d: entering cells_resort_particles %d\n",
                     this_node, global_flag));

  invalidate_ghosts();

  particle_invalidate_part_node();
  n_verlet_updates++;

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_exchange_and_sort_particles(global_flag);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_balance_particles(global_flag);
    break;
  case CELL_STRUCTURE_DOMDEC:
    dd_exchange_and_sort_particles(global_flag);
    break;
  }

#ifdef ADDITIONAL_CHECKS
  /* at the end of the day, everything should be consistent again */
  check_particle_consistency();
#endif

  ghost_communicator(&cell_structure.ghost_cells_comm);
  ghost_communicator(&cell_structure.exchange_ghosts_comm);

  resort_particles = 0;
  rebuild_verletlist = 1;

  on_resort_particles();

  CELL_TRACE(
      fprintf(stderr, "%d: leaving cells_resort_particles\n", this_node));
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
    /* this cell system doesn't need to react, just tell
       the others */
    on_boxl_change();
    break;
  }
}

/*************************************************/

void check_resort_particles() {
  const double skin2 = SQR(skin / 2.0);

  resort_particles =
      std::any_of(local_cells.particles().begin(),
                  local_cells.particles().end(), [&skin2](Particle const &p) {
                    return distance2(p.r.p, p.l.p_old) > skin2;
                  });

  announce_resort_particles();
}

/*************************************************/
void cells_update_ghosts() {
  /* if dd.use_vList is set, it so far means we want EXACT sorting of the
   * particles.*/
  if (dd.use_vList == 0)
    resort_particles = 1;

  if (resort_particles) {
#ifdef LEES_EDWARDS
    /* Communication step:  number of ghosts and ghost information */
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);
#else
    /* Communication step:  number of ghosts and ghost information */
    cells_resort_particles(CELL_NEIGHBOR_EXCHANGE);
#endif
  } else
    /* Communication step: ghost information */
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
}
