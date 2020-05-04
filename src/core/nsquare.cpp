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
 *  Implementation of nsquare.hpp.
 */

#include "nsquare.hpp"
#include "cells.hpp"

#include "AtomDecomposition.hpp"

#include <boost/optional.hpp>

#include <cassert>
boost::optional<AtomDecomposition> ad;

void nsq_topology_init(const boost::mpi::communicator &comm) {
  ad = AtomDecomposition(comm);

  cell_structure.m_local_cells = ad->local_cells();
  cell_structure.m_ghost_cells = ad->ghost_cells();

  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.particle_to_cell = [](const Particle &p) {
    return assert(ad), ad->particle_to_cell(p);
  };

  cell_structure.max_range = ad->max_range();
  cell_structure.exchange_ghosts_comm = ad->exchange_ghosts_comm();
  cell_structure.collect_ghost_force_comm = ad->collect_ghost_force_comm();
}

void nsq_exchange_particles(int global_flag, ParticleList *displaced_parts,
                            std::vector<Cell *> &modified_cells) {
  assert(displaced_parts);
  ad.value().resort(global_flag, *displaced_parts, modified_cells);
}
