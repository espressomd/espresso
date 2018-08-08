/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file topology.cpp
 *
 *  This file contains functions for handling the system topology.
 *
 *  For more information see topology.hpp
 *   */

#include "topology.hpp"
#include "particle_data.hpp"

std::vector<Molecule> topology;
int topo_part_info_synced = 0;

void realloc_topology(int size) {
  topology.resize(size);
  topo_part_info_synced = 0;
}

// Parallel function for synchronising topology and particle data
void sync_topo_part_info() {
  for (unsigned molid = 0; molid < topology.size(); molid++) {
    auto const &mol = topology.at(molid);
    for (auto const &pid : mol.part) {
      auto p = local_particles[pid];

      if (p) {
        p->p.mol_id = molid;
      }
    }
  }

  topo_part_info_synced = 1;
}
