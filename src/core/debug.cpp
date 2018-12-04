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
    Implements the malloc replacements as described in \ref debug.hpp
   "debug.hpp".
*/

#include "debug.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <cassert>

#ifdef ONEPART_DEBUG
int check_id = ONEPART_DEBUG_ID;
#endif

namespace {
  void check_node(int this_node, const Particle & p) {
    assert((cell_structure.position_to_node(p.r.p) == this_node) && "Particle on wrong node.");
  }

  void check_cell(const Cell * cell, const Particle & p) {
    assert((cell_structure.position_to_cell(p.r.p) == cell) && "Particle in wrong cell.");
  }

  /* For local particles we check that they are in the
   * particles index, and the index entry points to the
   * correct particle, and that the local particle is not
   * flagged as ghost. */
  void check_index_local_part(const Particle &p) {
    auto p_index = local_particles[p.identity()];
    assert((p_index) && "Local particle is missing from index");
    assert((&p == p_index) && "Wrong particle index entry.");
    assert((not p.l.ghost) && "Local particle is flagged as ghost.");
  }

  /* For ghosts we can only check that they are listed in
   * the index. But since they could also be on this node
   * as local parts or from other ghost cells, we can
   * not check if the entry is valid. */
  void check_index_ghost_part(const Particle &p) {
    auto p_index = local_particles[p.identity()];
    assert((p_index) && "Ghost particle is missing from index");
    assert((p_index->identity() == p.identity()) && "Wrong particle index entry.");
  }
}

void check_particle_consistency() {
  for(auto c : local_cells) {
    for(int i = 0; i < c->n; i++) {
      auto const&p = c->part[i];
      check_node(this_node, p);
      check_cell(c, p);
      check_index_local_part(p);
    }
  }

  for(auto c : local_cells) {
    for(int i = 0; i < c->n; i++) {
      auto const&p = c->part[i];
      check_index_ghost_part(p);
      }
  }
}

