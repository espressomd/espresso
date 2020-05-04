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
 *  Implements the malloc replacements as described in debug.hpp.
 */

#include "debug.hpp"

#include "cells.hpp"

#include <stdexcept>

void check_particle_consistency() {
  auto local_particles = cell_structure.local_particles();
  auto const cell_part_cnt = local_particles.size();

  auto const max_id = cell_structure.get_max_local_particle_id();

  for (auto const &p : local_particles) {
    auto const id = p.identity();

    if (id < 0 || id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    if (cell_structure.get_local_particle(id) != &p) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id */
  int local_part_cnt = 0;
  for (int n = 0; n < cell_structure.get_max_local_particle_id() + 1; n++) {
    if (cell_structure.get_local_particle(n) != nullptr) {
      local_part_cnt++;
      if (cell_structure.get_local_particle(n)->p.identity != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
      }
    }
  }

  if (local_part_cnt != cell_part_cnt) {
    throw std::runtime_error(
        std::to_string(cell_part_cnt) + " parts in cells but " +
        std::to_string(local_part_cnt) + " parts in local_particles");
  }
}

void check_particle_sorting() {
  for (auto cell : cell_structure.local_cells()) {
    for (auto const &p : cell->particles()) {
      if (cell_structure.particle_to_cell(p) != cell) {
        throw std::runtime_error("misplaced particle with id " +
                                 std::to_string(p.identity()));
      }
    }
  }
}
