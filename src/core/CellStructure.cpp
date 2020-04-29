/*
 * Copyright (C) 2010-2020 The ESPResSo project
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

#include "CellStructure.hpp"

#include <utils/contains.hpp>

void CellStructure::remove_particle(int id) {
  auto remove_all_bonds_to = [id](BondList &bl) {
    for (auto it = bl.begin(); it != bl.end();) {
      if (Utils::contains(it->partner_ids(), id)) {
        it = bl.erase(it);
      } else {
        std::advance(it, 1);
      }
    }
  };

  for (auto c : m_local_cells) {
    auto &parts = c->particles();

    for (auto it = parts.begin(); it != parts.end();) {
      if (it->identity() == id) {
        it = parts.erase(it);
        update_particle_index(id, nullptr);
        update_particle_index(parts);
      } else {
        remove_all_bonds_to(it->bonds());
        it++;
      }
    }
  }
}

Particle *CellStructure::add_local_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  if (sort_cell) {

    return std::addressof(
        append_indexed_particle(sort_cell->particles(), std::move(p)));
  }

  return {};
}

Particle *CellStructure::add_particle(Particle &&p) {
  auto const sort_cell = particle_to_cell(p);
  /* There is always at least one cell, so if the particle
   * does not belong to a cell on this node we can put it there. */
  auto cell = sort_cell ? sort_cell : local_cells()[0];

  /* If the particle isn't local a global resort may be
   * needed, otherwise a local resort if sufficient. */
  set_resort_particles(sort_cell ? Cells::RESORT_LOCAL : Cells::RESORT_GLOBAL);

  return std::addressof(
      append_indexed_particle(cell->particles(), std::move(p)));
}

int CellStructure::get_max_local_particle_id() const {
  auto it = std::find_if(m_particle_index.rbegin(), m_particle_index.rend(),
                         [](const Particle *p) { return p != nullptr; });

  return (it != m_particle_index.rend()) ? (*it)->identity() : -1;
}

void CellStructure::remove_all_particles() {
  for (auto c : m_local_cells) {
    c->particles().clear();
  }

  m_particle_index.clear();
}

/* Map the data parts flags from cells to those used internally
 * by the ghost communication */
static unsigned map_data_parts(unsigned data_parts) {
  using namespace Cells;

  /* clang-format off */
  return GHOSTTRANS_NONE
         | ((DATA_PART_PROPERTIES & data_parts) ? GHOSTTRANS_PROPRTS : 0u)
         | ((DATA_PART_POSITION & data_parts) ? GHOSTTRANS_POSITION : 0u)
         | ((DATA_PART_MOMENTUM & data_parts) ? GHOSTTRANS_MOMENTUM : 0u)
         | ((DATA_PART_FORCE & data_parts) ? GHOSTTRANS_FORCE : 0u)
         | ((DATA_PART_BONDS & data_parts) ? GHOSTTRANS_BONDS : 0u);
  /* clang-format on */
}

void CellStructure::ghosts_update(unsigned data_parts) {
  ghost_communicator(&exchange_ghosts_comm, map_data_parts(data_parts));
}
void CellStructure::ghosts_reduce_forces() {
  ghost_communicator(&collect_ghost_force_comm, GHOSTTRANS_FORCE);
}
