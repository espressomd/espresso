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

#include "AtomDecomposition.hpp"
#include "DomainDecomposition.hpp"

#include <utils/contains.hpp>

#include <algorithm>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

void CellStructure::check_particle_index() {
  auto const max_id = get_max_local_particle_id();

  for (auto const &p : local_particles()) {
    auto const id = p.identity();

    if (id < 0 || id > max_id) {
      throw std::runtime_error("Particle id out of bounds.");
    }

    if (get_local_particle(id) != &p) {
      throw std::runtime_error("Invalid local particle index entry.");
    }
  }

  /* checks: local particle id */
  int local_part_cnt = 0;
  for (int n = 0; n < get_max_local_particle_id() + 1; n++) {
    if (get_local_particle(n) != nullptr) {
      local_part_cnt++;
      if (get_local_particle(n)->p.identity != n) {
        throw std::runtime_error("local_particles part has corrupted id.");
      }
    }
  }

  if (local_part_cnt != local_particles().size()) {
    throw std::runtime_error(
        std::to_string(local_particles().size()) + " parts in cells but " +
        std::to_string(local_part_cnt) + " parts in local_particles");
  }
}

void CellStructure::check_particle_sorting() {
  for (auto cell : local_cells()) {
    for (auto const &p : cell->particles()) {
      if (particle_to_cell(p) != cell) {
        throw std::runtime_error("misplaced particle with id " +
                                 std::to_string(p.identity()));
      }
    }
  }
}

Cell *CellStructure::particle_to_cell(const Particle &p) {
  return decomposition().particle_to_cell(p);
}

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

  for (auto c : decomposition().local_cells()) {
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
  for (auto c : decomposition().local_cells()) {
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
#ifdef BOND_CONSTRAINT
         | ((DATA_PART_RATTLE & data_parts) ? GHOSTTRANS_RATTLE : 0u)
#endif
         | ((DATA_PART_BONDS & data_parts) ? GHOSTTRANS_BONDS : 0u);
  /* clang-format on */
}

void CellStructure::ghosts_count() {
  ghost_communicator(decomposition().exchange_ghosts_comm(),
                     GHOSTTRANS_PARTNUM);
}
void CellStructure::ghosts_update(unsigned data_parts) {
  ghost_communicator(decomposition().exchange_ghosts_comm(),
                     map_data_parts(data_parts));
}
void CellStructure::ghosts_reduce_forces() {
  ghost_communicator(decomposition().collect_ghost_force_comm(),
                     GHOSTTRANS_FORCE);
}
#ifdef BOND_CONSTRAINT
void CellStructure::ghosts_reduce_rattle_correction() {
  ghost_communicator(decomposition().collect_ghost_force_comm(),
                     GHOSTTRANS_RATTLE);
}
#endif

Utils::Span<Cell *> CellStructure::local_cells() {
  return decomposition().local_cells();
}

ParticleRange CellStructure::local_particles() {
  return Cells::particles(decomposition().local_cells());
}

ParticleRange CellStructure::ghost_particles() {
  return Cells::particles(decomposition().ghost_cells());
}

Utils::Vector3d CellStructure::max_cutoff() const {
  return decomposition().max_cutoff();
}

Utils::Vector3d CellStructure::max_range() const {
  return decomposition().max_range();
}

namespace {
/**
 * @brief Apply a @ref ParticleChange to a particle index.
 */
struct UpdateParticleIndexVisitor {
  CellStructure *cs;

  void operator()(RemovedParticle rp) const {
    cs->update_particle_index(rp.id, nullptr);
  }
  void operator()(ModifiedList mp) const { cs->update_particle_index(mp.pl); }
};
} // namespace

void CellStructure::resort_particles(int global_flag) {
  invalidate_ghosts();

  static std::vector<ParticleChange> diff;
  diff.clear();

  m_decomposition->resort(global_flag, diff);

  for (auto d : diff) {
    boost::apply_visitor(UpdateParticleIndexVisitor{this}, d);
  }

  m_rebuild_verlet_list = true;

#ifdef ADDITIONAL_CHECKS
  check_particle_index();
  check_particle_sorting();
#endif
}

void CellStructure::set_atom_decomposition(boost::mpi::communicator const &comm,
                                           BoxGeometry const &box) {
  set_particle_decomposition(std::make_unique<AtomDecomposition>(comm, box));
  m_type = CELL_STRUCTURE_NSQUARE;
}

void CellStructure::set_domain_decomposition(
    boost::mpi::communicator const &comm, double range, BoxGeometry const &box,
    LocalBox<double> const &local_geo) {
  set_particle_decomposition(
      std::make_unique<DomainDecomposition>(comm, range, box, local_geo));
  m_type = CELL_STRUCTURE_DOMDEC;
}
