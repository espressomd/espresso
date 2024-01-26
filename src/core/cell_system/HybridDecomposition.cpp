/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "cell_system/HybridDecomposition.hpp"

#include "cell_system/Cell.hpp"
#include "cell_system/CellStructure.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/sendrecv.hpp>

#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <set>
#include <utility>

HybridDecomposition::HybridDecomposition(boost::mpi::communicator comm,
                                         double cutoff_regular, double skin,
                                         std::function<bool()> get_ghost_flags,
                                         BoxGeometry const &box_geo,
                                         LocalBox const &local_box,
                                         std::set<int> n_square_types)
    : m_comm(std::move(comm)), m_box(box_geo), m_cutoff_regular(cutoff_regular),
      m_regular_decomposition(RegularDecomposition(
          m_comm, cutoff_regular + skin, m_box, local_box)),
      m_n_square(AtomDecomposition(m_comm, m_box)),
      m_n_square_types(std::move(n_square_types)),
      m_get_global_ghost_flags(std::move(get_ghost_flags)) {

  /* Vector containing cells of both child decompositions */
  m_local_cells = m_regular_decomposition.get_local_cells();
  auto local_cells_n_square = m_n_square.get_local_cells();
  std::copy(local_cells_n_square.begin(), local_cells_n_square.end(),
            std::back_inserter(m_local_cells));

  /* Vector containing ghost cells of both child decompositions */
  m_ghost_cells = m_regular_decomposition.get_ghost_cells();
  auto ghost_cells_n_square = m_n_square.get_ghost_cells();
  std::copy(ghost_cells_n_square.begin(), ghost_cells_n_square.end(),
            std::back_inserter(m_ghost_cells));

  /* Communicators that contain communications of both child decompositions */
  m_exchange_ghosts_comm = m_regular_decomposition.exchange_ghosts_comm();
  auto exchange_ghosts_comm_n_square = m_n_square.exchange_ghosts_comm();
  std::copy(exchange_ghosts_comm_n_square.communications.begin(),
            exchange_ghosts_comm_n_square.communications.end(),
            std::back_inserter(m_exchange_ghosts_comm.communications));

  m_collect_ghost_force_comm =
      m_regular_decomposition.collect_ghost_force_comm();
  auto collect_ghost_force_comm_n_square =
      m_n_square.collect_ghost_force_comm();
  std::copy(collect_ghost_force_comm_n_square.communications.begin(),
            collect_ghost_force_comm_n_square.communications.end(),
            std::back_inserter(m_collect_ghost_force_comm.communications));

  /* coupling between the child decompositions via neighborship relation */
  std::vector<Cell *> additional_reds = m_n_square.get_local_cells();
  std::copy(ghost_cells_n_square.begin(), ghost_cells_n_square.end(),
            std::back_inserter(additional_reds));
  for (auto &local_cell : m_regular_decomposition.local_cells()) {
    std::vector<Cell *> red_neighbors(local_cell->m_neighbors.red().begin(),
                                      local_cell->m_neighbors.red().end());
    std::vector<Cell *> black_neighbors(local_cell->m_neighbors.black().begin(),
                                        local_cell->m_neighbors.black().end());
    std::copy(additional_reds.begin(), additional_reds.end(),
              std::back_inserter(red_neighbors));
    local_cell->m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
  }
}

void HybridDecomposition::resort(bool global,
                                 std::vector<ParticleChange> &diff) {
  ParticleList displaced_parts;

  /* Check for n_square type particles in regular decomposition */
  for (auto &c : m_regular_decomposition.local_cells()) {
    for (auto it = c->particles().begin(); it != c->particles().end();) {
      /* Particle is in the right decomposition, i.e. has no n_square type */
      if (not is_n_square_type(it->type())) {
        std::advance(it, 1);
        continue;
      }

      /* else remove from current cell ... */
      auto p = std::move(*it);
      it = c->particles().erase(it);
      diff.emplace_back(ModifiedList{c->particles()});
      diff.emplace_back(RemovedParticle{p.id()});

      /* ... and insert into a n_square cell */
      auto const first_local_cell = m_n_square.get_local_cells()[0];
      first_local_cell->particles().insert(std::move(p));
      diff.emplace_back(ModifiedList{first_local_cell->particles()});
    }

    /* Now check for regular decomposition type particles in n_square */
    for (auto &c : m_n_square.local_cells()) {
      for (auto it = c->particles().begin(); it != c->particles().end();) {
        /* Particle is of n_square type */
        if (is_n_square_type(it->type())) {
          std::advance(it, 1);
          continue;
        }

        /* else remove from current cell ... */
        auto p = std::move(*it);
        it = c->particles().erase(it);
        diff.emplace_back(ModifiedList{c->particles()});
        diff.emplace_back(RemovedParticle{p.id()});

        /* ... and insert in regular decomposition */
        auto const target_cell = particle_to_cell(p);
        /* if particle belongs to this node insert it into correct cell */
        if (target_cell != nullptr) {
          target_cell->particles().insert(std::move(p));
          diff.emplace_back(ModifiedList{target_cell->particles()});
        }
        /* otherwise just put into regular decomposition */
        else {
          auto first_local_cell = m_regular_decomposition.get_local_cells()[0];
          first_local_cell->particles().insert(std::move(p));
          diff.emplace_back(ModifiedList{first_local_cell->particles()});
        }
      }
    }
  }

  /* now resort into correct cells within the respective decompositions */
  m_regular_decomposition.resort(global, diff);
  m_n_square.resort(global, diff);

  /* basically do CellStructure::ghost_count() */
  ghost_communicator(exchange_ghosts_comm(), m_box, GHOSTTRANS_PARTNUM);

  /* basically do CellStructure::ghost_update(unsigned data_parts) */
  ghost_communicator(exchange_ghosts_comm(), m_box,
                     map_data_parts(m_get_global_ghost_flags()));
}

std::size_t HybridDecomposition::count_particles(
    std::vector<Cell *> const &local_cells) const {
  std::size_t count_local = 0;
  std::size_t count_global = 0;
  for (auto const &cell : local_cells) {
    count_local += cell->particles().size();
  }
  boost::mpi::reduce(m_comm, count_local, count_global, std::plus<>{}, 0);
  return count_global;
}
