/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include "ParticleList.hpp"
#include "ghosts.hpp"
#include "ghosts/HaloExchange.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/sendrecv.hpp>

#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
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
          m_comm, cutoff_regular + skin, m_box, local_box, std::nullopt)),
      m_n_square(AtomDecomposition(m_comm, m_box)),
      m_n_square_types(std::move(n_square_types)),
      m_get_global_ghost_flags(std::move(get_ghost_flags)) {

  /* Vector containing cells of both child decompositions */
  m_local_cells = m_regular_decomposition.get_local_cells();
  auto local_cells_n_square = m_n_square.get_local_cells();
  std::ranges::copy(local_cells_n_square, std::back_inserter(m_local_cells));

  /* Vector containing ghost cells of both child decompositions */
  m_ghost_cells = m_regular_decomposition.get_ghost_cells();
  auto ghost_cells_n_square = m_n_square.get_ghost_cells();
  std::ranges::copy(ghost_cells_n_square, std::back_inserter(m_ghost_cells));

  /* coupling between the child decompositions via neighborship relation */
  std::vector<Cell *> additional_reds = m_n_square.get_local_cells();
  std::ranges::copy(ghost_cells_n_square, std::back_inserter(additional_reds));
  for (auto &local_cell : m_regular_decomposition.local_cells()) {
    std::vector<Cell *> red_neighbors(local_cell->m_neighbors.red().begin(),
                                      local_cell->m_neighbors.red().end());
    std::vector<Cell *> black_neighbors(local_cell->m_neighbors.black().begin(),
                                        local_cell->m_neighbors.black().end());
    std::ranges::copy(additional_reds, std::back_inserter(red_neighbors));
    local_cell->m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
  }

  m_halo_plan = make_halo_plan();

  /* classify local cells as interior or boundary.
   *
   * HybridDecomposition combines a RegularDecomposition (which already has
   * its own wrap-aware classification) with an AtomDecomposition n-square
   * child.  The combined cell set extends every regular local cell's
   * neighbor list with the n-square cells, so regular cells that were
   * interior under RegularDecomposition alone may now interact with
   * n-square ghost cells.  Determining which cells remain genuinely
   * interior after the coupling is non-trivial; we conservatively mark all
   * local cells as boundary.  The compute/comm overlap degenerates
   * gracefully to a no-op interior pass — identical to today's blocking
   * path.
   */
  GhostComm::mark_boundary_cells(HybridDecomposition::local_cells(),
                                 HybridDecomposition::ghost_cells());
  for (Cell *c : HybridDecomposition::local_cells()) {
    c->m_is_boundary = true;
  }
#ifdef ESPRESSO_ADDITIONAL_CHECKS
  assert(GhostComm::report_violations(
      GhostComm::validate_halo_plan(m_halo_plan,
                                    HybridDecomposition::local_cells(),
                                    HybridDecomposition::ghost_cells()),
      "HybridDecomposition"));
  // NOTE: validate_halo_plan_symmetry is NOT called here.
  // During checkpoint loading, decompositions are transiently rebuilt while
  // maximal_cutoff is rank-divergent (ranks may have different cell grids for a
  // brief window before the next consistent rebuild).  The transient plan is
  // never used — it is immediately replaced — so the asymmetry is harmless.
  // A construction-time collective all_to_all inside a ctor is also dangerous:
  // if one rank aborts the others block forever in the collective.
  // Symmetry is instead validated at FIRST USE of the plan in
  // halo_exchange_start (see GhostComm::halo_exchange_start in
  // HaloExchange.cpp).
#endif
}

GhostComm::HaloPlan HybridDecomposition::make_halo_plan() {
  // Build the combined plan: p2p neighbors and local copies come from the
  // regular child; the collective section comes from the n-square child.
  GhostComm::HaloPlan plan;
  plan.comm = m_comm; // use the world comm for HybridDecomposition's plan

  auto const *regular_plan = m_regular_decomposition.halo_plan();
  if (regular_plan) {
    plan.neighbors = regular_plan->neighbors;
    plan.local = regular_plan->local;
  }

  // Overlay the collective section from the n-square child.
  auto const *nsq_plan = m_n_square.halo_plan();
  if (nsq_plan && nsq_plan->collective) {
    plan.collective = nsq_plan->collective;
  }

  return plan;
}

void HybridDecomposition::resort(bool global,
                                 std::vector<ParticleChange> &diff) {
  ParticleList displaced_parts;

  /* Check for n_square type particles in regular decomposition */
  for (auto &cell_rd : m_regular_decomposition.local_cells()) {
    for (auto it = cell_rd->particles().begin();
         it != cell_rd->particles().end();) {
      /* Particle is in the right decomposition, i.e. has no n_square type */
      if (not is_n_square_type(it->type())) {
        std::advance(it, 1);
        continue;
      }

      /* else remove from current cell ... */
      auto p = std::move(*it);
      it = cell_rd->particles().erase(it);
      diff.emplace_back(ModifiedList{cell_rd->particles()});
      diff.emplace_back(RemovedParticle{p.id()});

      /* ... and insert into a n_square cell */
      auto const first_local_cell = m_n_square.get_local_cells()[0];
      first_local_cell->particles().insert(std::move(p));
      diff.emplace_back(ModifiedList{first_local_cell->particles()});
    }

    /* Now check for regular decomposition type particles in n_square */
    for (auto &cell_ns : m_n_square.local_cells()) {
      for (auto it = cell_ns->particles().begin();
           it != cell_ns->particles().end();) {
        /* Particle is of n_square type */
        if (is_n_square_type(it->type())) {
          std::advance(it, 1);
          continue;
        }

        /* else remove from current cell ... */
        auto p = std::move(*it);
        it = cell_ns->particles().erase(it);
        diff.emplace_back(ModifiedList{cell_ns->particles()});
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

  GhostComm::halo_exchange(
      m_halo_plan, m_box, GHOSTTRANS_PARTNUM,
      {GhostComm::Direction::Push, GhostComm::Combine::Overwrite});
  GhostComm::halo_exchange(
      m_halo_plan, m_box, map_data_parts(m_get_global_ghost_flags()),
      {GhostComm::Direction::Push, GhostComm::Combine::Overwrite});
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
