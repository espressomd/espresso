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
#include "cell_system/ParticleListOperations.hpp"

#include "particle_store/ParticleStore.hpp"

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

  /* Check for n_square type particles in regular decomposition. Iterate
   * committed rows by raw position and, for a misplaced particle, copy its live
   * row into the shared staging store (stage_row), mark the source row
   * pending-removed via drop_row (no swap-with-back; the raw index stays valid
   * so we advance unconditionally), and stage a reference to the staging row in
   * the target child cell (insert_staged_row). The next store rebuild --
   * triggered by m_commit_store below, BEFORE the child resorts -- emits
   * surviving rows then staged rows in permutation order, drops pending-removed
   * rows, and resets the staging store. */
  assert(m_migration_staging && "migration staging store not installed");
  auto &staging = *m_migration_staging.store;
  for (auto &cell_rd : m_regular_decomposition.local_cells()) {
    auto &store = cell_rd->store();
    // Iterate committed rows by RAW position: drop_row marks the row
    // pending-removed (no swap-with-back), so advance unconditionally.
    auto const rd_offset = cell_rd->offset();
    auto const rd_count = cell_rd->count();
    for (std::size_t index = 0u; index < rd_count; ++index) {
      auto const live_row = static_cast<int>(rd_offset + index);
      /* Skip a row a prior pass already dropped (marked pending-removed by
       * drop_row). The n_square-scan loop below is nested in this regular-cell
       * loop, so it re-scans every n_square cell once per regular cell; a row
       * already migrated in an earlier pass is invisible to live iteration and
       * must not be re-staged, or it would be duplicated once per outer pass.
       */
      if (store.is_pending_removal(live_row)) {
        continue;
      }
      auto const type = store.type(live_row);
      /* Particle is in the right decomposition, i.e. has no n_square type */
      if (not is_n_square_type(type)) {
        continue;
      }

      /* else remove from current cell ... */
      auto const staging_row = m_migration_staging.stage_row(live_row);
      CellParticleStorage::drop_row(*cell_rd, index);
      diff.emplace_back(ModifiedList{*cell_rd});
      diff.emplace_back(RemovedParticle{staging.id(staging_row)});

      /* ... and insert into a n_square cell */
      auto const first_local_cell = m_n_square.get_local_cells()[0];
      CellParticleStorage::insert_staged_row(*first_local_cell, staging,
                                             staging_row);
      diff.emplace_back(ModifiedList{*first_local_cell});
    }

    /* Now check for regular decomposition type particles in n_square */
    for (auto &cell_ns : m_n_square.local_cells()) {
      auto &store = cell_ns->store();
      // Iterate committed rows by RAW position: drop_row marks the row
      // pending-removed (no swap-with-back), so advance unconditionally.
      auto const ns_offset = cell_ns->offset();
      auto const ns_count = cell_ns->count();
      for (std::size_t index = 0u; index < ns_count; ++index) {
        auto const live_row = static_cast<int>(ns_offset + index);
        /* Skip a row a prior pass already dropped (marked pending-removed by
         * drop_row). This loop is re-entered once per regular cell (it is
         * nested in the regular-cell loop above), so an rd-type particle
         * already migrated out of this n_square cell in an earlier pass is
         * invisible to live iteration and must not be re-staged, or it would be
         * duplicated once per outer pass. */
        if (store.is_pending_removal(live_row)) {
          continue;
        }
        auto const type = store.type(live_row);
        /* Particle is of n_square type */
        if (is_n_square_type(type)) {
          continue;
        }

        /* else remove from current cell ... */
        auto const staging_row = m_migration_staging.stage_row(live_row);
        CellParticleStorage::drop_row(*cell_ns, index);
        diff.emplace_back(ModifiedList{*cell_ns});
        diff.emplace_back(RemovedParticle{staging.id(staging_row)});

        /* ... and insert in regular decomposition. The home cell is decided
         * from the staged row's position (read through a view over it). */
        auto const target_cell =
            particle_to_cell(staging.make_view(staging_row));
        /* if particle belongs to this node insert it into correct cell */
        if (target_cell != nullptr) {
          CellParticleStorage::insert_staged_row(*target_cell, staging,
                                                 staging_row);
          diff.emplace_back(ModifiedList{*target_cell});
        }
        /* otherwise just put into regular decomposition */
        else {
          auto first_local_cell = m_regular_decomposition.get_local_cells()[0];
          CellParticleStorage::insert_staged_row(*first_local_cell, staging,
                                                 staging_row);
          diff.emplace_back(ModifiedList{*first_local_cell});
        }
      }
    }
  }

  /* The type-based moves above staged staging-row references into their target
   * child cells (insert_staged_row). Commit now so the child resorts below
   * iterate the correct committed cell contents -- otherwise a particle moved
   * into a cell would be invisible to that cell's own resort, changing the
   * final placement/order. The commit copies each staged row into a committed
   * row and resets the shared staging store, so the child resorts start with a
   * clean staging store. */
  if (m_commit_store) {
    m_commit_store();
  }

  /* now resort into correct cells within the respective decompositions */
  m_regular_decomposition.resort(global, diff);
  m_n_square.resort(global, diff);

  /* The child resorts staged migrated/new particles into cells but did not
   * commit them to store rows. Commit now so the internal ghost communications
   * below see committed rows/columns (the PARTNUM step still uses
   * Cell::size() = rows+staged for downstream ghost layers, but the DATA step
   * reads committed views). */
  if (m_commit_store) {
    m_commit_store();
  }

  GhostComm::halo_exchange(
      m_halo_plan, m_box, GHOSTTRANS_PARTNUM,
      {GhostComm::Direction::Push, GhostComm::Combine::Overwrite});
  /* Committing the just-staged ghosts before the DATA transfer, same reason. */
  if (m_commit_store) {
    m_commit_store();
  }
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
