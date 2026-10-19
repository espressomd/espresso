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

#include "cell_system/AtomDecomposition.hpp"

#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "particle_store/MigrationPack.hpp"

#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_to_all.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

void AtomDecomposition::configure_neighbors() {
  std::vector<Cell *> red_neighbors;
  std::vector<Cell *> black_neighbors;

  /* distribute force calculation work  */
  for (int n = 0; n < m_comm.size(); n++) {
    if (m_comm.rank() == n) {
      continue;
    }

    if (n < m_comm.rank()) {
      red_neighbors.push_back(&cells.at(n));
    } else {
      black_neighbors.push_back(&cells.at(n));
    }
  }

  local().m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
}

GhostComm::HaloPlan AtomDecomposition::make_halo_plan() {
  using GhostComm::CollectivePattern;
  using GhostComm::CollectiveSection;
  using GhostComm::HaloPlan;

  HaloPlan plan;
  plan.comm = m_comm;

  // Single rank: no communication needed; collective section is None.
  if (m_comm.size() == 1) {
    plan.collective = CollectiveSection{CollectivePattern::None, {}};
    return plan;
  }

  // One cell pointer per rank: cells[root] is the cell owned by that root.
  // The engine uses op.direction to pick Broadcast (Push) or ReduceSum
  // (Reduce) at run time, so we store Broadcast as the canonical marker that
  // this section is active.  run_collective reads op.direction to decide which
  // MPI collective to invoke.
  std::vector<Cell *> cell_ptrs;
  cell_ptrs.reserve(static_cast<std::size_t>(m_comm.size()));
  for (int n = 0; n < m_comm.size(); ++n) {
    cell_ptrs.push_back(std::addressof(cells.at(static_cast<std::size_t>(n))));
  }
  plan.collective =
      CollectiveSection{CollectivePattern::Broadcast, std::move(cell_ptrs)};
  return plan;
}

void AtomDecomposition::configure_comms() {
  m_halo_plan = make_halo_plan();
  // NOTE: validation is deferred to the constructor, AFTER mark_cells() has
  // populated local_cells()/ghost_cells(). Validating here would check empty
  // spans (vacuously) since mark_cells() runs later.
}

void AtomDecomposition::mark_cells() {
  m_local_cells.resize(1, std::addressof(local()));
  m_ghost_cells.clear();
  for (int n = 0; n < m_comm.size(); n++) {
    if (n != m_comm.rank()) {
      m_ghost_cells.push_back(std::addressof(cells.at(n)));
    }
  }
}

void AtomDecomposition::resort(bool global_flag,
                               std::vector<ParticleChange> &diff) {
  auto &store = local().store();
  // Fold positions of all committed rows (write-through views).
  for (auto &p : local().particles()) {
    Utils::Vector3d position = p.pos();
    Utils::Vector3i image_box = p.image_box();
    m_box.fold_position(position, image_box);
    p.pos() = position;
    p.image_box() = image_box;

    p.pos_at_last_verlet_update() = p.pos();
  }

  /* Local updates are a NoOp for this decomposition. */
  if (not global_flag) {
    return;
  }

  // A mis-owned particle is copied out of the live store into a staging-store
  // row (stage_row); its staging-row index goes into the target rank's per-rank
  // bucket, which is then packed per-field (MigrationPack::pack_rows) and
  // exchanged as a byte buffer via all_to_all.
  assert(m_migration_staging && "migration staging store not installed");
  auto &staging = *m_migration_staging.store;

  // Sort displaced particles into per-rank staging-row buckets (iterate the
  // committed rows by RAW position; drop_row marks the row pending-removed --
  // the range keeps its order, no swap-with-back -- so we advance
  // unconditionally).
  std::vector<std::vector<int>> send_rows(m_comm.size());
  auto const row_offset = local().offset();
  auto const row_count = local().count();
  for (std::size_t index = 0u; index < row_count; ++index) {
    auto const live_row = static_cast<int>(row_offset + index);
    // A pending-removed row is dead; never migrate it (see the same guard in
    // RegularDecomposition::resort). Ownership here is id-based and a removal
    // does not change the id, so this cannot currently fire -- it keeps the
    // three decompositions' resort loops on the same rule.
    if (store.is_pending_removal(live_row)) {
      continue;
    }
    auto const id = store.id(live_row);
    auto const target_node = id_to_rank(id);
    if (target_node != m_comm.rank()) {
      diff.emplace_back(RemovedParticle{id});
      send_rows.at(target_node)
          .push_back(m_migration_staging.stage_row(live_row));
      CellParticleStorage::drop_row(local(), index);
    }
  }

  // Pack each rank's bucket into a byte buffer and exchange.
  std::vector<std::vector<char>> send_buf(m_comm.size());
  for (int n = 0; n < m_comm.size(); ++n) {
    MigrationPack::pack_rows(staging, send_rows[static_cast<std::size_t>(n)],
                             send_buf[static_cast<std::size_t>(n)]);
  }
  std::vector<std::vector<char>> recv_buf(m_comm.size());
  boost::mpi::all_to_all(m_comm, send_buf, recv_buf);

  diff.emplace_back(ModifiedList{local()});

  // Unpack the received buffers (in rank order) into fresh staging rows, then
  // stage each staging row into the local cell as a row reference: the next
  // store rebuild copies it into a committed row. The staging store is NOT
  // cleared here -- the staged row references must remain valid until
  // CellStructure commits them (ensure_particle_store_synchronized), which then
  // resets the staging store.
  for (auto const &buffer : recv_buf) {
    if (buffer.empty()) {
      continue;
    }
    std::uint64_t count = 0u;
    std::memcpy(&count, buffer.data(), sizeof(count));
    if (count == 0u) {
      continue;
    }
    auto const first_row =
        m_migration_staging.reserve_rows(static_cast<int>(count));
    MigrationPack::unpack_rows(staging, first_row, buffer);
    for (int k = 0; k < static_cast<int>(count); ++k) {
      CellParticleStorage::insert_staged_row(local(), staging, first_row + k);
    }
  }
}

AtomDecomposition::AtomDecomposition(BoxGeometry const &box_geo)
    : m_box(box_geo) {}

AtomDecomposition::AtomDecomposition(boost::mpi::communicator comm,
                                     BoxGeometry const &box_geo)
    : m_comm(std::move(comm)), cells(m_comm.size()), m_box(box_geo) {
  /* create communicators */
  configure_comms();
  /* configure neighbor relations */
  configure_neighbors();
  /* fill local and ghost cell lists */
  mark_cells();
  /* classify local cells as interior or boundary.
   *
   * AtomDecomposition has no spatial locality: the single local cell
   * interacts with every other rank's cell and there is no subset of
   * particles whose force contributions are guaranteed to arrive before
   * the velocity update.  Interior is therefore always empty and all local
   * cells are boundary.  mark_boundary_cells() handles the ghost-neighbour
   * case (multi-rank); the explicit loop below catches the single-rank case
   * where there are no ghost cells and no neighbours at all.
   */
  GhostComm::mark_boundary_cells(AtomDecomposition::local_cells(),
                                 AtomDecomposition::ghost_cells());
  for (Cell *c : AtomDecomposition::local_cells()) {
    c->m_is_boundary = true;
  }
#ifdef ESPRESSO_ADDITIONAL_CHECKS
  // Validate now that local_cells()/ghost_cells() are populated by
  // mark_cells().
  assert(GhostComm::report_violations(
      GhostComm::validate_halo_plan(m_halo_plan,
                                    AtomDecomposition::local_cells(),
                                    AtomDecomposition::ghost_cells()),
      "AtomDecomposition"));
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

Utils::Vector3d AtomDecomposition::max_cutoff() const {
  return Utils::Vector3d::broadcast(std::numeric_limits<double>::infinity());
}

Utils::Vector3d AtomDecomposition::max_range() const { return max_cutoff(); }
