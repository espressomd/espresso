/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlan.hpp"

#include <functional>
#include <span>
#include <unordered_map>
#include <unordered_set>

namespace GhostComm {

/**
 * @brief Classify each local cell as interior or boundary.
 *
 * A local cell is *boundary* for two independent reasons:
 *
 *   Source 1 — geometric wrap rules: the cell has a ghost neighbor (rule a)
 *     or participates in a periodic wrap-around neighborship on a single MPI
 *     rank (rule b, captured by the optional @p wrap_predicate).  These rules
 *     serve future Euclidean-distance uses: wrap-neighborship must mark a cell
 *     boundary even when it is not exported by the current plan.
 *
 *   Source 2 — plan membership: the cell appears as a send source in any
 *     NeighborComm or as a LocalComm src.  Plans for geometries such as
 *     Lees-Edwards fully-connected boundaries and ELC periodicity-change paths
 *     can export cells that the geometric rules never see.  Those cells must
 *     still be boundary so that the integrator step-2 / force-reduce overlap
 *     does not update their velocities before remote force contributions
 *     arrive.
 *
 * Both sources are complementary: keep both to satisfy present and future
 * callers.  @ref mark_plan_cells_boundary applies source 2.
 *
 * All other local cells are *interior*.
 *
 * The function first resets every local cell to interior (idempotent on
 * repeated calls / plan rebuild), then marks boundary cells.
 *
 * Call this right after `m_halo_plan = make_halo_plan()` in each
 * decomposition, where `local_cells()` and `ghost_cells()` are already
 * populated.
 *
 * @param local_cells  Local cell pointer span.
 * @param ghost_cells  Ghost cell pointer span.
 * @param wrap_predicate  Optional predicate <tt>(Cell const *a, Cell const *b)
 *   -> bool</tt> that returns @c true when the neighbor relation a->b crosses
 *   a periodic box boundary (both cells are local).
 *   When omitted, no additional wrapping boundary is detected.
 */
inline void mark_boundary_cells(
    std::span<Cell *const> local_cells, std::span<Cell *const> ghost_cells,
    std::function<bool(Cell const *, Cell const *)> wrap_predicate = nullptr) {
  // Build a set of ghost ParticleList pointers for O(1) lookup.
  std::unordered_set<ParticleList const *> ghost_set;
  ghost_set.reserve(ghost_cells.size());
  for (Cell *c : ghost_cells) {
    ghost_set.insert(&c->particles());
  }

  // Reset all local cells to interior first (idempotent on rebuild).
  for (Cell *c : local_cells) {
    c->m_is_boundary = false;
  }

  // Mark a cell as boundary iff any of its neighbors is a ghost (rule a) or
  // the wrap predicate fires for that neighbor pair (rule b).
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      if (ghost_set.count(&n->particles()) ||
          (wrap_predicate && wrap_predicate(c, n))) {
        c->m_is_boundary = true;
        break;
      }
    }
  }
}

/**
 * @brief Mark plan-exported local cells as boundary (source 2, see
 * @ref mark_boundary_cells).
 *
 * Iterates the halo plan and marks every local cell that appears as a send
 * source.  Specifically: every `NeighborComm::send[k].cell` and every
 * `LocalComm::src`.  Recv/dst targets are ghost cells by construction and
 * are not touched.
 *
 * Call this immediately after @ref mark_boundary_cells so that the combined
 * invariant "interior ⇒ not exported by the plan" holds for any plan shape.
 *
 * @param plan        The halo plan produced by `make_halo_plan()`.
 * @param local_cells Local cell pointer span (used to build the ParticleList ->
 *                    Cell reverse map in O(n)).
 */
inline void mark_plan_cells_boundary(HaloPlan const &plan,
                                     std::span<Cell *const> local_cells) {
  // Build a reverse map: ParticleList* -> Cell* for O(1) lookup.
  std::unordered_map<ParticleList const *, Cell *> pl_to_cell;
  pl_to_cell.reserve(local_cells.size());
  for (Cell *c : local_cells) {
    pl_to_cell[&c->particles()] = c;
  }

  // Mark every NeighborComm send source boundary.
  for (auto const &nc : plan.neighbors) {
    for (auto const &sr : nc.send) {
      auto it = pl_to_cell.find(sr.cell);
      if (it != pl_to_cell.end()) {
        it->second->m_is_boundary = true;
      }
    }
  }

  // Mark every LocalComm src boundary.
  for (auto const &lc : plan.local) {
    auto it = pl_to_cell.find(lc.src);
    if (it != pl_to_cell.end()) {
      it->second->m_is_boundary = true;
    }
  }
}

} // namespace GhostComm
