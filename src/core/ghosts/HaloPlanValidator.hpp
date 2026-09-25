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

#include <span>
#include <string>
#include <vector>

namespace GhostComm {

/**
 * @brief Validate a HaloPlan for correctness.
 *
 * Checks:
 * 1. Coverage: every ghost cell that a local cell references is filled.
 *    Point-to-point ghosts appear as a recv/local.dst target exactly once; no
 *    p2p target lies outside the ghost set; none appears twice. Ghosts filled
 *    by the collective broadcast/reduce section (AtomDecomposition,
 *    HybridDecomposition) are covered by that section instead of a recv/dst
 *    target. A ghost that no local cell references (e.g. the halo-layer cells
 *    on a single MPI rank, where the plan is intentionally empty) carries no
 *    physics and is not required to be filled.
 * 2. Neighborship-match: every local cell's ghost neighbor is covered, either
 *    as a recv/dst target or by the collective section.
 * 3. Peer-uniqueness: each peer value appears in at most one NeighborComm.
 * 4. Shape: each NeighborComm has send.size() == recv.size().
 * 5. Interior/boundary consistency: a cell marked interior has no ghost
 *    neighbor (enforced by mark_boundary_cells()).
 * 6. Overlap-safety invariant: an interior cell must not appear as a
 *    NeighborComm send source or a LocalComm src.  This is the exact
 *    precondition that the integrator step-2 / force-reduce overlap (Task
 *    5.3) relies on: interior cells receive no reduce contributions and
 *    therefore do not require the reduce to complete before their velocity
 *    is updated.  LocalComm.dst and NeighborComm.recv are ghost cells by
 *    construction (guaranteed by the dst-in-ghost-set coverage check #1),
 *    so they are not checked here.
 *
 * @returns a human-readable violation string per problem; empty = valid.
 */
std::vector<std::string> validate_halo_plan(HaloPlan const &plan,
                                            std::span<Cell *const> local_cells,
                                            std::span<Cell *const> ghost_cells);

/**
 * @brief Cross-rank symmetry check for a HaloPlan.
 *
 * Uses a collective all-to-all exchange of per-rank send-counts so that every
 * rank participates regardless of whether the neighbor sets are symmetric.
 * This avoids the deadlock that a naive per-neighbor isend/irecv pattern
 * would cause when a rank sends to a peer that has no matching recv posted.
 *
 * Invariant: for every peer P, my recv-count from P must equal P's send-count
 * to me (i.e. my_recv_from[P] == peers_send_to_me[P]).
 *
 * @returns a human-readable violation string per mismatch; empty = valid.
 */
std::vector<std::string> validate_halo_plan_symmetry(HaloPlan const &plan);

/**
 * @brief Print violations to stderr and return whether the list was empty.
 *
 * Intended for use inside assert(): a bare
 * assert(validate_halo_plan(...).empty()) aborts without showing WHICH
 * invariant failed, which makes CI failures undebuggable.  Wrap the call:
 * assert(GhostComm::report_violations(validate_halo_plan(...), "context")).
 */
bool report_violations(std::vector<std::string> const &violations,
                       char const *context);

} // namespace GhostComm
