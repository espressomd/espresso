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
/** \file
 *  Local (single-rank) validator for GhostComm::HaloPlan.
 *  See HaloPlanValidator.hpp.
 */

#include "ghosts/HaloPlanValidator.hpp"

#include <boost/mpi/collectives/all_to_all.hpp>

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace GhostComm {

std::vector<std::string>
validate_halo_plan(HaloPlan const &plan, std::span<Cell *const> local_cells,
                   std::span<Cell *const> ghost_cells) {
  std::vector<std::string> violations;

  // Build the set of expected ghost ParticleList pointers.
  std::unordered_set<ParticleList const *> ghost_set;
  for (Cell *c : ghost_cells) {
    ghost_set.insert(&c->particles());
  }

  // Collective-covered cells: AtomDecomposition (collective-only) and
  // HybridDecomposition (neighbors + collective) fill some/all ghosts via the
  // n-square broadcast/reduce section, NOT via neighbors/local. The section's
  // `cells` vector identifies exactly which ParticleLists it covers (one per
  // rank; the entry for the local rank is this rank's own cell, the rest are
  // ghost copies of every other rank's owned cell). Ghost cells found here are
  // covered even though they are not point-to-point recv/dst targets.
  bool const has_collective =
      plan.collective.has_value() &&
      plan.collective->pattern != CollectivePattern::None;
  std::unordered_set<ParticleList const *> collective_set;
  if (has_collective) {
    for (ParticleList const *pl : plan.collective->cells) {
      collective_set.insert(pl);
    }
  }

  // Set of ghost cells that some local cell actually interacts with, i.e. that
  // appear in a local cell's neighbor stencil. Only these ghosts must be
  // filled by the plan; a ghost that no local cell ever references carries no
  // physics and needs no communication. This happens on a single MPI rank:
  // make_halo_plan() returns an empty plan there because init_cell_interactions
  // wires periodic neighbourships directly between local cells (see
  // RegularDecomposition), so the halo-layer cells that mark_cells() still
  // classifies as ghosts are never referenced and must not be flagged as
  // uncovered. Double-fill and out-of-ghost-set targets stay strict below, and
  // the neighborship-match check still requires every *referenced* ghost to be
  // covered, so this cannot mask a real missing-communication defect.
  std::unordered_set<ParticleList const *> referenced_ghosts;
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      ParticleList const *pl = &n->particles();
      if (ghost_set.contains(pl)) {
        referenced_ghosts.insert(pl);
      }
    }
  }

  // Check peer-uniqueness and shape; accumulate recv/dst fill counts.
  std::unordered_map<ParticleList const *, int> fill_count;
  std::unordered_set<int> seen_peers;

  for (auto const &nc : plan.neighbors) {
    // Peer-uniqueness
    if (!seen_peers.insert(nc.peer).second) {
      std::ostringstream oss;
      oss << "peer " << nc.peer << " appears in more than one NeighborComm";
      violations.emplace_back(oss.str());
    }

    // Shape: send.size() == recv.size()
    if (nc.send.size() != nc.recv.size()) {
      std::ostringstream oss;
      oss << "NeighborComm peer=" << nc.peer
          << " has send.size()=" << nc.send.size()
          << " != recv.size()=" << nc.recv.size();
      violations.emplace_back(oss.str());
    }

    // Accumulate recv fill counts; check targets are in ghost set.
    for (ParticleList const *pl : nc.recv) {
      if (not ghost_set.contains(pl)) {
        std::ostringstream oss;
        oss << "NeighborComm peer=" << nc.peer
            << " recv target is not a ghost cell";
        violations.emplace_back(oss.str());
      }
      ++fill_count[pl];
    }
  }

  // Accumulate local.dst fill counts; check targets are in ghost set.
  for (auto const &lc : plan.local) {
    ParticleList const *pl = lc.dst;
    if (not ghost_set.contains(pl)) {
      violations.emplace_back("LocalComm dst target is not a ghost cell");
    }
    ++fill_count[pl];
  }

  // Coverage: every ghost that a local cell references must be filled.
  // Point-to-point ghosts must appear exactly once as a recv/dst target; a
  // double-fill is always a defect. Ghosts covered by the collective section
  // are filled by the broadcast/reduce instead, so a zero point-to-point
  // fill-count is fine for them (and they do not contribute to fill_count, so
  // they cannot trigger the double-fill path). A ghost that no local cell
  // references (e.g. the halo-layer cells on a single MPI rank, where the plan
  // is intentionally empty) carries no physics and is not required to be
  // filled, so a zero fill-count is fine for it too.
  for (Cell *c : ghost_cells) {
    ParticleList const *pl = &c->particles();
    auto it = fill_count.find(pl);
    int count = (it != fill_count.end()) ? it->second : 0;
    if (count == 0) {
      if (not collective_set.contains(pl) and referenced_ghosts.contains(pl)) {
        violations.emplace_back(
            "ghost cell is never filled (missing recv/dst)");
      }
    } else if (count > 1) {
      std::ostringstream oss;
      oss << "ghost cell is filled " << count << " times (expected 1)";
      violations.emplace_back(oss.str());
    }
  }

  // Neighborship-match: every local cell's ghost neighbor must be covered,
  // either as a point-to-point recv/dst target or by the collective section.
  for (Cell *c : local_cells) {
    for (Cell *n : c->neighbors().all()) {
      ParticleList const *pl = &n->particles();
      if (not ghost_set.contains(pl)) {
        continue; // not a ghost neighbor
      }
      if (collective_set.contains(pl)) {
        continue; // covered by the collective broadcast/reduce section
      }
      auto it = fill_count.find(pl);
      int count = (it != fill_count.end()) ? it->second : 0;
      if (count == 0) {
        violations.emplace_back(
            "local cell has ghost neighbor that is not a covered recv/dst "
            "target (referenced-but-uncommunicated ghost)");
      }
    }
  }

  // Interior/boundary consistency: a local cell marked interior
  // (is_boundary()==false) must have no ghost neighbor.
  for (Cell *c : local_cells) {
    if (c->is_boundary()) {
      continue; // boundary cells are expected to have ghost neighbors
    }
    for (Cell *n : c->neighbors().all()) {
      if (ghost_set.contains(&n->particles())) {
        violations.emplace_back("interior cell has a ghost neighbor");
        break; // one violation per cell is sufficient
      }
    }
  }

  // Overlap-safety invariant (Task 5.3 precondition): an interior cell must
  // not appear as a send source in any NeighborComm, and must not appear as
  // a LocalComm src.  Send sources are real cells whose data gets copied to a
  // peer's ghost; LocalComm srcs are the matching self-copy sources.  If an
  // interior cell were in either, the force-reduce step would overwrite its
  // ghost copy after the interior velocity update — violating the overlap
  // correctness assumption.
  //
  // Build a set of ParticleList* that belong to interior local cells.
  std::unordered_set<ParticleList const *> interior_set;
  for (Cell *c : local_cells) {
    if (!c->is_boundary()) {
      interior_set.insert(&c->particles());
    }
  }
  if (!interior_set.empty()) {
    for (auto const &nc : plan.neighbors) {
      for (auto const &sr : nc.send) {
        if (interior_set.count(sr.cell)) {
          std::ostringstream oss;
          oss << "interior cell appears as NeighborComm send source (peer="
              << nc.peer << ") — overlap-safety invariant violated";
          violations.emplace_back(oss.str());
        }
      }
    }
    for (auto const &lc : plan.local) {
      if (interior_set.count(lc.src)) {
        violations.emplace_back(
            "interior cell appears as LocalComm src — overlap-safety invariant "
            "violated");
      }
    }
  }

  return violations;
}

std::vector<std::string> validate_halo_plan_symmetry(HaloPlan const &plan) {
  std::vector<std::string> violations;

  int const n = plan.comm.size();
  int const me = plan.comm.rank();

  // Build per-rank send and recv count arrays (size n, default 0).
  std::vector<int> my_send_to(n, 0);
  std::vector<int> my_recv_from(n, 0);
  for (auto const &nc : plan.neighbors) {
    my_send_to[nc.peer] = static_cast<int>(nc.send.size());
    my_recv_from[nc.peer] = static_cast<int>(nc.recv.size());
  }

  // Collective all-to-all: element j of my_send_to goes to rank j.
  // After this call peers_send_to_me[j] == rank j's send count toward me.
  std::vector<int> peers_send_to_me(n, 0);
  boost::mpi::all_to_all(plan.comm, my_send_to, peers_send_to_me);

  // Check invariant: what I expect to receive from j == what j sends me.
  for (int j = 0; j < n; ++j) {
    if (j == me)
      continue;
    if (my_recv_from[j] != peers_send_to_me[j]) {
      std::ostringstream oss;
      oss << "symmetry mismatch with peer " << j << ": I expect to recv "
          << my_recv_from[j] << " items but peer sends " << peers_send_to_me[j]
          << " items";
      violations.emplace_back(oss.str());
    }
  }

  return violations;
}

bool report_violations(std::vector<std::string> const &violations,
                       char const *context) {
  for (auto const &v : violations) {
    std::cerr << "halo plan validation (" << context << "): " << v << "\n";
  }
  return violations.empty();
}

} // namespace GhostComm
