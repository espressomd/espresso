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

#include "ghosts/HaloPlan.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_to_all.hpp>

#include <cassert>
#include <cstddef>
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

  // One cell pointer per rank: cells[root] is the ParticleList for that root.
  // The engine uses op.direction to pick Broadcast (Push) or ReduceSum
  // (Reduce) at run time, so we store Broadcast as the canonical marker that
  // this section is active.  run_collective reads op.direction to decide which
  // MPI collective to invoke.
  std::vector<ParticleList *> cell_ptrs;
  cell_ptrs.reserve(static_cast<std::size_t>(m_comm.size()));
  for (int n = 0; n < m_comm.size(); ++n) {
    cell_ptrs.push_back(&cells.at(static_cast<std::size_t>(n)).particles());
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
  for (auto &p : local().particles()) {
    m_box.fold_position(p.pos(), p.image_box());

    p.pos_at_last_verlet_update() = p.pos();
  }

  /* Local updates are a NoOp for this decomposition. */
  if (not global_flag) {
    return;
  }

  /* Sort displaced particles by the node they belong to. */
  std::vector<std::vector<Particle>> send_buf(m_comm.size());
  for (auto it = local().particles().begin();
       it != local().particles().end();) {
    auto const target_node = id_to_rank(it->id());
    if (target_node != m_comm.rank()) {
      diff.emplace_back(RemovedParticle{it->id()});
      send_buf.at(target_node).emplace_back(std::move(*it));
      it = local().particles().erase(it);
    } else {
      ++it;
    }
  }

  /* Exchange particles */
  std::vector<std::vector<Particle>> recv_buf(m_comm.size());
  boost::mpi::all_to_all(m_comm, send_buf, recv_buf);

  diff.emplace_back(ModifiedList{local().particles()});

  /* Add new particles belonging to this node */
  for (auto &parts : recv_buf) {
    for (auto &p : parts) {
      local().particles().insert(std::move(p));
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
