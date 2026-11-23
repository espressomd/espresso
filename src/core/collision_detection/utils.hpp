/*
 * Copyright (C) 2011-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_COLLISION_DETECTION

#include "CollisionPair.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "bonds.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "virtual_sites.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/utility.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace CollisionDetection {

inline auto &get_part(CellStructure &cell_structure, int id) {
  auto const p = cell_structure.get_local_particle(id);

  if (not p) {
    throw std::runtime_error("Could not handle collision because particle " +
                             std::to_string(id) + " was not found.");
  }

  return *p;
}

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
inline void place_vs_and_relate_to_particle(
    CellStructure &cell_structure, BoxGeometry const &box_geo,
    int const part_type_vs, double const min_global_cut,
    int const current_vs_pid, Utils::Vector3d const &pos, int const relate_to) {
  Particle new_part;
  new_part.id() = current_vs_pid;
  new_part.pos() = pos;
  auto p_vs = cell_structure.add_particle(std::move(new_part));
  vs_relate_to(*p_vs, get_part(cell_structure, relate_to), box_geo,
               min_global_cut);
  p_vs->type() = part_type_vs;
}
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE

inline auto gather_collision_queue(std::vector<CollisionPair> const &local) {
  auto global = local;
  if (::comm_cart.size() > 1) {
    Utils::Mpi::gather_buffer(global, ::comm_cart);
    boost::mpi::broadcast(::comm_cart, global, 0);
  }
  return global;
}

/**
 * @brief Bind the centers of the colliding particles.
 *
 * The local queues are gathered into the global collision queue, because
 * the global bond storage must be modified identically on every rank, while
 * a collision across a rank boundary may be queued on only one of the ranks,
 * or on both. Each pair is ordered as (smaller id, larger id) and duplicates
 * are removed, so that exactly one bond is created per colliding pair,
 * owned by the particle with the smaller id.
 *
 * The bond is also registered in the flat bond lists of the rank where its
 * owner is a non-ghost particle, which avoids rebuilding these lists. This
 * matches the rank that evaluates the bond after the next rebuild.
 *
 * @param[in] collision_queue  Local queue.
 * @return The global collision queue, identical on all ranks, without
 *         duplicates, and with each pair ordered as (smaller id, larger id).
 */
inline auto add_bind_centers(std::vector<CollisionPair> const &collision_queue,
                             System::System &system, int bond_id) {
  auto &cell_structure = *system.cell_structure;
  // gather and put every pair in a fixed order: (smaller id, larger id)
  auto global_queue = gather_collision_queue(collision_queue);
  for (auto &c : global_queue) {
    if (c.first > c.second) {
      std::swap(c.first, c.second);
    }
  }
  std::ranges::sort(global_queue);
  auto const [first, last] = std::ranges::unique(global_queue);
  global_queue.erase(first, last);

  for (auto const &c : global_queue) {
    // bond list on all ranks: the same single bond, owned by c.first
    ::add_bond(system, bond_id, {c.first, c.second});
    // local list for this step: only the rank where c.first is a real
    // (non-ghost) particle adds it
    auto const *p = cell_structure.get_local_particle(c.first);
    if (p and not p->is_ghost()) {
      cell_structure.add_new_bond(bond_id, {c.first, c.second});
    }
  }
  return global_queue;
}

} // namespace CollisionDetection

#endif // ESPRESSO_COLLISION_DETECTION
