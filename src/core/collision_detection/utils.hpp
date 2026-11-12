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
 * The local queue is gathered first, so that every rank creates the bond
 * entries of all collisions: the rank that detected a collision across a
 * rank boundary only holds a ghost copy of the bond partner, and an entry
 * written there never reaches the rank owning the partner.
 *
 * @param[in,out] collision_queue  Local queue. Pairs are reordered such
 *                                 that the first particle is not a ghost.
 * @return The global collision queue, identical on all ranks.
 */
inline auto add_bind_centers(std::vector<CollisionPair> &collision_queue,
                             System::System &system, int bond_id) {
  auto &cell_structure = *system.cell_structure;
  for (auto &c : collision_queue) {
    // Ensure that the bond is associated with the non-ghost particle
    if (cell_structure.get_local_particle(c.first)->is_ghost()) {
      std::swap(c.first, c.second);
    }
    // The bond is evaluated on the rank that detected the collision.
    cell_structure.add_new_bond(bond_id, {c.first, c.second});
  }

  // Because MPI rank 1's queue containing (@c p1_on_rank_1, @c p2_on_rank_2)
  // doesn't guarantee that the same pair (with or without swapped order) is
  // also queued on the MPI rank 2, both ranks may independently create
  // the bond. This is orthogonal to bond storage and not addressed here.
  auto global_collision_queue = gather_collision_queue(collision_queue);
  for (auto const &c : global_collision_queue) {
    // writes the entries on whichever participants are locally known
    ::add_bond(system, bond_id, {c.first, c.second});
  }
  return global_collision_queue;
}

} // namespace CollisionDetection

#endif // ESPRESSO_COLLISION_DETECTION
