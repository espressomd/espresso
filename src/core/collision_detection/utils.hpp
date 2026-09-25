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

#include <cassert>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace CollisionDetection {

// get_local_particle returns a by-value view, so this must return by value too
// -- returning a reference would dangle into the local optional. A Particle is
// a 16-byte handle; by-value is cheap.
inline Particle get_part(CellStructure &cell_structure, int id) {
  auto p = cell_structure.get_local_particle(id);

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
  // Build the new virtual site into a staging-store row via a view, then hand
  // it to add_particle (which stages the underlying row).
  auto new_part = cell_structure.make_new_particle_view();
  new_part.id() = current_vs_pid;
  new_part.pos() = pos;
  // add_particle returns a by-value view (nullopt if not local); the caller
  // always creates the VS on a rank that owns it, so it resolves here.
  auto p_vs = cell_structure.add_particle(std::move(new_part));
  assert(p_vs.has_value());
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

inline void add_bind_centers(std::vector<CollisionPair> &collision_queue,
                             System::System &system, int bond_id) {
  for (auto &c : collision_queue) {
    // Ensure that the bond is associated with the non-ghost particle
    if (system.cell_structure->get_local_particle(c.first)->is_ghost()) {
      std::swap(c.first, c.second);
    }

    // Because MPI rank 1's queue containing (@c p1_on_rank_1, @c p2_on_rank_2)
    // doesn't guarantee that the same pair (with or without swapped order) is
    // also queued on the MPI rank 2.
    // Once we change bond storage, some syncing has to be done.
    assert(use_one_sided_bond_storage);
    ::add_bond(system, bond_id, {c.first, c.second});
    system.cell_structure->add_new_bond(bond_id, {c.first, c.second});
  }
}

} // namespace CollisionDetection

#endif // ESPRESSO_COLLISION_DETECTION
