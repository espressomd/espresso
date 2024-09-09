/*
 * Copyright (C) 2011-2024 The ESPResSo project
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

#include <config/config.hpp>

#ifdef COLLISION_DETECTION
#ifdef VIRTUAL_SITES_RELATIVE

#include "BindAtPointOfCollision.hpp"
#include "CollisionPair.hpp"
#include "utils.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"
#include "virtual_sites.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/all_compare.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/operations.hpp>

#include <cassert>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CollisionDetection {

void BindAtPointOfCollision::initialize(System::System &system) {
  // Validate distance
  if (distance <= 0.) {
    throw std::domain_error("Parameter 'distance' must be > 0");
  }
  // Cache square of cutoff
  distance_sq = Utils::sqr(distance);

  // Check vs placement parameter
  if (vs_placement < 0. or vs_placement > 1.) {
    throw std::domain_error("Parameter 'vs_placement' must be between 0 and 1");
  }

  // Check if bond exists
  assert(system.bonded_ias->contains(bond_vs));
  // The bond between the virtual sites can be pair or triple
  auto const n_partners = number_of_partners(*system.bonded_ias->at(bond_vs));
  if (n_partners != 1 and n_partners != 2) {
    throw std::runtime_error("The bond type to be used for binding virtual "
                             "sites needs to be a pair bond");
  }

  // Create particle types
  if (part_type_vs < 0) {
    throw std::domain_error("Collision detection particle type for virtual "
                            "sites needs to be >=0");
  }
  system.nonbonded_ias->make_particle_type_exist(part_type_vs);
}

void BindAtPointOfCollision::handle_collisions(
    System::System &system, std::vector<CollisionPair> &local_collision_queue) {
  auto &cell_structure = *system.cell_structure;
  auto const min_global_cut = system.get_min_global_cut();
  auto const &box_geo = *system.box_geo;

  add_bind_centers(local_collision_queue, cell_structure, bond_centers);

  // Gather the global collision queue, because only one node has a collision
  // across node boundaries in its queue.
  // The other node might still have to change particle properties on its
  // non-ghost particle
  auto global_collision_queue = gather_collision_queue(local_collision_queue);

  // Synchornize max_seen_part
  auto const global_max_seen_particle = boost::mpi::all_reduce(
      ::comm_cart, cell_structure.get_max_local_particle_id(),
      boost::mpi::maximum<int>());

  int current_vs_pid = global_max_seen_particle + 1;

  // Iterate over global collision queue
  for (auto &c : global_collision_queue) {

    // Get particle pointers
    Particle *p1 = cell_structure.get_local_particle(c.first);
    Particle *p2 = cell_structure.get_local_particle(c.second);

    // Only nodes take part in particle creation and binding
    // that see both particles

    // If we cannot access both particles, both are ghosts,
    // or one is ghost and one is not accessible
    // we only increase the counter for the ext id to use based on the
    // number of particles created by other nodes
    if (((!p1 or p1->is_ghost()) and (!p2 or p2->is_ghost())) or !p1 or !p2) {
      // Increase local counters
      current_vs_pid += 2;
      continue;
    }
    // We consider the pair because one particle is local to the node and
    // the other is local or ghost.

    // Enable rotation on the particles to which vs will be attached
    p1->set_can_rotate_all_axes();
    p2->set_can_rotate_all_axes();

    // Positions of the virtual sites
    auto const vec21 = box_geo.get_mi_vector(p1->pos(), p2->pos());
    auto const pos1 = p1->pos() - vec21 * vs_placement;
    auto const pos2 = p1->pos() - vec21 * (1. - vs_placement);

    auto handle_particle = [&](Particle *p, Utils::Vector3d const &pos) {
      if (not p->is_ghost()) {
        place_vs_and_relate_to_particle(cell_structure, box_geo, part_type_vs,
                                        min_global_cut, current_vs_pid, pos,
                                        p->id());
        // Particle storage locations may have changed due to
        // added particle
        p1 = cell_structure.get_local_particle(c.first);
        p2 = cell_structure.get_local_particle(c.second);
      }
    };

    // place virtual sites on the node where the base particle is not a ghost
    handle_particle(p1, pos1);
    // Increment counter
    current_vs_pid++;

    handle_particle(p2, pos2);
    // Increment counter
    current_vs_pid++;

    // Create bonds
    auto const n_partners = number_of_partners(*system.bonded_ias->at(bond_vs));
    if (n_partners == 1) {
      // Create bond between the virtual particles
      const int bondG[] = {current_vs_pid - 2};
      // Only add bond if vs was created on this node
      if (auto p = cell_structure.get_local_particle(current_vs_pid - 1))
        p->bonds().insert({bond_vs, bondG});
    }
    if (n_partners == 2) {
      // Create 1st bond between the virtual particles
      const int bondG[] = {c.first, c.second};
      // Only add bond if vs was created on this node
      if (auto p = cell_structure.get_local_particle(current_vs_pid - 1))
        p->bonds().insert({bond_vs, bondG});
      if (auto p = cell_structure.get_local_particle(current_vs_pid - 2))
        p->bonds().insert({bond_vs, bondG});
    }
  } // Loop over all collisions in the queue

#ifdef ADDITIONAL_CHECKS
  assert(Utils::Mpi::all_compare(::comm_cart, current_vs_pid) &&
         "Nodes disagree about current_vs_pid");
#endif

  // If any node had a collision, all nodes need to resort
  if (not global_collision_queue.empty()) {
    cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
    cell_structure.update_ghosts_and_resort_particle(
        Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS);
    system.update_used_propagations();
  }
}

} // namespace CollisionDetection

#endif // VIRTUAL_SITES_RELATIVE
#endif // COLLISION_DETECTION
