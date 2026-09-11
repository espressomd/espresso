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

#include <config/config.hpp>

#ifdef ESPRESSO_COLLISION_DETECTION
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE

#include "CollisionPair.hpp"
#include "GlueToSurface.hpp"
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
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CollisionDetection {

void GlueToSurface::initialize(System::System &system) {
  // Validate distance
  if (distance <= 0.) {
    throw std::domain_error("Parameter 'distance' must be > 0");
  }
  // Cache square of cutoff
  distance_sq = Utils::sqr(distance);

  assert(system.bonded_ias->contains(bond_centers));
  if (number_of_partners(*system.bonded_ias->at(bond_centers)) != 1) {
    throw std::runtime_error("The bond type to be used for binding particle "
                             "centers needs to be a pair bond");
  }
  assert(system.bonded_ias->contains(bond_vs));
  if (number_of_partners(*system.bonded_ias->at(bond_vs)) != 1) {
    throw std::runtime_error("The bond type to be used for binding virtual "
                             "sites needs to be a pair bond");
  }

  if (part_type_vs < 0) {
    throw std::domain_error("Collision detection particle type for virtual "
                            "sites needs to be >=0");
  }
  system.nonbonded_ias->make_particle_type_exist(part_type_vs);

  if (part_type_to_be_glued < 0) {
    throw std::domain_error("Collision detection particle type to be glued "
                            "needs to be >=0");
  }
  system.nonbonded_ias->make_particle_type_exist(part_type_to_be_glued);

  if (part_type_to_attach_vs_to < 0) {
    throw std::domain_error("Collision detection particle type to attach "
                            "the virtual site to needs to be >=0");
  }
  system.nonbonded_ias->make_particle_type_exist(part_type_to_attach_vs_to);

  if (part_type_after_glueing < 0) {
    throw std::domain_error("Collision detection particle type after gluing "
                            "needs to be >=0");
  }
  system.nonbonded_ias->make_particle_type_exist(part_type_after_glueing);
}

void GlueToSurface::handle_collisions(
    System::System &system, std::vector<CollisionPair> &local_collision_queue) {
  auto &cell_structure = *system.cell_structure;
  auto const min_global_cut = system.get_min_global_cut();
  auto const &box_geo = *system.box_geo;
  // Note that the glue to surface mode adds bonds between the centers
  // but does so later in the process. This is needed to guarantee that
  // a particle can only be glued once, even if queued twice in a single
  // time step

  // Gather the global collision queue, because only one node has a collision
  // across node boundaries in its queue.
  // The other node might still have to change particle properties on its
  // non-ghost particle
  auto global_collision_queue = gather_collision_queue(local_collision_queue);

  // Sync max_seen_part
  auto const global_max_seen_particle = boost::mpi::all_reduce(
      ::comm_cart, cell_structure.get_max_local_particle_id(),
      boost::mpi::maximum<int>());

  int current_vs_pid = global_max_seen_particle + 1;

  // Iterate over global collision queue
  for (auto const &[pid1, pid2] : global_collision_queue) {

    // get_local_particle returns by-value views. p1/p2 are held across
    // place_vs_and_relate_to_particle, which calls add_particle -> a store
    // rebuild that RENUMBERS every row and bumps the generation, so they must
    // be re-resolved by id afterwards.
    std::optional<Particle> p1 = cell_structure.get_local_particle(pid1);
    std::optional<Particle> p2 = cell_structure.get_local_particle(pid2);
    auto const resolved_generation =
        cell_structure.particle_store().generation();

    // Only nodes take part in particle creation and binding
    // that see both particles

    // If we cannot access both particles, both are ghosts,
    // or one is ghost and one is not accessible
    // we only increase the counter for the ext id to use based on the
    // number of particles created by other nodes
    if (((!p1 or p1->is_ghost()) and (!p2 or p2->is_ghost())) or !p1 or !p2) {
      // For glue to surface, we have only one vs
      current_vs_pid++;

      if (p1 and p1->type() == part_type_to_be_glued) {
        p1->type() = part_type_after_glueing;
      }
      if (p2 and p2->type() == part_type_to_be_glued) {
        p2->type() = part_type_after_glueing;
      }
      continue;
    }
    // If particles are made inert by a type change on collision:
    // We skip the pair if one of the particles has already reacted
    // but we still increase the particle counters, as other nodes
    // can not always know whether or not a vs is placed
    if (part_type_after_glueing != part_type_to_be_glued) {
      if ((p1->type() == part_type_after_glueing) or
          (p2->type() == part_type_after_glueing)) {
        current_vs_pid++;
        continue;
      }
    }

    double ratio = -1.;
    auto const vec21 = box_geo.get_mi_vector(Utils::Vector3d(p1->pos()),
                                             Utils::Vector3d(p2->pos()));
    auto const dist = vec21.norm();

    // Find out, which is the particle to be glued.
    if ((p1->type() == part_type_to_be_glued) and
        (p2->type() == part_type_to_attach_vs_to)) {
      ratio = 1. - dist_glued_part_to_vs / dist;
    } else if ((p2->type() == part_type_to_be_glued) and
               (p1->type() == part_type_to_attach_vs_to)) {
      ratio = dist_glued_part_to_vs / dist;
    }
    assert(ratio != -1.);
    auto const pos = Utils::Vector3d(p2->pos()) + vec21 * ratio;
    // Capture the base-particle id and its ghost flag BY VALUE before
    // add_particle. add_particle triggers a store rebuild that renumbers rows,
    // so reading id()/is_ghost() from a view held across that call would alias
    // whatever now occupies the old row. Snapshotting the scalar values
    // decouples the decision from row identity.
    auto const attach_vs_to =
        (p1->type() == part_type_to_attach_vs_to) ? *p1 : *p2;
    auto const attach_vs_to_id = attach_vs_to.id();
    auto const attach_vs_to_is_ghost = attach_vs_to.is_ghost();

    // Add a bond between the centers of the colliding particles
    // The bond is placed on the node that has p1
    if (!p1->is_ghost()) {
      const int bondG[] = {pid2};
      get_part(cell_structure, pid1).bonds().insert({bond_centers, bondG});
    }

    // Change type of particle being attached, to make it inert
    if (p1->type() == part_type_to_be_glued) {
      p1->type() = part_type_after_glueing;
    }
    if (p2->type() == part_type_to_be_glued) {
      p2->type() = part_type_after_glueing;
    }

    if (attach_vs_to_is_ghost) {
      current_vs_pid++;
    } else {
      // Debug canary: guard the held p1/p2 views before the add invalidates
      // them (see ParticleStoreGuard). attach_vs_to_id was snapshotted above.
      ParticleStoreGuard::assert_generation(cell_structure.particle_store(),
                                            resolved_generation);
      // VS placement happens on the node that has p1
      place_vs_and_relate_to_particle(cell_structure, box_geo, part_type_vs,
                                      min_global_cut, current_vs_pid, pos,
                                      attach_vs_to_id);
      // Particle storage locations changed due to the added particle:
      // re-resolve both base-particle views by id. The re-resolved views are
      // used immediately below with no further topology change in between, so
      // no fresh generation stamp is needed (the per-iteration stamp above is
      // const and guards the views held ACROSS the add).
      p1 = cell_structure.get_local_particle(pid1);
      p2 = cell_structure.get_local_particle(pid2);
      current_vs_pid++;
    }
    // Create bond between the virtual particles
    auto const glued_pid =
        (p1->type() == part_type_after_glueing) ? p1->id() : p2->id();
    int const bondG[] = {current_vs_pid - 1};
    get_part(cell_structure, glued_pid).bonds().insert({bond_vs, bondG});
  } // Loop over all collisions in the queue

#ifdef ESPRESSO_ADDITIONAL_CHECKS
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

#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#endif // ESPRESSO_COLLISION_DETECTION
