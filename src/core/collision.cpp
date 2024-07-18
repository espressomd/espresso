/*
 * Copyright (C) 2011-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef COLLISION_DETECTION

#include "collision.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "system/System.hpp"
#include "virtual_sites.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/all_compare.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/serialization.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace boost {
namespace serialization {
template <typename Archive>
void serialize(Archive &ar, CollisionPair &c, const unsigned int) {
  ar & c.pp1;
  ar & c.pp2;
}
} // namespace serialization
} // namespace boost

namespace {
Particle &get_part(CellStructure &cell_structure, int id) {
  auto const p = cell_structure.get_local_particle(id);

  if (not p) {
    throw std::runtime_error("Could not handle collision because particle " +
                             std::to_string(id) + " was not found.");
  }

  return *p;
}
} // namespace

void CollisionDetection::initialize() {
  // If mode is OFF, no further checks
  if (mode == CollisionModeType::OFF) {
    return;
  }
  // Validate distance
  if (mode != CollisionModeType::OFF) {
    if (distance <= 0.) {
      throw std::domain_error("Parameter 'distance' must be > 0");
    }

    // Cache square of cutoff
    distance_sq = Utils::sqr(distance);
  }

#ifndef VIRTUAL_SITES_RELATIVE
  // The collision modes involving virtual sites also require the creation of a
  // bond between the colliding
  // If we don't have virtual sites, virtual site binding isn't possible.
  if ((mode == CollisionModeType::BIND_VS) or
      (mode == CollisionModeType::GLUE_TO_SURF)) {
    throw std::runtime_error("collision modes based on virtual sites require "
                             "the VIRTUAL_SITES_RELATIVE feature");
  }
#endif

#ifdef VIRTUAL_SITES
  // Check vs placement parameter
  if (mode == CollisionModeType::BIND_VS) {
    if (vs_placement < 0. or vs_placement > 1.) {
      throw std::domain_error(
          "Parameter 'vs_placement' must be between 0 and 1");
    }
  }
#endif

  auto &system = get_system();
  auto &bonded_ias = *system.bonded_ias;
  auto &nonbonded_ias = *system.nonbonded_ias;

  // Check if bond exists
  assert(mode != CollisionModeType::BIND_CENTERS or
         bonded_ias.contains(bond_centers));
  assert(mode != CollisionModeType::BIND_VS or bonded_ias.contains(bond_vs));

  // If the bond type to bind particle centers is not a pair bond...
  // Check that the bonds have the right number of partners
  if ((mode == CollisionModeType::BIND_CENTERS) and
      (number_of_partners(*bonded_ias.at(bond_centers)) != 1)) {
    throw std::runtime_error("The bond type to be used for binding particle "
                             "centers needs to be a pair bond");
  }

  // The bond between the virtual sites can be pair or triple
  if ((mode == CollisionModeType::BIND_VS) and
      (number_of_partners(*bonded_ias.at(bond_vs)) != 1 and
       number_of_partners(*bonded_ias.at(bond_vs)) != 2)) {
    throw std::runtime_error("The bond type to be used for binding virtual "
                             "sites needs to be a pair bond");
  }

  // Create particle types

  if (mode == CollisionModeType::BIND_VS) {
    if (vs_particle_type < 0) {
      throw std::domain_error("Collision detection particle type for virtual "
                              "sites needs to be >=0");
    }
    nonbonded_ias.make_particle_type_exist(vs_particle_type);
  }

  if (mode == CollisionModeType::GLUE_TO_SURF) {
    if (vs_particle_type < 0) {
      throw std::domain_error("Collision detection particle type for virtual "
                              "sites needs to be >=0");
    }
    nonbonded_ias.make_particle_type_exist(vs_particle_type);

    if (part_type_to_be_glued < 0) {
      throw std::domain_error("Collision detection particle type to be glued "
                              "needs to be >=0");
    }
    nonbonded_ias.make_particle_type_exist(part_type_to_be_glued);

    if (part_type_to_attach_vs_to < 0) {
      throw std::domain_error("Collision detection particle type to attach "
                              "the virtual site to needs to be >=0");
    }
    nonbonded_ias.make_particle_type_exist(part_type_to_attach_vs_to);

    if (part_type_after_glueing < 0) {
      throw std::domain_error("Collision detection particle type after gluing "
                              "needs to be >=0");
    }
    nonbonded_ias.make_particle_type_exist(part_type_after_glueing);
  }

  system.on_short_range_ia_change();
}

/** @brief Calculate position of vs for GLUE_TO_SURFACE mode.
 *  Returns id of particle to bind vs to.
 */
static auto const &glue_to_surface_calc_vs_pos(
    Particle const &p1, Particle const &p2, BoxGeometry const &box_geo,
    CollisionDetection const &collision_params, Utils::Vector3d &pos) {
  double ratio;
  auto const vec21 = box_geo.get_mi_vector(p1.pos(), p2.pos());
  auto const dist = vec21.norm();

  // Find out, which is the particle to be glued.
  if ((p1.type() == collision_params.part_type_to_be_glued) and
      (p2.type() == collision_params.part_type_to_attach_vs_to)) {
    ratio = 1. - collision_params.dist_glued_part_to_vs / dist;
  } else if ((p2.type() == collision_params.part_type_to_be_glued) and
             (p1.type() == collision_params.part_type_to_attach_vs_to)) {
    ratio = collision_params.dist_glued_part_to_vs / dist;
  } else {
    throw std::runtime_error("This should never be thrown. Bug.");
  }
  pos = p2.pos() + vec21 * ratio;
  if (p1.type() == collision_params.part_type_to_attach_vs_to)
    return p1;

  return p2;
}

static void bind_at_point_of_collision_calc_vs_pos(
    Particle const &p1, Particle const &p2, BoxGeometry const &box_geo,
    CollisionDetection const &collision_params, Utils::Vector3d &pos1,
    Utils::Vector3d &pos2) {
  auto const vec21 = box_geo.get_mi_vector(p1.pos(), p2.pos());
  pos1 = p1.pos() - vec21 * collision_params.vs_placement;
  pos2 = p1.pos() - vec21 * (1. - collision_params.vs_placement);
}

#ifdef VIRTUAL_SITES_RELATIVE
static void place_vs_and_relate_to_particle(
    CellStructure &cell_structure, BoxGeometry const &box_geo,
    CollisionDetection const &collision_params, double const min_global_cut,
    int const current_vs_pid, Utils::Vector3d const &pos, int const relate_to) {
  Particle new_part;
  new_part.id() = current_vs_pid;
  new_part.pos() = pos;
  auto p_vs = cell_structure.add_particle(std::move(new_part));
  vs_relate_to(*p_vs, get_part(cell_structure, relate_to), box_geo,
               min_global_cut);
  p_vs->type() = collision_params.vs_particle_type;
}

static void bind_at_poc_create_bond_between_vs(
    CellStructure &cell_structure, BondedInteractionsMap const &bonded_ias,
    CollisionDetection const &collision_params, int const current_vs_pid,
    CollisionPair const &c) {
  switch (number_of_partners(*bonded_ias.at(collision_params.bond_vs))) {
  case 1: {
    // Create bond between the virtual particles
    const int bondG[] = {current_vs_pid - 2};
    // Only add bond if vs was created on this node
    if (auto p = cell_structure.get_local_particle(current_vs_pid - 1))
      p->bonds().insert({collision_params.bond_vs, bondG});
    break;
  }
  case 2: {
    // Create 1st bond between the virtual particles
    const int bondG[] = {c.pp1, c.pp2};
    // Only add bond if vs was created on this node
    if (auto p = cell_structure.get_local_particle(current_vs_pid - 1))
      p->bonds().insert({collision_params.bond_vs, bondG});
    if (auto p = cell_structure.get_local_particle(current_vs_pid - 2))
      p->bonds().insert({collision_params.bond_vs, bondG});
    break;
  }
  }
}

static void glue_to_surface_bind_part_to_vs(
    Particle const *const p1, Particle const *const p2,
    int const vs_pid_plus_one, CollisionDetection const &collision_params,
    CellStructure &cell_structure) {
  // Create bond between the virtual particles
  const int bondG[] = {vs_pid_plus_one - 1};

  if (p1->type() == collision_params.part_type_after_glueing) {
    get_part(cell_structure, p1->id())
        .bonds()
        .insert({collision_params.bond_vs, bondG});
  } else {
    get_part(cell_structure, p2->id())
        .bonds()
        .insert({collision_params.bond_vs, bondG});
  }
}

#endif // VIRTUAL_SITES_RELATIVE

static auto gather_collision_queue(std::vector<CollisionPair> const &local) {
  auto global = local;
  if (comm_cart.size() > 1) {
    Utils::Mpi::gather_buffer(global, ::comm_cart);
    boost::mpi::broadcast(::comm_cart, global, 0);
  }
  return global;
}

void CollisionDetection::handle_collisions(CellStructure &cell_structure) {
  // Note that the glue to surface mode adds bonds between the centers
  // but does so later in the process. This is needed to guarantee that
  // a particle can only be glued once, even if queued twice in a single
  // time step
  auto const bind_centers = mode != CollisionModeType::OFF and
                            mode != CollisionModeType::GLUE_TO_SURF;
  if (bind_centers) {
    for (auto &c : local_collision_queue) {
      // put the bond to the non-ghost particle; at least one partner always is
      if (cell_structure.get_local_particle(c.pp1)->is_ghost()) {
        std::swap(c.pp1, c.pp2);
      }

      const int bondG[] = {c.pp2};

      get_part(cell_structure, c.pp1).bonds().insert({bond_centers, bondG});
    }
  }

#ifdef VIRTUAL_SITES_RELATIVE
  auto &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const min_global_cut = system.get_min_global_cut();
  if ((mode == CollisionModeType::BIND_VS) ||
      (mode == CollisionModeType::GLUE_TO_SURF)) {
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
    for (auto &c : global_collision_queue) {

      // Get particle pointers
      Particle *p1 = cell_structure.get_local_particle(c.pp1);
      Particle *p2 = cell_structure.get_local_particle(c.pp2);

      // Only nodes take part in particle creation and binding
      // that see both particles

      // If we cannot access both particles, both are ghosts,
      // or one is ghost and one is not accessible
      // we only increase the counter for the ext id to use based on the
      // number of particles created by other nodes
      if (((!p1 or p1->is_ghost()) and (!p2 or p2->is_ghost())) or !p1 or !p2) {
        // Increase local counters
        if (mode == CollisionModeType::BIND_VS) {
          current_vs_pid++;
        }
        // For glue to surface, we have only one vs
        current_vs_pid++;
        if (mode == CollisionModeType::GLUE_TO_SURF) {
          if (p1 and p1->type() == part_type_to_be_glued) {
            p1->type() = part_type_after_glueing;
          }
          if (p2 and p2->type() == part_type_to_be_glued) {
            p2->type() = part_type_after_glueing;
          }
        } // mode glue to surface

      } else { // We consider the pair because one particle
               // is local to the node and the other is local or ghost
        // If we are in the two vs mode
        // Virtual site related to first particle in the collision
        if (mode == CollisionModeType::BIND_VS) {
          Utils::Vector3d pos1, pos2;

          // Enable rotation on the particles to which vs will be attached
          p1->set_can_rotate_all_axes();
          p2->set_can_rotate_all_axes();

          // Positions of the virtual sites
          bind_at_point_of_collision_calc_vs_pos(*p1, *p2, box_geo, *this, pos1,
                                                 pos2);

          auto handle_particle = [&](Particle *p, Utils::Vector3d const &pos) {
            if (not p->is_ghost()) {
              place_vs_and_relate_to_particle(cell_structure, box_geo, *this,
                                              min_global_cut, current_vs_pid,
                                              pos, p->id());
              // Particle storage locations may have changed due to
              // added particle
              p1 = cell_structure.get_local_particle(c.pp1);
              p2 = cell_structure.get_local_particle(c.pp2);
            }
          };

          // place virtual sites on the node where the base particle is not a
          // ghost
          handle_particle(p1, pos1);
          // Increment counter
          current_vs_pid++;

          handle_particle(p2, pos2);
          // Increment counter
          current_vs_pid++;
          // Create bonds between the vs.

          bind_at_poc_create_bond_between_vs(cell_structure, *system.bonded_ias,
                                             *this, current_vs_pid, c);
        } // mode VS

        if (mode == CollisionModeType::GLUE_TO_SURF) {
          // If particles are made inert by a type change on collision:
          // We skip the pair if one of the particles has already reacted
          // but we still increase the particle counters, as other nodes
          // can not always know whether or not a vs is placed
          if (part_type_after_glueing != part_type_to_be_glued) {
            if ((p1->type() == part_type_after_glueing) ||
                (p2->type() == part_type_after_glueing)) {
              current_vs_pid++;
              continue;
            }
          }

          Utils::Vector3d pos;
          Particle const &attach_vs_to =
              glue_to_surface_calc_vs_pos(*p1, *p2, box_geo, *this, pos);

          // Add a bond between the centers of the colliding particles
          // The bond is placed on the node that has p1
          if (!p1->is_ghost()) {
            const int bondG[] = {c.pp2};
            get_part(cell_structure, c.pp1)
                .bonds()
                .insert({bond_centers, bondG});
          }

          // Change type of particle being attached, to make it inert
          if (p1->type() == part_type_to_be_glued) {
            p1->type() = part_type_after_glueing;
          }
          if (p2->type() == part_type_to_be_glued) {
            p2->type() = part_type_after_glueing;
          }

          // Vs placement happens on the node that has p1
          if (!attach_vs_to.is_ghost()) {
            place_vs_and_relate_to_particle(cell_structure, box_geo, *this,
                                            min_global_cut, current_vs_pid, pos,
                                            attach_vs_to.id());
            // Particle storage locations may have changed due to
            // added particle
            p1 = cell_structure.get_local_particle(c.pp1);
            p2 = cell_structure.get_local_particle(c.pp2);
            current_vs_pid++;
          } else { // Just update the books
            current_vs_pid++;
          }
          glue_to_surface_bind_part_to_vs(p1, p2, current_vs_pid, *this,
                                          cell_structure);
        }
      } // we considered the pair
    } // Loop over all collisions in the queue
#ifdef ADDITIONAL_CHECKS
    if (!Utils::Mpi::all_compare(comm_cart, current_vs_pid)) {
      throw std::runtime_error("Nodes disagree about current_vs_pid");
    }
#endif

    // If any node had a collision, all nodes need to resort
    if (not global_collision_queue.empty()) {
      cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
      cell_structure.update_ghosts_and_resort_particle(
          Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS);
    }
    system.update_used_propagations();
  } // are we in one of the vs_based methods
#endif // defined VIRTUAL_SITES_RELATIVE

  clear_queue();
}

#endif // COLLISION_DETECTION
