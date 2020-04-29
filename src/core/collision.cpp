/*
 * Copyright (C) 2011-2019 The ESPResSo project
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

#include "collision.hpp"

#ifdef COLLISION_DETECTION
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "virtual_sites/VirtualSitesRelative.hpp"

#include <utils/mpi/all_compare.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/algorithm/clamp.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/serialization.hpp>

#include <vector>

/// Data type holding the info about a single collision
typedef struct {
  int pp1; // 1st particle id
  int pp2; // 2nd particle id
} collision_struct;

namespace boost {
namespace serialization {
template <typename Archive>
void serialize(Archive &ar, collision_struct &c, const unsigned int) {
  ar &c.pp1;
  ar &c.pp2;
}
} // namespace serialization
} // namespace boost

/// During force calculation, colliding particles are recorded in the queue.
/// The queue is processed after force calculation, when it is safe to add
/// particles.
static std::vector<collision_struct> local_collision_queue;

/// Parameters for collision detection
Collision_parameters collision_params;

namespace {
Particle &get_part(int id) {
  auto const p = cell_structure.get_local_particle(id);

  if (not p) {
    throw std::runtime_error("Could not handle collision because particle " +
                             std::to_string(id) + " was not found.");
  }

  return *p;
}
} // namespace

/** @brief Return true if a bond between the centers of the colliding particles
 * needs to be placed. At this point, all modes need this */
inline bool bind_centers() {
  // Note that the glue to surface mode adds bonds between the centers
  // but does so later in the process. This is needed to guarantee that
  // a particle can only be glued once, even if queued twice in a single
  // time step
  return collision_params.mode != COLLISION_MODE_OFF &&
         collision_params.mode != COLLISION_MODE_GLUE_TO_SURF;
}

bool validate_collision_parameters() {
  // If mode is OFF, no further checks
  if (collision_params.mode == COLLISION_MODE_OFF) {
    return true;
  }
  // Validate distance
  if (collision_params.mode != COLLISION_MODE_OFF) {
    if (collision_params.distance <= 0.) {
      runtimeErrorMsg() << "collision_detection distance must be >0";
      return false;
    }

    // Cache square of cutoff
    collision_params.distance2 = Utils::sqr(collision_params.distance);

    if (collision_params.distance > min_global_cut) {
      runtimeErrorMsg() << "The minimum global cutoff (System.min_global_cut) "
                           "must be larger or equal the collision detection "
                           "distance.";
    }
  }

#ifndef VIRTUAL_SITES_RELATIVE
  // The collision modes involving virtual sites also require the creation of a
  // bond between the colliding
  // If we don't have virtual sites, virtual site binding isn't possible.
  if ((collision_params.mode & COLLISION_MODE_VS) ||
      (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    runtimeErrorMsg() << "Virtual sites based collision modes require "
                         "the VIRTUAL_SITES feature";
    return false;
  }
#endif

// Check vs placement parameter
#ifdef VIRTUAL_SITES
  if (collision_params.mode & COLLISION_MODE_VS) {
    if ((collision_params.vs_placement < 0) ||
        (collision_params.vs_placement > 1)) {
      runtimeErrorMsg() << "The collision detection vs_placement parameter "
                           "needs to be between 0 and 1.";
      return false;
    }
  }
#endif

  // Check if bonded ia exist
  if ((collision_params.mode & COLLISION_MODE_BOND) &&
      (collision_params.bond_centers >= bonded_ia_params.size())) {
    runtimeErrorMsg() << "The bond type to be used for binding particle "
                         "centers does not exist";
    return false;
  }

  if ((collision_params.mode & COLLISION_MODE_VS) &&
      (collision_params.bond_vs >= bonded_ia_params.size())) {
    runtimeErrorMsg()
        << "The bond type to be used for binding virtual sites does not exist";
    return false;
  }

  if ((collision_params.mode & COLLISION_MODE_BOND) &&
      collision_params.bond_centers == -1) {
    runtimeErrorMsg() << "The bond_centers parameter is unknown. Did you add "
                         "the interaction using system.bonded_inter.add?";
    return false;
  }

  // If the bond type to bind particle centers is not a pair bond...
  // Check that the bonds have the right number of partners
  if ((collision_params.mode & COLLISION_MODE_BOND) &&
      (bonded_ia_params[collision_params.bond_centers].num != 1)) {
    runtimeErrorMsg() << "The bond type to be used for binding particle "
                         "centers needs to be a pair bond";
    return false;
  }

  // The bond between the virtual sites can be pair or triple
  if ((collision_params.mode & COLLISION_MODE_VS) &&
      !(bonded_ia_params[collision_params.bond_vs].num == 1 ||
        bonded_ia_params[collision_params.bond_vs].num == 2)) {
    runtimeErrorMsg() << "The bond type to be used for binding virtual sites "
                         "needs to be a pair or three-particle bond";
    return false;
  }

  if (collision_params.mode & COLLISION_MODE_BIND_THREE_PARTICLES) {
    if (collision_params.bond_three_particles +
            collision_params.three_particle_angle_resolution >
        bonded_ia_params.size()) {
      runtimeErrorMsg()
          << "Insufficient bonds defined for three particle binding.";
      return false;
    }

    for (int i = collision_params.bond_three_particles;
         i < collision_params.bond_three_particles +
                 collision_params.three_particle_angle_resolution;
         i++) {
      if (bonded_ia_params[i].num != 2) {
        runtimeErrorMsg()
            << "The bonds for three particle binding need to be angle bonds.";
        return false;
      }
    }
  }

  // Create particle types

  if (collision_params.mode & COLLISION_MODE_VS) {
    if (collision_params.vs_particle_type < 0) {
      runtimeErrorMsg() << "Collision detection particle type for virtual "
                           "sites needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.vs_particle_type);
  }

  if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF) {
    if (collision_params.vs_particle_type < 0) {
      runtimeErrorMsg() << "Collision detection particle type for virtual "
                           "sites needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.vs_particle_type);

    if (collision_params.part_type_to_be_glued < 0) {
      runtimeErrorMsg()
          << "Collision detection particle type to be glued needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.part_type_to_be_glued);

    if (collision_params.part_type_to_attach_vs_to < 0) {
      runtimeErrorMsg() << "Collision detection particle type to attach the "
                           "virtual site to  needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.part_type_to_attach_vs_to);

    if (collision_params.part_type_after_glueing < 0) {
      runtimeErrorMsg()
          << "Collision detection particle type after gluing needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.part_type_after_glueing);
  }

  recalc_forces = true;
  rebuild_verletlist = true;

  return true;
}

void prepare_local_collision_queue() { local_collision_queue.clear(); }

void queue_collision(const int part1, const int part2) {
  local_collision_queue.push_back({part1, part2});
}

/** @brief Calculate position of vs for GLUE_TO_SURFACE mode
 *    Returns id of particle to bind vs to */
const Particle &glue_to_surface_calc_vs_pos(const Particle &p1,
                                            const Particle &p2,
                                            Utils::Vector3d &pos) {
  double c;
  auto const vec21 = get_mi_vector(p1.r.p, p2.r.p, box_geo);
  const double dist_betw_part = vec21.norm();

  // Find out, which is the particle to be glued.
  if ((p1.p.type == collision_params.part_type_to_be_glued) &&
      (p2.p.type == collision_params.part_type_to_attach_vs_to)) {
    c = 1 - collision_params.dist_glued_part_to_vs / dist_betw_part;
  } else if ((p2.p.type == collision_params.part_type_to_be_glued) &&
             (p1.p.type == collision_params.part_type_to_attach_vs_to)) {
    c = collision_params.dist_glued_part_to_vs / dist_betw_part;
  } else {
    throw std::runtime_error("This should never be thrown. Bug.");
  }
  for (int i = 0; i < 3; i++) {
    pos[i] = p2.r.p[i] + vec21[i] * c;
  }
  if (p1.p.type == collision_params.part_type_to_attach_vs_to)
    return p1;

  return p2;
}

void bind_at_point_of_collision_calc_vs_pos(const Particle *const p1,
                                            const Particle *const p2,
                                            Utils::Vector3d &pos1,
                                            Utils::Vector3d &pos2) {
  auto const vec21 = get_mi_vector(p1->r.p, p2->r.p, box_geo);
  pos1 = p1->r.p - vec21 * collision_params.vs_placement;
  pos2 = p1->r.p - vec21 * (1. - collision_params.vs_placement);
}

// Considers three particles for three_particle_binding and performs
// the binding if the criteria are met //
void coldet_do_three_particle_bond(Particle &p, Particle &p1, Particle &p2) {
  // If p1 and p2 are not closer or equal to the cutoff distance, skip
  // p1:
  if (get_mi_vector(p.r.p, p1.r.p, box_geo).norm() > collision_params.distance)
    return;
  // p2:
  if (get_mi_vector(p.r.p, p2.r.p, box_geo).norm() > collision_params.distance)
    return;

  // Check, if there already is a three-particle bond centered on p
  // with p1 and p2 as partners. If so, skip this triplet.
  // Note that the bond partners can appear in any order.

  /* Check if a bond is a bond placed by the collision detection */
  auto is_collision_bond = [](BondView const &bond) {
    return (bond.bond_id() >= collision_params.bond_three_particles) and
           (bond.bond_id() <=
            collision_params.bond_three_particles +
                collision_params.three_particle_angle_resolution);
  };
  /* Check if the bond is between the particles we are currently considering */
  auto has_same_partners = [id1 = p1.identity(),
                            id2 = p2.identity()](BondView const &bond) {
    auto const partner_ids = bond.partner_ids();

    return ((partner_ids[0] == id1) && (partner_ids[1] == id2)) ||
           ((partner_ids[0] == id2) && (partner_ids[1] == id1));
  };

  if (boost::algorithm::any_of(p.bonds(), [=](auto const &bond) {
        return is_collision_bond(bond) and has_same_partners(bond);
      })) {
    return;
  }

  // If we are still here, we need to create angular bond
  // First, find the angle between the particle p, p1 and p2

  /* vector from p to p1 */
  auto const vec1 = get_mi_vector(p.r.p, p1.r.p, box_geo).normalize();
  /* vector from p to p2 */
  auto const vec2 = get_mi_vector(p.r.p, p2.r.p, box_geo).normalize();

  auto const cosine =
      boost::algorithm::clamp(vec1 * vec2, -TINY_COS_VALUE, TINY_COS_VALUE);

  // Bond angle
  auto const phi = acos(cosine);

  // We find the bond id by dividing the range from 0 to pi in
  // three_particle_angle_resolution steps and by adding the id
  // of the bond for zero degrees.
  auto const bond_id = static_cast<int>(
      floor(phi / M_PI *
                (collision_params.three_particle_angle_resolution - 1) +
            0.5) +
      collision_params.bond_three_particles);

  // Create the bond

  // First, fill bond data structure
  const std::array<int, 2> bondT = {p1.p.identity, p2.p.identity};

  p.bonds().insert({bond_id, bondT});
}

#ifdef VIRTUAL_SITES_RELATIVE
void place_vs_and_relate_to_particle(const int current_vs_pid,
                                     const Utils::Vector3d &pos,
                                     int relate_to) {
  Particle new_part;
  new_part.p.identity = current_vs_pid;
  new_part.r.p = pos;
  auto p_vs = cell_structure.add_particle(std::move(new_part));

  local_vs_relate_to(*p_vs, get_part(relate_to));

  p_vs->p.is_virtual = true;
  p_vs->p.type = collision_params.vs_particle_type;
}

void bind_at_poc_create_bond_between_vs(const int current_vs_pid,
                                        const collision_struct &c) {
  switch (bonded_ia_params[collision_params.bond_vs].num) {
  case 1: {
    // Create bond between the virtual particles
    const int bondG[] = {current_vs_pid - 2};
    // Only add bond if vs was created on this node
    if (cell_structure.get_local_particle(current_vs_pid - 1))
      get_part(current_vs_pid - 1)
          .bonds()
          .insert({collision_params.bond_vs, bondG});
    break;
  }
  case 2: {
    // Create 1st bond between the virtual particles
    const int bondG[] = {c.pp1, c.pp2};
    // Only add bond if vs was created on this node
    if (cell_structure.get_local_particle(current_vs_pid - 1))
      get_part(current_vs_pid - 1)
          .bonds()
          .insert({collision_params.bond_vs, bondG});
    if (cell_structure.get_local_particle(current_vs_pid - 2))
      get_part(current_vs_pid - 2)
          .bonds()
          .insert({collision_params.bond_vs, bondG});
    break;
  }
  }
}

void glue_to_surface_bind_part_to_vs(const Particle *const p1,
                                     const Particle *const p2,
                                     const int vs_pid_plus_one,
                                     const collision_struct &c) {
  // Create bond between the virtual particles
  const int bondG[] = {vs_pid_plus_one - 1};

  if (p1->p.type == collision_params.part_type_after_glueing) {
    get_part(p1->p.identity).bonds().insert({collision_params.bond_vs, bondG});
  } else {
    get_part(p2->p.identity).bonds().insert({collision_params.bond_vs, bondG});
  }
}

#endif

std::vector<collision_struct> gather_global_collision_queue() {
  std::vector<collision_struct> res = local_collision_queue;
  Utils::Mpi::gather_buffer(res, comm_cart);
  boost::mpi::broadcast(comm_cart, res, 0);

  return res;
}

static void three_particle_binding_do_search(Cell *basecell, Particle &p1,
                                             Particle &p2) {
  auto handle_cell = [&p1, &p2](Cell *c) {
    for (auto &P : c->particles()) {
      // Skip collided particles themselves
      if ((P.p.identity == p1.p.identity) || (P.p.identity == p2.p.identity)) {
        continue;
      }

      // We need all cyclical permutations, here (bond is placed on 1st
      // particle, order of bond partners does not matter, so we don't need
      // non-cyclic permutations).
      // coldet_do_three_particle_bond checks the bonding criterion and if
      // the involved particles are not already bonded before it binds them.
      if (!P.l.ghost) {
        coldet_do_three_particle_bond(P, p1, p2);
      }

      if (!p1.l.ghost) {
        coldet_do_three_particle_bond(p1, P, p2);
      }

      if (!p2.l.ghost) {
        coldet_do_three_particle_bond(p2, P, p1);
      }
    }
  };

  /* Search the base cell ... */
  handle_cell(basecell);

  /* ... and all the neighbors. */
  for (auto &n : basecell->neighbors().all()) {
    handle_cell(n);
  }
}

// Goes through the collision queue and for each pair in it
// looks for a third particle by using the domain decomposition
// cell system. If found, it performs three particle binding
void three_particle_binding_domain_decomposition(
    const std::vector<collision_struct> &gathered_queue) {

  for (auto &c : gathered_queue) {
    // If we have both particles, at least as ghosts, Get the corresponding cell
    // indices
    if (cell_structure.get_local_particle(c.pp1) &&
        cell_structure.get_local_particle(c.pp2)) {
      Particle &p1 = *cell_structure.get_local_particle(c.pp1);
      Particle &p2 = *cell_structure.get_local_particle(c.pp2);
      auto cell1 = find_current_cell(p1);
      auto cell2 = find_current_cell(p2);

      if (cell1)
        three_particle_binding_do_search(cell1, p1, p2);
      if (cell2 && cell1 != cell2)
        three_particle_binding_do_search(cell2, p1, p2);

    } // If local particles exist
  }   // Loop over total collisions
}

// Handle the collisions stored in the queue
void handle_collisions() {

  if (collision_params.exception_on_collision) {
    for (auto &c : local_collision_queue) {
      runtimeWarningMsg() << "Collision between particles "
                          << std::min(c.pp1, c.pp2) << " and "
                          << std::max(c.pp1, c.pp2);
    }
  }

  // Note that the glue to surface mode adds bonds between the centers
  // but does so later in the process. This is needed to guarantee that
  // a particle can only be glued once, even if queued twice in a single
  // time step
  if (bind_centers()) {
    for (auto &c : local_collision_queue) {
      // put the bond to the non-ghost particle; at least one partner always is
      if (cell_structure.get_local_particle(c.pp1)->l.ghost) {
        std::swap(c.pp1, c.pp2);
      }

      const int bondG[] = {c.pp2};

      get_part(c.pp1).bonds().insert({collision_params.bond_centers, bondG});
    }
  }

// Virtual sites based collision schemes
#ifdef VIRTUAL_SITES_RELATIVE
  if ((collision_params.mode & COLLISION_MODE_VS) ||
      (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    // Gather the global collision queue, because only one node has a collision
    // across node boundaries in its queue.
    // The other node might still have to change particle properties on its
    // non-ghost particle
    auto gathered_queue = gather_global_collision_queue();

    // Sync max_seen_part
    auto const global_max_seen_particle = boost::mpi::all_reduce(
        comm_cart, cell_structure.get_max_local_particle_id(),
        boost::mpi::maximum<int>());

    int current_vs_pid = global_max_seen_particle + 1;

    // Iterate over global collision queue
    for (auto &c : gathered_queue) {

      // Get particle pointers
      Particle *p1 = cell_structure.get_local_particle(c.pp1);
      Particle *p2 = cell_structure.get_local_particle(c.pp2);

      // Only nodes take part in particle creation and binding
      // that see both particles

      // If we cannot access both particles, both are ghosts,
      // or one is ghost and one is not accessible
      // we only increase the counter for the ext id to use based on the
      // number of particles created by other nodes
      if (((!p1 or p1->l.ghost) and (!p2 or p2->l.ghost)) or !p1 or !p2) {
        // Increase local counters
        if (collision_params.mode & COLLISION_MODE_VS) {
          current_vs_pid++;
        }
        // For glue to surface, we have only one vs
        current_vs_pid++;
        if (collision_params.mode == COLLISION_MODE_GLUE_TO_SURF) {
          if (p1)
            if (p1->p.type == collision_params.part_type_to_be_glued) {
              p1->p.type = collision_params.part_type_after_glueing;
            }
          if (p2)
            if (p2->p.type == collision_params.part_type_to_be_glued) {
              p2->p.type = collision_params.part_type_after_glueing;
            }
        } // mode glue to surface

      } else { // We consider the pair because one particle
               // is local to the node and the other is local or ghost
        // If we are in the two vs mode
        // Virtual site related to first particle in the collision
        if (collision_params.mode & COLLISION_MODE_VS) {
          Utils::Vector3d pos1, pos2;

          // Enable rotation on the particles to which vs will be attached
          p1->p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;
          p2->p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;

          // Positions of the virtual sites
          bind_at_point_of_collision_calc_vs_pos(p1, p2, pos1, pos2);

          auto handle_particle = [&](Particle *p, Utils::Vector3d const &pos) {
            if (not p->l.ghost) {
              place_vs_and_relate_to_particle(current_vs_pid, pos,
                                              p->identity());
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

          bind_at_poc_create_bond_between_vs(current_vs_pid, c);
        } // mode VS

        if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF) {
          // If particles are made inert by a type change on collision:
          // We skip the pair if one of the particles has already reacted
          // but we still increase the particle counters, as other nodes
          // can not always know whether or not a vs is placed
          if (collision_params.part_type_after_glueing !=
              collision_params.part_type_to_be_glued) {
            if ((p1->p.type == collision_params.part_type_after_glueing) ||
                (p2->p.type == collision_params.part_type_after_glueing)) {
              current_vs_pid++;
              continue;
            }
          }

          Utils::Vector3d pos;
          const Particle &attach_vs_to =
              glue_to_surface_calc_vs_pos(*p1, *p2, pos);

          // Add a bond between the centers of the colliding particles
          // The bond is placed on the node that has p1
          if (!p1->l.ghost) {
            const int bondG[] = {c.pp2};
            get_part(c.pp1).bonds().insert(
                {collision_params.bond_centers, bondG});
          }

          // Change type of particle being attached, to make it inert
          if (p1->p.type == collision_params.part_type_to_be_glued) {
            p1->p.type = collision_params.part_type_after_glueing;
          }
          if (p2->p.type == collision_params.part_type_to_be_glued) {
            p2->p.type = collision_params.part_type_after_glueing;
          }

          // Vs placement happens on the node that has p1
          if (!attach_vs_to.l.ghost) {
            place_vs_and_relate_to_particle(current_vs_pid, pos,
                                            attach_vs_to.identity());
            // Particle storage locations may have changed due to
            // added particle
            p1 = cell_structure.get_local_particle(c.pp1);
            p2 = cell_structure.get_local_particle(c.pp2);
            current_vs_pid++;
          } else { // Just update the books
            current_vs_pid++;
          }
          glue_to_surface_bind_part_to_vs(p1, p2, current_vs_pid, c);
        }
      } // we considered the pair
    }   // Loop over all collisions in the queue
#ifdef ADDITIONAL_CHECKS
    if (!Utils::Mpi::all_compare(comm_cart, current_vs_pid)) {
      throw std::runtime_error("Nodes disagree about current_vs_pid");
    }
#endif

    // If any node had a collision, all nodes need to resort
    if (!gathered_queue.empty()) {
      cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
      cells_update_ghosts(Cells::DATA_PART_PROPERTIES | Cells::DATA_PART_BONDS);
    }
  }    // are we in one of the vs_based methods
#endif // defined VIRTUAL_SITES_RELATIVE

  // three-particle-binding part
  if (collision_params.mode & (COLLISION_MODE_BIND_THREE_PARTICLES)) {
    auto gathered_queue = gather_global_collision_queue();
    three_particle_binding_domain_decomposition(gathered_queue);
  } // if TPB

  local_collision_queue.clear();
}

#endif
