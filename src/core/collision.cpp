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
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "rotation.hpp"
#include "virtual_sites/VirtualSitesRelative.hpp"

#include <utils/mpi/all_compare.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/algorithm/clamp.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/serialization.hpp>

#include <vector>

/// Data type holding the info about a single collision
typedef struct {
  std::vector<int> pp; // particle id vector
  std::vector<int> idx; // index in collision_params vectors
} collision_struct;

namespace boost {
namespace serialization {
template <typename Archive>
void serialize(Archive &ar, collision_struct &c, const unsigned int) {
  ar &c.pp[0];
  ar &c.pp[1];
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
  return assert(local_particles[id]), *local_particles[id];
}
} // namespace


bool validate_collision_parameters() {
  // If collision detection is deactivated, no further checks
  if (!collision_params.active) {
    return true;
  }

  // Validate distance
  if (collision_params.distance <= 0.) {
    runtimeErrorMsg() << "collision_detection distance must be > 0.";
    return false;
  }

  // Cache square of cutoff
  collision_params.distance_cutoff = Utils::sqr(collision_params.distance);

  if (collision_params.distance > min_global_cut) {
    runtimeErrorMsg() << "The minimum global cutoff (System.min_global_cut) "
                         "must be larger or equal the collision detection "
                         "distance.";
    return false;
  }

  if (collision_params.rate < 0.) {
    runtimeErrorMsg() << "The rate at which a collision will form a bond needs "
                         "to be lager than 0.";
    return false;
  }

  if (collision_params.particle_type.size() <= 1) {
    runtimeErrorMsg() << "At least two, not necessary different, particle types "
                         "are needed fot the collision detection.";
    return false;
  }

  if (collision_params.particle_type.size() < 
      collision_params.particle_type_after_collision.size()) {
    runtimeErrorMsg() << "The size of the vs_particle_type_after_collision array "
                         "needs to be equal or less than the size of the "
                         "particle_type array.";
    return false;
  }

#ifndef VIRTUAL_SITES_RELATIVE
  // Check if VIRTUAL_SITES is needed
  if (!collision_params.vs_particle_type.empty()) {
    runtimeErrorMsg() << "Virtual sites based collisions require "
                         "the VIRTUAL_SITES feature";
    return false;
  }
#endif

  // Check vs parameters
  if (collision_params.particle_type.size() < 
      collision_params.vs_particle_type.size()) {
    runtimeErrorMsg() << "The size of the vs_particle_type array needs to be "
                         "equal or less than the size of the particle_type array.";
    return false;
  }

  if (collision_params.vs_particle_type.size() != 
      collision_params.distance_vs_particle.size()) {
    runtimeErrorMsg() << "Every virtual site created needs a distance to "
                         "its real particle specified (vs_particle_type array "
                         "does not match distance_vs_particle array).";
    return false;
  }

  for (auto d: collision_params.distance_vs_particle) {
    if(d < 0.0 or d > 1.0){
      runtimeErrorMsg() << "The distance of the vitual sites amd the particle "
                           "needs to be between 0 and 1.";
      return false;
    }
  }

  if (collision_params.bond_type < 0) {
    runtimeErrorMsg() << "The collision bond type needs to be defined.";
    return false;
  }

  if (collision_params.bond_type >= bonded_ia_params.size()) {
    runtimeErrorMsg() << "The bond type to be used for binding particle "
                         "centers does not exist.";
    return false;
  }

  if (!collision_params.vs_particle_type.size() >= 1) {
    if (collision_params.vs_bond_type < 0) {
      runtimeErrorMsg() << "The collision bond type of the virtual sites needs "
                           "to be defined when creating virtual sites.";
      return false;
    }
    if (collision_params.vs_bond_type >= bonded_ia_params.size()) {
    runtimeErrorMsg() << "The bond type to be used for binding the virtual"
                         "particles does not exist.";
    return false;
    }
  }

  recalc_forces = true;
  rebuild_verletlist = true;

  return true;
}

void prepare_local_collision_queue() { local_collision_queue.clear(); }

void queue_collision(const std::vector<int> particles) {
  collision_struct c;
  //indices of the vector entries
  std::vector<int> indices;

  for(auto part = particles.begin(); part != particles.end(); part++){
    auto it = std::find(collision_params.particle_type.begin(),
                        collision_params.particle_type.end(), local_particles[*part]->p.type);

    indices.push_back(std::distance(collision_params.particle_type.begin(), it));
  }
  c.pp = particles;
  c.idx = indices;
  local_collision_queue.push_back(c);
}

void bind_at_point_of_collision_calc_vs_pos(std::vector<Particle*> p,
                                            std::vector<int> indices,
                                            std::vector<Utils::Vector3d> &pos) {
  auto const vec = get_mi_vector(p[0]->r.p, p[1]->r.p, box_geo);
  pos[0] = p[0]->r.p - vec * collision_params.distance_vs_particle[indices[0]];
  if(pos.size()>1)
    pos[1] = p[1]->r.p + vec * collision_params.distance_vs_particle[indices[1]];
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

  // Iterate over existing bonds of p

  if (p.bl.e) {
    int b = 0;
    while (b < p.bl.n) {
      int size = bonded_ia_params[p.bl.e[b]].num;

      if (size == 2) {
        // Check if the bond type is within the range used by the collision
        // detection,
        if ((p.bl.e[b] >= collision_params.bond_three_particles) &
            (p.bl.e[b] <=
             collision_params.bond_three_particles +
                 collision_params.three_particle_angle_resolution)) {
          // check, if p1 and p2 are the bond partners, (in any order)
          // if yes, skip triplet
          if (((p.bl.e[b + 1] == p1.p.identity) &&
               (p.bl.e[b + 2] == p2.p.identity)) ||
              ((p.bl.e[b + 1] == p2.p.identity) &&
               (p.bl.e[b + 2] == p1.p.identity)))
            return;
        } // if bond type
      }   // if size==2

      // Go to next bond
      b += size + 1;
    } // bond loop
  }   // if bond list defined

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
  const Utils::Vector3i bondT = {bond_id, p1.p.identity, p2.p.identity};

  local_add_particle_bond(p, bondT);
}

#ifdef VIRTUAL_SITES_RELATIVE
void place_vs_and_relate_to_particle(const int current_vs_pid,
                                     const Utils::Vector3d &pos,
                                     int relate_to) {
  added_particle(current_vs_pid);
  Particle new_part;
  new_part.p.identity = current_vs_pid;
  new_part.r.p = pos;
  auto p_vs = append_indexed_particle(cell_structure.local_cells().cell[0],
                                      std::move(new_part));

  local_vs_relate_to(*p_vs, get_part(relate_to));

  p_vs->p.is_virtual = true;
  p_vs->p.type = collision_params.vs_particle_type[0];//TODO
}

void bind_at_poc_create_bond_between_vs(const int current_vs_pid/*,
                                        const collision_struct &c */) {
  switch (bonded_ia_params[collision_params.vs_bond_type].num) {
  case 1: {
    // Create bond between the virtual particles
    const int bondG[] = {collision_params.vs_bond_type, current_vs_pid - 2};
    // Only add bond if vs was created on this node
    if (local_particles[current_vs_pid - 1])
      local_add_particle_bond(get_part(current_vs_pid - 1), bondG);
    break;
  }
// TODO
/*  case 2: {
    // Create 1st bond between the virtual particles
    const int bondG[] = {collision_params.bond_vs, c.pp1, c.pp2};
    // Only add bond if vs was created on this node
    if (local_particles[current_vs_pid - 1])
      local_add_particle_bond(get_part(current_vs_pid - 1), bondG);
    if (local_particles[current_vs_pid - 2])
      local_add_particle_bond(get_part(current_vs_pid - 2), bondG);
    break;
  }*/
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
    for (int p_id = 0; p_id < c->n; p_id++) {
      auto &P = c->part[p_id];

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
    if ((local_particles[c.pp[0]]) && (local_particles[c.pp[1]])) {
      Particle &p1 = *local_particles[c.pp[0]];
      Particle &p2 = *local_particles[c.pp[1]];
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
                          << std::min(c.pp[0], c.pp[1]) << " and "
                          << std::max(c.pp[0], c.pp[1]);
    }
  }

  // Note that the glue to surface mode adds bonds between the centers
  // but does so later in the process. This is needed to guarantee that
  // a particle can only be glued once, even if queued twice in a single
  // time step
  for (auto &c : local_collision_queue) {
    // put the bond to the non-ghost particle; at least one partner always is
    if (local_particles[c.pp[0]]->l.ghost) {
      std::swap(c.pp[0], c.pp[1]);
    }
    int bondG[2];
    bondG[0] = collision_params.bond_type;
    bondG[1] = c.pp[1];
    local_add_particle_bond(get_part(c.pp[0]), bondG);

    // Change the particle types
    for(int i=0; i<c.pp.size(); i++){
      if(collision_params.particle_type_after_collision[c.idx[i]])
        local_particles[c.pp[i]]->p.type = collision_params.particle_type_after_collision[c.idx[i]];
    }
  }

// Virtual sites based collision schemes
#ifdef VIRTUAL_SITES_RELATIVE
  if (!collision_params.vs_particle_type.empty()) {
    // Gather the global collision queue, because only one node has a collision
    // across node boundaries in its queue.
    // The other node might still have to change particle properties on its
    // non-ghost particle
    auto gathered_queue = gather_global_collision_queue();

    // Sync max_seen_part
    MPI_Allreduce(MPI_IN_PLACE, &max_seen_particle, 1, MPI_INT, MPI_MAX,
                  comm_cart);

    // Make sure, the local_particles array is long enough
    realloc_local_particles(max_seen_particle);

    int current_vs_pid = max_seen_particle + 1;

    // Iterate over global collision queue
    for (auto &c : gathered_queue) {
        if(collision_params.vs_particle_type[c.idx[0]] ||
           collision_params.vs_particle_type[c.idx[1]]){

        // Get particle pointers
        std::vector<Particle *> particles = {local_particles[c.pp[0]], local_particles[c.pp[1]]};

        // If we cannot access both particles, both are ghosts,
        // ore one is ghost and one is not accessible
        // we only increase the counter for the ext id to use based on the
        // number of particles created by other nodes
        if (((!particles[0] or particles[0]->l.ghost) and 
             (!particles[1] or particles[1]->l.ghost)) or 
            !particles[0] or !particles[1]) {
          // Increase local counters
          added_particle(current_vs_pid);
          current_vs_pid++;
          if(collision_params.vs_particle_type.size() > 1){
            // only increase the particle counter further if more than one vs is created
            added_particle(current_vs_pid);
            current_vs_pid++;
          }
          // change ghost particle types
          for(int i = 0; i<particles.size(); i++){
            if (particles[i])
              if (collision_params.particle_type_after_collision[c.idx[i]]) {
                particles[i]->p.type = collision_params.particle_type_after_collision[c.idx[i]];
              }
          }

        } else { // We consider the pair because one particle
                 // is local to the node and the other is local or ghost

          if(!collision_params.vs_particle_type[c.idx[0]]){
            std::swap(particles[0],particles[1]);
            std::swap(c.idx[0],c.idx[1]);
            std::swap(c.pp[0],c.pp[1]);
          }

          // Check how many virtual sides we need
          std::vector<Utils::Vector3d> pos;
          for(auto index : c.idx)
            if(collision_params.vs_particle_type[index])
              pos.push_back({0.0,0.0,0.0});

          // Create the virtual sides
          if(pos.size() <= 2){
            bind_at_point_of_collision_calc_vs_pos(particles, c.idx, pos);
          }else{
            throw std::runtime_error("Unknown behavior of the colliding paticles. Bug.");
          }

          auto handle_particle = [&](Particle *p, Utils::Vector3d const &pos) {
            if (not p->l.ghost) {
              place_vs_and_relate_to_particle(current_vs_pid, pos,
                                              p->identity());
              // Particle storage locations may have changed due to
              // added particle
              particles[0] = local_particles[c.pp[0]];
              particles[1] = local_particles[c.pp[1]];
            } else {
              added_particle(current_vs_pid);
            }
          };

          for (int i=0; i<pos.size(); i++){
            // place virtual sites on the node where the base particle is not a
            // ghost
            handle_particle(particles[i], pos[i]);
            // Increment counter
            current_vs_pid++;
          }


          if(pos.size() == 2){
            // Create a bonde between the virtual sides
            bind_at_poc_create_bond_between_vs(current_vs_pid);
          }else if(pos.size() == 1){
            int bond[] = {collision_params.vs_bond_type, current_vs_pid - 1};
            local_add_particle_bond(*particles[0], bond);
          }
        } // we considered the pair
      }
    }   // Loop over all collisions in the queue
#ifdef ADDITIONAL_CHECKS
    if (!Utils::Mpi::all_compare(comm_cart, current_vs_pid)) {
      throw std::runtime_error("Nodes disagree about current_vs_pid");
    }
    if (!Utils::Mpi::all_compare(comm_cart, max_seen_particle)) {
      throw std::runtime_error("Nodes disagree about max_seen_particle");
    }
#endif

    // If any node had a collision, all nodes need to resort
    if (!gathered_queue.empty()) {
      set_resort_particles(Cells::RESORT_GLOBAL);
      cells_update_ghosts(GHOSTTRANS_PROPRTS | GHOSTTRANS_BONDS);
    }
  }    // are we in one of the vs_based methods
#endif // defined VIRTUAL_SITES_RELATIVE
  /* TODO
  // three-particle-binding part
  if (collision_params.mode & (COLLISION_MODE_BIND_THREE_PARTICLES)) {
    auto gathered_queue = gather_global_collision_queue();
    three_particle_binding_domain_decomposition(gathered_queue);
  } // if TPB
  */
  local_collision_queue.clear();
}

#endif
