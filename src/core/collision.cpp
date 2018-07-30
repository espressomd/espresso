/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>

#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "particle_data.hpp"
#include "collision.hpp"
#include "initialize.hpp"
#include "interaction_data.hpp"
#include "rotation.hpp"
#include "virtual_sites/VirtualSitesRelative.hpp"

#ifdef COLLISION_DETECTION_DEBUG
#define TRACE(a) a
#else
#define TRACE(a)
#endif

#ifdef COLLISION_DETECTION

/// Data type holding the info about a single collision
typedef struct {
  int pp1; // 1st particle id
  int pp2; // 2nd particle id
} collision_struct;

// During force calculation, colliding particles are recorded in the queue
// The queue is processed after force calculation, when it is save to add
// particles
static std::vector<collision_struct> local_collision_queue;

/// Parameters for collision detection
Collision_parameters collision_params;

/** @brief Return true if a bond between the centers of the colliding particles
 * needs to be placed. At this point, all modes need this */
inline bool bind_centers() {
  return collision_params.mode != COLLISION_MODE_OFF;
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
    if (collision_params.distance > min_global_cut) {
      runtimeErrorMsg() << "The minimum global cutoff (System.min_global_cut) "
                           "must be larger or equal the collision detection "
                           "distance.";
    }
  }

#ifndef VIRTUAL_SITES_RELATIVE
  // The collision modes involving virutal istes also requires the creation of a
  // bond between the colliding
  // If we don't have virtual sites, virtual site binding isn't possible.
  if ((collision_params.mode & COLLISION_MODE_VS) ||
      (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    runtimeErrorMsg() << "Virtual sites based collisoin modes modes require "
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

  // For vs based methods, Binding so far only works on a single cpu
  //  if ((collision_params.mode & COLLISION_MODE_VS) ||(collision_params.mode &
  //  COLLISION_MODE_GLUE_TO_SURF))
  //    if (n_nodes != 1) {
  //      runtimeErrorMsg() << "Virtual sites based collision modes only work on
  //      a single node.";
  //      return false;
  //    }
  //
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
          << "Collision detection particle type after glueing needs to be >=0";
      return false;
    }
    if (this_node == 0)
      make_particle_type_exist(collision_params.part_type_after_glueing);
  }

  recalc_forces = 1;
  rebuild_verletlist = 1;

  return true;
}

//* Allocate memory for the collision queue /
void prepare_local_collision_queue() { local_collision_queue.clear(); }

void queue_collision(const int part1, const int part2) {
  local_collision_queue.push_back({part1, part2});
}

/** @brief Calculate position of vs for GLUE_TO_SURFACE mode
*    Reutnrs id of particle to bind vs to */
int glue_to_surface_calc_vs_pos(const Particle *const p1,
                                const Particle *const p2, Vector3d& pos) {
  int bind_vs_to_pid;
  double vec21[3];
  double c;
  get_mi_vector(vec21, p1->r.p, p2->r.p);
  const double dist_betw_part = sqrt(sqrlen(vec21));

  // Find out, which is the particle to be glued.
  if ((p1->p.type == collision_params.part_type_to_be_glued) &&
      (p2->p.type == collision_params.part_type_to_attach_vs_to)) {
    c = 1 - collision_params.dist_glued_part_to_vs / dist_betw_part;
    bind_vs_to_pid = p2->p.identity;
  } else if ((p2->p.type == collision_params.part_type_to_be_glued) &&
             (p1->p.type == collision_params.part_type_to_attach_vs_to)) {
    c = collision_params.dist_glued_part_to_vs / dist_betw_part;
    bind_vs_to_pid = p1->p.identity;
  } else {
    throw std::runtime_error("This should never be thrown. Bug.");
  }
  for (int i = 0; i < 3; i++) {
    pos[i] = p2->r.p[i] + vec21[i] * c;
  }
  return bind_vs_to_pid;
}

void bind_at_point_of_collision_calc_vs_pos(const Particle *const p1,
                                            const Particle *const p2,
                                            Vector3d& pos1, Vector3d& pos2) {
  double vec21[3];
  get_mi_vector(vec21, p1->r.p, p2->r.p);
  for (int i = 0; i < 3; i++) {
    pos1[i] = p1->r.p[i] - vec21[i] * collision_params.vs_placement;
    pos2[i] = p1->r.p[i] - vec21[i] * (1. - collision_params.vs_placement);
  }
}

// Considers three particles for three_particle_binding and performs
// the binding if the criteria are met //
void coldet_do_three_particle_bond(Particle &p, Particle &p1, Particle &p2) {
  double vec21[3];
  // If p1 and p2 are not closer or equal to the cutoff distance, skip
  // p1:
  get_mi_vector(vec21, p.r.p, p1.r.p);
  if (sqrt(sqrlen(vec21)) > collision_params.distance)
    return;
  // p2:
  get_mi_vector(vec21, p.r.p, p2.r.p);
  if (sqrt(sqrlen(vec21)) > collision_params.distance)
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
  double cosine = 0.0;

  double vec1[3], vec2[3];
  /* vector from p to p1 */
  get_mi_vector(vec1, p.r.p, p1.r.p);
  // Normalize
  double dist2 = sqrlen(vec1);
  double d1i = 1.0 / sqrt(dist2);
  for (int j = 0; j < 3; j++)
    vec1[j] *= d1i;

  /* vector from p to p2 */
  get_mi_vector(vec2, p.r.p, p2.r.p);
  // normalize
  dist2 = sqrlen(vec2);
  double d2i = 1.0 / sqrt(dist2);
  for (int j = 0; j < 3; j++)
    vec2[j] *= d2i;

  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);

  // Handle case where cosine is nearly 1 or nearly -1
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;

  // Bond angle
  double phi = acos(cosine);

  // We find the bond id by dividing the range from 0 to pi in
  // three_particle_angle_resolution steps and by adding the id
  // of the bond for zero degrees.
  int bond_id =
      floor(phi / M_PI *
                (collision_params.three_particle_angle_resolution - 1) +
            0.5) +
      collision_params.bond_three_particles;

  // Create the bond

  // First, fill bond data structure
  int bondT[3];
  bondT[0] = bond_id;
  bondT[1] = p1.p.identity;
  bondT[2] = p2.p.identity;
  local_change_bond(p.p.identity, bondT, 0);
}

#ifdef VIRTUAL_SITES_RELATIVE
void place_vs_and_relate_to_particle(const int current_vs_pid,
                                     const Vector3d& pos,
                                     const int relate_to,
                                     const Vector3d& initial_pos) {

  // The virtual site is placed at initial_pos which will be in the local
  // node's domain. It will then be moved to its final position.
  // A resort occurs after vs-based collisions anyway, which will move the vs
  // into the right cell.
  added_particle(current_vs_pid);
  local_place_particle(current_vs_pid, initial_pos.data(), 1);
  local_particles[current_vs_pid]->r.p=pos;
  local_vs_relate_to(current_vs_pid, relate_to);

  (local_particles[max_seen_particle])->p.is_virtual = 1;
  (local_particles[max_seen_particle])->p.type =
      collision_params.vs_particle_type;
}

void bind_at_poc_create_bond_between_vs(const int current_vs_pid,
                                        const collision_struct &c) {
  int bondG[3];

  switch (bonded_ia_params[collision_params.bond_vs].num) {
  case 1: {
    // Create bond between the virtual particles
    bondG[0] = collision_params.bond_vs;
    bondG[1] = current_vs_pid - 2;
    local_change_bond(current_vs_pid - 1, bondG, 0);
    break;
  }
  case 2: {
    // Create 1st bond between the virtual particles
    bondG[0] = collision_params.bond_vs;
    bondG[1] = c.pp1;
    bondG[2] = c.pp2;
    local_change_bond(current_vs_pid - 1, bondG, 0);
    local_change_bond(current_vs_pid - 2, bondG, 0);
    break;
  }
  }
}

void glue_to_surface_bind_part_to_vs(const Particle *const p1,
                                     const Particle *const p2,
                                     const int vs_pid_plus_one,
                                     const collision_struct &c) {
  int bondG[3];
  // Create bond between the virtual particles
  bondG[0] = collision_params.bond_vs;
  bondG[1] = vs_pid_plus_one - 1;
  if (p1->p.type == collision_params.part_type_after_glueing) {
    local_change_bond(p1->p.identity, bondG, 0);
  } else {
    local_change_bond(p2->p.identity, bondG, 0);
  }
}

#endif

std::vector<collision_struct> gather_global_collision_queue() {
  std::vector<collision_struct> res;

  std::vector<int> displacements(n_nodes); // offsets into collisions

  // Initialize number of collisions gathered from all processors
  std::vector<int> counts(n_nodes);
  // Total number of collisions
  int total_collisions;
  int local_queue_size = local_collision_queue.size();
  MPI_Allreduce(&local_queue_size, &total_collisions, 1, MPI_INT, MPI_SUM,
                comm_cart);

  if (total_collisions == 0)
    return res;

  // Gather number of collisions
  MPI_Allgather(&local_queue_size, 1, MPI_INT, &(counts[0]), 1, MPI_INT,
                comm_cart);

  // initialize displacement information for all nodes
  displacements[0] = 0;

  // Find where to place collision information for each processor
  std::vector<int> byte_counts(n_nodes);
  for (int k = 1; k < n_nodes; k++)
    displacements[k] =
        displacements[k - 1] + (counts[k - 1]) * sizeof(collision_struct);

  for (int k = 0; k < n_nodes; k++)
    byte_counts[k] = counts[k] * sizeof(collision_struct);

  // Allocate mem for the new collision info

  res.resize(total_collisions);

  // Gather collision informtion from all nodes and send it to all nodes
  MPI_Allgatherv(&(local_collision_queue[0]), byte_counts[this_node], MPI_BYTE,
                 &(res[0]), &(byte_counts[0]), &(displacements[0]), MPI_BYTE,
                 comm_cart);

  return res;
}

// this looks in all local particles for a particle close to those in a
// 2-particle collision. If it finds them, it performs three particle binding
void three_particle_binding_full_search(
    const std::vector<collision_struct> &gathered_queue) {
  // Handler checking a single particle against collision queue. Used for
  // local particles and ghosts
  auto handle_single_particle = [&gathered_queue](Particle &p) {
    for (auto &c : gathered_queue) {
      Particle *p1 = local_particles[c.pp1];
      Particle *p2 = local_particles[c.pp2];

      // Check, whether p is equal to one of the particles in the
      // collision. If so, skip
      if ((p.p.identity == p1->p.identity) ||
          (p.p.identity == p2->p.identity)) {
        continue;
      }

      // The following checks,
      // if the particle p is closer that the cutoff from p1 and/or p2.
      // If yes, three particle bonds are created on all particles
      // which have two other particles within the cutoff distance,
      // unless such a bond already exists

      // We need all cyclical permutations, here
      // (bond is placed on 1st particle, order of bond partners
      // does not matter, so we don't neet non-cyclic permutations):
      coldet_do_three_particle_bond(p, *p1, *p2);
      coldet_do_three_particle_bond(*p1, p, *p2);
      coldet_do_three_particle_bond(*p2, p, *p1);
    }
  };

  // Iterate over all local particles
  for (auto &p : local_cells.particles()) {
    handle_single_particle(p);
  };

  // And ghosts
  for (auto &p : ghost_cells.particles()) {
    handle_single_particle(p);
  }
}

static void three_particle_binding_dd_do_search(
    Cell *basecell, Particle &p1, Particle &p2) {
  int basecellno = std::distance(&cells[0], basecell);
  for (int n = 0; n < 27; ++n) {
    Cell &cell = cells[dd_full_shell_neigh(basecellno, n)];
    for (int pno = 0; pno < cell.n; ++pno) {
      Particle &P = cell.part[pno];
      // Skip collided particles themselves
      if ((P.p.identity == p1.p.identity) ||
          (P.p.identity == p2.p.identity)) {
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
    if ((local_particles[c.pp1] != NULL) && (local_particles[c.pp2] != NULL)) {

      Particle &p1 = *local_particles[c.pp1];
      Particle &p2 = *local_particles[c.pp2];
      auto cell1 = cell_structure.position_to_cell(p1.r.p.data());
      auto cell2 = cell_structure.position_to_cell(p2.r.p.data());

      three_particle_binding_dd_do_search(cell1, p1, p2);
      if (cell1 != cell2)
        three_particle_binding_dd_do_search(cell2, p1, p2);

    } // If local particles exist

  } // Loop over total collisions
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

  if (bind_centers()) {
    for (auto &c : local_collision_queue) {
      // put the bond to the non-ghost particle; at least one partner always is
      if (local_particles[c.pp1]->l.ghost) {
        std::swap(c.pp1, c.pp2);
      }
      int bondG[2];
      bondG[0] = collision_params.bond_centers;
      bondG[1] = c.pp2;
      local_change_bond(c.pp1, bondG, 0);
    }
  }

// Virtual sites based collision schemese
#ifdef VIRTUAL_SITES_RELATIVE
  if ((collision_params.mode & COLLISION_MODE_VS) ||
      (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF)) {
    // Gather the global collision queue, because only one node hsa a collision
    // across node boundaries in its cueue.
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

      // Get particle pointers
      Particle *const p1 = local_particles[c.pp1];
      Particle *const p2 = local_particles[c.pp2];

      // If we have none of the two partic;es, only increase the counter for the
      // next id to use
      if (p1 == NULL and p2 == NULL) {
        if (collision_params.mode & COLLISION_MODE_VS) {
          added_particle(current_vs_pid);
          current_vs_pid++;
        }
        // For glue to surface, we have only one vs
        added_particle(current_vs_pid);
        current_vs_pid++;

      } else { // We have at least one

        // Calculate initial position for new vs, which is in the local node's
        // domain
        // Vs is moved afterwards and resorted after all collision s are handled
        // Use position of non-ghost colliding particle.
        Vector3d initial_pos;
        if (p1->l.ghost)
          initial_pos = p2->r.p;
        else
          initial_pos= p1->r.p;

        // If we are in the two vs mode
        // Virtual site related to first particle in the collision
        if (collision_params.mode & COLLISION_MODE_VS) {
          Vector3d pos1, pos2;

          // Enable rotation on the particles to which vs will be attached
          p1->p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;
          p2->p.rotation = ROTATION_X | ROTATION_Y | ROTATION_Z;

          // The vs placement is done by the node on which p1 is non-ghost
          if (!p1->l.ghost) {
            bind_at_point_of_collision_calc_vs_pos(p1, p2, pos1, pos2);
            place_vs_and_relate_to_particle(current_vs_pid, pos1, c.pp1,
                                            initial_pos);
            current_vs_pid++;
            place_vs_and_relate_to_particle(current_vs_pid, pos2, c.pp2,
                                            initial_pos);
            current_vs_pid++;
            bind_at_poc_create_bond_between_vs(current_vs_pid, c);
          } else { // Just update the books
            added_particle(current_vs_pid);
            current_vs_pid++;

            added_particle(current_vs_pid);
            current_vs_pid++;
          }
        }
        if (collision_params.mode & COLLISION_MODE_GLUE_TO_SURF) {
          Vector3d pos;
          const int pid = glue_to_surface_calc_vs_pos(p1, p2, pos);

          // Change type of partilce being attached, to make it inert
          if (p1->p.type == collision_params.part_type_to_be_glued) {
            p1->p.type = collision_params.part_type_after_glueing;
          }
          if (p2->p.type == collision_params.part_type_to_be_glued) {
            p2->p.type = collision_params.part_type_after_glueing;
          }

          // Vs placement happens on the node that has p1
          if (!p1->l.ghost) {
            place_vs_and_relate_to_particle(current_vs_pid, pos, pid,
                                            initial_pos);
            current_vs_pid++;
          } else { // Just update the books
            added_particle(current_vs_pid);
            current_vs_pid++;
          }
          glue_to_surface_bind_part_to_vs(p1, p2, current_vs_pid, c);
        }
      } // Loop over all collisions in the queue

      // If any node had a collision, all nodes need to do on_particle_change
      // and resort

      if (gathered_queue.size() > 0) {
        on_particle_change();
        announce_resort_particles();
        cells_update_ghosts();
      }
    }  // total_collision>0
  }    // are we in one of the vs_based methods
#endif // defined VIRTUAL_SITES_RELATIVE

  // three-particle-binding part
  if (collision_params.mode & (COLLISION_MODE_BIND_THREE_PARTICLES)) {
    auto gathered_queue = gather_global_collision_queue();

    // If we don't have domain decomposition, we need to do a full sweep over
    // all
    // particles in the system. (slow)
    if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
      three_particle_binding_full_search(gathered_queue);
    } // if cell structure != domain decomposition
    else {
      three_particle_binding_domain_decomposition(gathered_queue);
    } // If we have doamin decomposition

  } // if TPB

  local_collision_queue.clear();
}

#endif
