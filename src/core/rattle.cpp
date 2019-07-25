/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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

#include "rattle.hpp"

int n_rigidbonds = 0;

#ifdef BOND_CONSTRAINT

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

/** \name Private functions */
/************************************************************/
/*@{*/

/** Calculates the corrections required for each of the particle coordinates
    according to the RATTLE algorithm. Invoked from \ref correct_pos_shake()*/
void compute_pos_corr_vec(int *repeat_, const ParticleRange &particles);

/** Positional Corrections are added to the current particle positions. Invoked
 * from \ref correct_pos_shake() */
void app_pos_correction(const ParticleRange &particles);

/** Transfers temporarily the current forces from f.f[3] of the \ref Particle
    structure to r.p_old[3] location and also initializes velocity correction
    vector. Invoked from \ref correct_vel_shake()*/
void transfer_force_init_vel(const ParticleRange &particles,
                             const ParticleRange &ghost_particles);

/** Calculates corrections of the  current particle velocities according to
   RATTLE
    algorithm. Invoked from \ref correct_vel_shake()*/
void compute_vel_corr_vec(int *repeat_, const ParticleRange &particles);

/** Velocity corrections are added to the current particle velocities. Invoked
   from
    \ref correct_vel_shake()*/
void apply_vel_corr(const ParticleRange &particles);

/**Invoked from \ref correct_vel_shake(). Put back the forces from r.p_old to
 * f.f*/
void revert_force(const ParticleRange &particles,
                  const ParticleRange &ghost_particles);

/**For debugging purpose--prints the bond lengths between particles that have
rigid_bonds*/
void print_bond_len();

/*@}*/

/*Initialize old positions (particle positions at previous time step)
  of the particles*/
void save_old_pos(const ParticleRange &particles,
                  const ParticleRange &ghost_particles) {
  auto save_pos = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.r.p_old[j] = p.r.p[j];
  };

  for (auto &p : particles)
    save_pos(p);

  for (auto &p : ghost_particles)
    save_pos(p);
}

/**Initialize the correction vector. The correction vector is stored in f.f of
 * particle structure. */
void init_correction_vector(const ParticleRange &particles) {
  auto reset_force = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.f.f[j] = 0.0;
  };

  for (auto &p : particles)
    reset_force(p);

  for (auto &p : ghost_cells.particles())
    reset_force(p);
}

/**Compute positional corrections*/
void compute_pos_corr_vec(int *repeat_, const ParticleRange &particles) {
  Bonded_ia_parameters *ia_params;
  int j, k, cnt = -1;
  Particle *p1, *p2;

  for (auto &p : particles) {
    p1 = &p;
    k = 0;
    while (k < p1->bl.n) {
      ia_params = &bonded_ia_params[p1->bl.e[k++]];
      if (ia_params->type == BONDED_IA_RIGID_BOND) {
        cnt++;
        p2 = local_particles[p1->bl.e[k++]];
        if (!p2) {
          runtimeErrorMsg() << "rigid bond broken between particles "
                            << p1->p.identity << " and " << p1->bl.e[k - 1]
                            << " (particles not stored on the same node)";
          return;
        }

        auto const r_ij = get_mi_vector(p1->r.p, p2->r.p, box_geo);
        auto const r_ij2 = r_ij.norm2();

        if (fabs(1.0 - r_ij2 / ia_params->p.rigid_bond.d2) >
            ia_params->p.rigid_bond.p_tol) {
          auto const r_ij_t = get_mi_vector(p1->r.p_old, p2->r.p_old, box_geo);
          auto const r_ij_dot = r_ij_t * r_ij;
          auto const G = 0.50 * (ia_params->p.rigid_bond.d2 - r_ij2) /
                         r_ij_dot / (p1->p.mass + p2->p.mass);

          auto const pos_corr = G * r_ij_t;
          p1->f.f += pos_corr * p2->p.mass;
          p2->f.f -= pos_corr * p1->p.mass;

          /*Increase the 'repeat' flag by one */
          *repeat_ = *repeat_ + 1;
        }
      } else
        /* skip bond partners of nonrigid bond */
        k += ia_params->num;
    } // while loop
  }   // for i loop
}

/**Apply corrections to each particle**/
void app_pos_correction(const ParticleRange &particles) {
  /*Apply corrections*/
  for (auto &p : particles) {
    for (int j = 0; j < 3; j++) {
      p.r.p[j] += p.f.f[j];
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

void correct_pos_shake(const ParticleRange &particles) {
  int repeat_, cnt = 0;
  int repeat = 1;

  while (repeat != 0 && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector(cell_structure.local_cells().particles());
    repeat_ = 0;
    compute_pos_corr_vec(&repeat_, cell_structure.local_cells().particles());
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    app_pos_correction(particles);
    /**Ghost Positions Update*/
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    if (this_node == 0)
      MPI_Reduce(&repeat_, &repeat, 1, MPI_INT, MPI_SUM, 0, comm_cart);
    else
      MPI_Reduce(&repeat_, nullptr, 1, MPI_INT, MPI_SUM, 0, comm_cart);
    MPI_Bcast(&repeat, 1, MPI_INT, 0, comm_cart);

    cnt++;
  } // while(repeat) loop
  if (cnt >= SHAKE_MAX_ITERATIONS) {
    runtimeErrorMsg() << "RATTLE failed to converge after " << cnt
                      << " iterations";
  }

  check_resort_particles();
}

/**The forces are transfered temporarily from f.f member of particle structure
   to r.p_old,
    which is idle now and initialize the velocity correction vector to zero at
   f.f[3]
    of Particle structure*/
void transfer_force_init_vel(const ParticleRange &particles,
                             const ParticleRange &ghost_particles) {
  auto copy_reset = [](Particle &p) {
    for (int j = 0; j < 3; j++) {
      p.r.p_old[j] = p.f.f[j];
      p.f.f[j] = 0.0;
    }
  };

  for (auto &p : particles)
    copy_reset(p);

  for (auto &p : ghost_particles)
    copy_reset(p);
}

/** Velocity correction vectors are computed*/
void compute_vel_corr_vec(int *repeat_, const ParticleRange &particles) {
  Bonded_ia_parameters *ia_params;
  int j, k;
  Particle *p1, *p2;

  for (auto &p : particles) {
    p1 = &p;
    k = 0;
    while (k < p1->bl.n) {
      ia_params = &bonded_ia_params[p1->bl.e[k++]];
      if (ia_params->type == BONDED_IA_RIGID_BOND) {
        p2 = local_particles[p1->bl.e[k++]];
        if (!p2) {
          runtimeErrorMsg() << "rigid bond broken between particles "
                            << p1->p.identity << " and " << p1->bl.e[k - 1]
                            << " (particles not stored on the same node)";
          return;
        }

        auto const v_ij = p1->m.v - p2->m.v;
        auto const r_ij = get_mi_vector(p1->r.p, p2->r.p, box_geo);

        auto const v_proj = v_ij * r_ij;
        if (std::abs(v_proj) > ia_params->p.rigid_bond.v_tol) {
          auto const K =
              v_proj / ia_params->p.rigid_bond.d2 / (p1->p.mass + p2->p.mass);

          auto const vel_corr = K * r_ij;

          p1->f.f -= vel_corr * p2->p.mass;
          p2->f.f += vel_corr * p1->p.mass;

          *repeat_ = *repeat_ + 1;
        }
      } else
        k += ia_params->num;
    } // while loop
  }   // for i loop
}

/**Apply velocity corrections*/
void apply_vel_corr(const ParticleRange &particles) {
  /*Apply corrections*/
  for (auto &p : particles) {
    for (int j = 0; j < 3; j++) {
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

/**Put back the forces from r.p_old to f.f*/
void revert_force(const ParticleRange &particles,
                  const ParticleRange &ghost_particles) {
  auto revert = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.f.f[j] = p.r.p_old[j];
  };

  for (auto &p : particles)
    revert(p);

  for (auto &p : ghost_particles)
    revert(p);
}

void correct_vel_shake(CellStructure &cell_structure) {
  int repeat_, repeat = 1, cnt = 0;
  /**transfer the current forces to r.p_old of the particle structure so that
  velocity corrections can be stored temporarily at the f.f[3] of the particle
  structure  */
  auto particles = cell_structure.local_cells().particles();
  auto ghost_particles = cell_structure.ghost_cells().particles();

  transfer_force_init_vel(particles, ghost_particles);
  while (repeat != 0 && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector(particles);
    repeat_ = 0;
    compute_vel_corr_vec(&repeat_, cell_structure.local_cells().particles());
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    apply_vel_corr(particles);
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    if (this_node == 0)
      MPI_Reduce(&repeat_, &repeat, 1, MPI_INT, MPI_SUM, 0, comm_cart);
    else
      MPI_Reduce(&repeat_, nullptr, 1, MPI_INT, MPI_SUM, 0, comm_cart);

    MPI_Bcast(&repeat, 1, MPI_INT, 0, comm_cart);
    cnt++;
  }

  if (cnt >= SHAKE_MAX_ITERATIONS) {
    fprintf(stderr,
            "%d: VEL CORRECTIONS IN RATTLE failed to converge after %d "
            "iterations !!\n",
            this_node, cnt);
    errexit();
  }
  /**Puts back the forces from r.p_old to f.f[3]*/
  revert_force(particles, ghost_particles);
}

/*****************************************************************************
 *   setting parameters
 *****************************************************************************/
int rigid_bond_set_params(int bond_type, double d, double p_tol, double v_tol) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.rigid_bond.d2 = d * d;
  bonded_ia_params[bond_type].p.rigid_bond.p_tol = 2.0 * p_tol;
  bonded_ia_params[bond_type].p.rigid_bond.v_tol = v_tol;
  bonded_ia_params[bond_type].type = BONDED_IA_RIGID_BOND;
  bonded_ia_params[bond_type].num = 1;
  n_rigidbonds += 1;
  mpi_bcast_ia_params(bond_type, -1);
  mpi_bcast_parameter(FIELD_RIGIDBONDS);

  return ES_OK;
}

#endif
