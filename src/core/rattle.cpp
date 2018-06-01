/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#include "domain_decomposition.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

int n_rigidbonds = 0;

#ifdef BOND_CONSTRAINT

/** \name Private functions */
/************************************************************/
/*@{*/

/** Calculates the corrections required for each of the particle coordinates
    according to the RATTLE algorithm. Invoked from \ref correct_pos_shake()*/
void compute_pos_corr_vec(int *repeat_);

/** Positional Corrections are added to the current particle positions. Invoked
 * from \ref correct_pos_shake() */
void app_pos_correction();

/** Transfers temporarily the current forces from f.f[3] of the \ref Particle
    structure to r.p_old[3] location and also intializes velocity correction
    vector. Invoked from \ref correct_vel_shake()*/
void transfer_force_init_vel();

/** Calculates corrections of the  current particle velocities according to
   RATTLE
    algorithm. Invoked from \ref correct_vel_shake()*/
void compute_vel_corr_vec(int *repeat_);

/** Velocity corrections are added to the current particle velocities. Invoked
   from
    \ref correct_vel_shake()*/
void apply_vel_corr();

/**Invoked from \ref correct_vel_shake(). Put back the forces from r.p_old to
 * f.f*/
void revert_force();

/**For debugging purpose--prints the bond lengths between particles that have
rigid_bonds*/
void print_bond_len();

/*@}*/

/*Initialize old positions (particle positions at previous time step)
  of the particles*/
void save_old_pos() {
  auto save_pos = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.r.p_old[j] = p.r.p[j];
  };

  for (auto &p : local_cells.particles())
    save_pos(p);

  for (auto &p : ghost_cells.particles())
    save_pos(p);
}

/**Initialize the correction vector. The correction vector is stored in f.f of
 * particle strcuture. */
void init_correction_vector() {
  auto reset_force = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.f.f[j] = 0.0;
  };

  for (auto &p : local_cells.particles())
    reset_force(p);

  for (auto &p : ghost_cells.particles())
    reset_force(p);
}

/**Compute positional corrections*/
void compute_pos_corr_vec(int *repeat_) {
  Bonded_ia_parameters *ia_params;
  int j, k, cnt = -1;
  Particle *p1, *p2;
  double r_ij_t[3], r_ij[3], r_ij_dot, G, pos_corr, r_ij2;

  for (auto &p : local_cells.particles()) {
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

        get_mi_vector(r_ij, p1->r.p, p2->r.p);
        r_ij2 = sqrlen(r_ij);
        if (fabs(1.0 - r_ij2 / ia_params->p.rigid_bond.d2) >
            ia_params->p.rigid_bond.p_tol) {
          get_mi_vector(r_ij_t, p1->r.p_old, p2->r.p_old);
          r_ij_dot = scalar(r_ij_t, r_ij);
          G = 0.50 * (ia_params->p.rigid_bond.d2 - r_ij2) / r_ij_dot;
#ifdef MASS
          G /= ((*p1).p.mass + (*p2).p.mass);
#else
          G /= 2;
#endif
          for (j = 0; j < 3; j++) {
            pos_corr = G * r_ij_t[j];
            p1->f.f[j] += pos_corr * (*p2).p.mass;
            p2->f.f[j] -= pos_corr * (*p1).p.mass;
          }
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
void app_pos_correction() {
  /*Apply corrections*/
  for (auto &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      p.r.p[j] += p.f.f[j];
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

void correct_pos_shake() {
  int repeat_, cnt = 0;
  int repeat = 1;

  while (repeat != 0 && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector();
    repeat_ = 0;
    compute_pos_corr_vec(&repeat_);
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    app_pos_correction();
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
void transfer_force_init_vel() {
  auto copy_reset = [](Particle &p) {
    for (int j = 0; j < 3; j++) {
      p.r.p_old[j] = p.f.f[j];
      p.f.f[j] = 0.0;
    }
  };

  for (auto &p : local_cells.particles())
    copy_reset(p);

  for (auto &p : ghost_cells.particles())
    copy_reset(p);
}

/** Velocity correction vectors are computed*/
void compute_vel_corr_vec(int *repeat_) {
  Bonded_ia_parameters *ia_params;
  int j, k;
  Particle *p1, *p2;
  double v_ij[3], r_ij[3], K, vel_corr;

  for (auto &p : local_cells.particles()) {
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

        vecsub(p1->m.v, p2->m.v, v_ij);
        get_mi_vector(r_ij, p1->r.p, p2->r.p);
        if (fabs(scalar(v_ij, r_ij)) > ia_params->p.rigid_bond.v_tol) {
          K = scalar(v_ij, r_ij) / ia_params->p.rigid_bond.d2;
#ifdef MASS
          K /= ((*p1).p.mass + (*p2).p.mass);
#else
          K /= 2.0;
#endif
          for (j = 0; j < 3; j++) {
            vel_corr = K * r_ij[j];
            p1->f.f[j] -= vel_corr * (*p2).p.mass;
            p2->f.f[j] += vel_corr * (*p1).p.mass;
          }
          *repeat_ = *repeat_ + 1;
        }
      } else
        k += ia_params->num;
    } // while loop
  }   // for i loop
}

/**Apply velocity corrections*/
void apply_vel_corr() {
  /*Apply corrections*/
  for (auto &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

/**Put back the forces from r.p_old to f.f*/
void revert_force() {
  auto revert = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.f.f[j] = p.r.p_old[j];
  };

  for (auto &p : local_cells.particles())
    revert(p);

  for (auto &p : ghost_cells.particles())
    revert(p);
}

void correct_vel_shake() {
  int repeat_, repeat = 1, cnt = 0;
  /**transfer the current forces to r.p_old of the particle structure so that
  velocity corrections can be stored temporarily at the f.f[3] of the particle
  structure  */
  transfer_force_init_vel();
  while (repeat != 0 && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector();
    repeat_ = 0;
    compute_vel_corr_vec(&repeat_);
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    apply_vel_corr();
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    if (this_node == 0)
      MPI_Reduce(&repeat_, &repeat, 1, MPI_INT, MPI_SUM, 0, comm_cart);
    else
      MPI_Reduce(&repeat_, nullptr, 1, MPI_INT, MPI_SUM, 0, comm_cart);

    MPI_Bcast(&repeat, 1, MPI_INT, 0, comm_cart);
    cnt++;
  }

  if (cnt >= SHAKE_MAX_ITERATIONS) {
    fprintf(stderr, "%d: VEL CORRECTIONS IN RATTLE failed to converge after %d "
                    "iterations !!\n",
            this_node, cnt);
    errexit();
  }
  /**Puts back the forces from r.p_old to f.f[3]*/
  revert_force();
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
  bonded_ia_params[bond_type].p.rigid_bond.v_tol = v_tol * time_step;
  bonded_ia_params[bond_type].type = BONDED_IA_RIGID_BOND;
  bonded_ia_params[bond_type].num = 1;
  n_rigidbonds += 1;
  mpi_bcast_ia_params(bond_type, -1);
  mpi_bcast_parameter(FIELD_RIGIDBONDS);

  return ES_OK;
}

#endif
