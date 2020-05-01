/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "rattle.hpp"

int n_rigidbonds = 0;

#ifdef BOND_CONSTRAINT

#include "CellStructure.hpp"
#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"

#include <utils/constants.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <cmath>

/** \name Private functions */
/************************************************************/
/*@{*/

/** Positional Corrections are added to the current particle positions. Invoked
 * from \ref correct_pos_shake() */
static void app_pos_correction(const ParticleRange &particles);

/** Transfers temporarily the current forces from f.f[3] of the \ref Particle
    structure to r.p_old[3] location and also initializes velocity correction
    vector. Invoked from \ref correct_vel_shake()*/
static void transfer_force_init_vel(const ParticleRange &particles,
                                    const ParticleRange &ghost_particles);

/** Calculates corrections of the  current particle velocities according to
   RATTLE
    algorithm. Invoked from \ref correct_vel_shake()*/
static void compute_vel_corr_vec(int *repeat_, CellStructure &cs);

/** Velocity corrections are added to the current particle velocities. Invoked
   from
    \ref correct_vel_shake()*/
static void apply_vel_corr(const ParticleRange &particles);

/**Invoked from \ref correct_vel_shake(). Put back the forces from r.p_old to
 * f.f*/
static void revert_force(const ParticleRange &particles,
                         const ParticleRange &ghost_particles);

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
static void init_correction_vector(const ParticleRange &local_particles,
                                   const ParticleRange &ghost_particles) {
  auto reset_force = [](Particle &p) {
    for (int j = 0; j < 3; j++)
      p.f.f[j] = 0.0;
  };

  for (auto &p : local_particles)
    reset_force(p);

  for (auto &p : ghost_particles)
    reset_force(p);
}

/**
 * @brief Add the position correction to particles.
 *
 * The position correction is accumulated in the forces
 * of the particles so that it can be reduced over the ghosts.
 *
 * @param ia_params Parameters
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool add_pos_corr_vec(Bonded_ia_parameters const &ia_params,
                             Particle &p1, Particle &p2) {
  auto const r_ij = get_mi_vector(p1.r.p, p2.r.p, box_geo);
  auto const r_ij2 = r_ij.norm2();

  if (fabs(1.0 - r_ij2 / ia_params.p.rigid_bond.d2) >
      ia_params.p.rigid_bond.p_tol) {
    auto const r_ij_t = get_mi_vector(p1.r.p_old, p2.r.p_old, box_geo);
    auto const r_ij_dot = r_ij_t * r_ij;
    auto const G = 0.50 * (ia_params.p.rigid_bond.d2 - r_ij2) / r_ij_dot /
                   (p1.p.mass + p2.p.mass);

    auto const pos_corr = G * r_ij_t;
    p1.f.f += pos_corr * p2.p.mass;
    p2.f.f -= pos_corr * p1.p.mass;

    return true;
  }

  return false;
}
/**Compute positional corrections*/
static void compute_pos_corr_vec(int *repeat_, CellStructure &cs) {
  for (auto &p1 : cs.local_particles()) {
    cs.execute_bond_handler(p1, [repeat_](Particle &p1, int bond_id,
                                          Utils::Span<Particle *> partners) {
      auto const &iaparams = bonded_ia_params[bond_id];

      if (iaparams.type == BONDED_IA_RIGID_BOND) {
        auto const corrected = add_pos_corr_vec(iaparams, p1, *partners[0]);
        if (corrected)
          *repeat_ += 1;
      }

      /* Rigid bonds cannot break */
      return false;
    });
  }
}

/**Apply corrections to each particle**/
static void app_pos_correction(const ParticleRange &particles) {
  /*Apply corrections*/
  for (auto &p : particles) {
    for (int j = 0; j < 3; j++) {
      p.r.p[j] += p.f.f[j];
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

/** Calculates the corrections required for each of the particle coordinates
    according to the RATTLE algorithm. Invoked from \ref correct_pos_shake()*/
void correct_pos_shake(CellStructure &cs) {
  cells_update_ghosts(Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES);

  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  int cnt = 0;
  bool repeat = true;

  while (repeat && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector(particles, ghost_particles);
    int repeat_ = 0;
    compute_pos_corr_vec(&repeat_, cs);
    cell_structure.ghosts_reduce_forces();

    app_pos_correction(particles);
    /**Ghost Positions Update*/
    cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);

    repeat = boost::mpi::all_reduce(comm_cart, (repeat_ > 0),
                                    std::logical_or<bool>());

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
static void transfer_force_init_vel(const ParticleRange &particles,
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

/**
 * @brief Add the velocity correction to particles.
 *
 * The position correction is accumulated in the forces
 * of the particles so that it can be reduced over the ghosts.
 *
 * @param ia_params Parameters
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool add_vel_corr_vec(Bonded_ia_parameters const &ia_params,
                             Particle &p1, Particle &p2) {
  auto const v_ij = p1.m.v - p2.m.v;
  auto const r_ij = get_mi_vector(p1.r.p, p2.r.p, box_geo);

  auto const v_proj = v_ij * r_ij;
  if (std::abs(v_proj) > ia_params.p.rigid_bond.v_tol) {
    auto const K = v_proj / ia_params.p.rigid_bond.d2 / (p1.p.mass + p2.p.mass);

    auto const vel_corr = K * r_ij;

    p1.f.f -= vel_corr * p2.p.mass;
    p2.f.f += vel_corr * p1.p.mass;

    return true;
  }

  return false;
}

/** Velocity correction vectors are computed*/
static void compute_vel_corr_vec(int *repeat_, CellStructure &cs) {
  for (auto &p1 : cs.local_particles()) {
    cs.execute_bond_handler(p1, [repeat_](Particle &p1, int bond_id,
                                          Utils::Span<Particle *> partners) {
      auto const &iaparams = bonded_ia_params[bond_id];

      if (iaparams.type == BONDED_IA_RIGID_BOND) {
        auto const corrected = add_vel_corr_vec(iaparams, p1, *partners[0]);
        if (corrected)
          *repeat_ += 1;
      }

      /* Rigid bonds can not break */
      return false;
    });
  }
}

/**Apply velocity corrections*/
static void apply_vel_corr(const ParticleRange &particles) {
  /*Apply corrections*/
  for (auto &p : particles) {
    for (int j = 0; j < 3; j++) {
      p.m.v[j] += p.f.f[j];
    }
    /**Completed for one particle*/
  } // for i loop
}

/**Put back the forces from r.p_old to f.f*/
static void revert_force(const ParticleRange &particles,
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

void correct_vel_shake(CellStructure &cs) {
  cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);

  /**transfer the current forces to r.p_old of the particle structure so that
  velocity corrections can be stored temporarily at the f.f[3] of the particle
  structure  */
  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  transfer_force_init_vel(particles, ghost_particles);

  bool repeat = true;
  int cnt = 0;
  while (repeat && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector(particles, ghost_particles);
    int repeat_ = 0;
    compute_vel_corr_vec(&repeat_, cs);
    cell_structure.ghosts_reduce_forces();
    apply_vel_corr(particles);
    cs.ghosts_update(Cells::DATA_PART_MOMENTUM);

    repeat = boost::mpi::all_reduce(comm_cart, (repeat_ > 0),
                                    std::logical_or<bool>());
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
