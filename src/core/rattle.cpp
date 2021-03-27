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

#ifdef BOND_CONSTRAINT

#include "CellStructure.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "interactions.hpp"

#include <utils/constants.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/range/algorithm.hpp>

#include <cmath>

void save_old_pos(const ParticleRange &particles,
                  const ParticleRange &ghost_particles) {
  auto save_pos = [](Particle &p) {
    p.r.p_last_timestep = p.r.p;
  };

  boost::for_each(particles, save_pos);
  boost::for_each(ghost_particles, save_pos);
}

/** Initialize the velocity correction vectors. The correction vectors are
 *  stored in @ref ParticleForce::f "Particle::f::f".
 */
static void init_correction_vector(const ParticleRange &particles,
                                   const ParticleRange &ghost_particles) {
  auto reset_force = [](Particle &p) {
    p.rattle.correction.fill(0);
  };

  boost::for_each(particles, reset_force);
  boost::for_each(ghost_particles, reset_force);
}

/**
 * @brief Add the position correction to particles.
 *
 * @param ia_params Parameters
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool add_pos_corr_vec(RigidBond const &ia_params, Particle &p1,
                             Particle &p2) {
  auto const r_ij = get_mi_vector(p1.r.p, p2.r.p, box_geo);
  auto const r_ij2 = r_ij.norm2();

  if (fabs(1.0 - r_ij2 / ia_params.d2) > ia_params.p_tol) {
    auto const r_ij_t =
        get_mi_vector(p1.r.p_last_timestep, p2.r.p_last_timestep, box_geo);
    auto const r_ij_dot = r_ij_t * r_ij;
    auto const G =
        0.50 * (ia_params.d2 - r_ij2) / r_ij_dot / (p1.p.mass + p2.p.mass);

    auto const pos_corr = G * r_ij_t;
    p1.rattle.correction += pos_corr * p2.p.mass;
    p2.rattle.correction -= pos_corr * p1.p.mass;

    return true;
  }

  return false;
}

/** Compute position corrections */
static void compute_pos_corr_vec(int *repeat_, CellStructure &cs) {
  cs.bond_loop(
      [repeat_](Particle &p1, int bond_id, Utils::Span<Particle *> partners) {
        auto const &iaparams = bonded_ia_params[bond_id];

        if (auto const *bond = boost::get<RigidBond>(&iaparams)) {
          auto const corrected = add_pos_corr_vec(*bond, p1, *partners[0]);
          if (corrected)
            *repeat_ += 1;
        }

        /* Rigid bonds cannot break */
        return false;
      });
}

/** Apply position corrections */
static void app_pos_correction(const ParticleRange &particles) {
  for (auto &p : particles) {
      p.r.p += p.rattle.correction;
      p.m.v += p.rattle.correction;
  }
}

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
    cell_structure.ghosts_reduce_rattle_correction();

    app_pos_correction(particles);
    /* Ghost Positions Update */
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
static bool add_vel_corr_vec(RigidBond const &ia_params, Particle &p1,
                             Particle &p2) {
  auto const v_ij = p1.m.v - p2.m.v;
  auto const r_ij = get_mi_vector(p1.r.p, p2.r.p, box_geo);

  auto const v_proj = v_ij * r_ij;
  if (std::abs(v_proj) > ia_params.v_tol) {
    auto const K = v_proj / ia_params.d2 / (p1.p.mass + p2.p.mass);

    auto const vel_corr = K * r_ij;

    p1.rattle.correction -= vel_corr * p2.p.mass;
    p2.rattle.correction += vel_corr * p1.p.mass;

    return true;
  }

  return false;
}

/** Compute velocity correction vectors */
static void compute_vel_corr_vec(int *repeat_, CellStructure &cs) {
  cs.bond_loop(
      [repeat_](Particle &p1, int bond_id, Utils::Span<Particle *> partners) {
        auto const &iaparams = bonded_ia_params[bond_id];

        if (auto const *bond = boost::get<RigidBond>(&iaparams)) {
          auto const corrected = add_vel_corr_vec(*bond, p1, *partners[0]);
          if (corrected)
            *repeat_ += 1;
        }

        /* Rigid bonds cannot break */
        return false;
      });
}

/** Apply velocity corrections */
static void apply_vel_corr(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.m.v += p.rattle.correction;
  }
}

void correct_vel_shake(CellStructure &cs) {
  cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);

  /* transfer the current forces to r.p_last_timestep of the particle structure
   * so that velocity corrections can be stored temporarily at the f.f member of
   * the particle structure */
  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  bool repeat = true;
  int cnt = 0;
  while (repeat && cnt < SHAKE_MAX_ITERATIONS) {
    init_correction_vector(particles, ghost_particles);
    int repeat_ = 0;
    compute_vel_corr_vec(&repeat_, cs);
    cell_structure.ghosts_reduce_rattle_correction();

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
}

#endif
