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
#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/range/algorithm.hpp>

/**
 * @brief copy current position
 *
 * @param particles particle range
 * @param ghost_particles ghost particle range
 */
void save_old_position(const ParticleRange &particles,
                       const ParticleRange &ghost_particles) {
  auto save_pos = [](Particle &p) { p.r.p_last_timestep = p.r.p; };

  boost::for_each(particles, save_pos);
  boost::for_each(ghost_particles, save_pos);
}

/**
 * @brief reset correction vectors to zero
 *
 * @param particles particle range
 * @param ghost_particles ghost particle range
 */
static void init_correction_vector(const ParticleRange &particles,
                                   const ParticleRange &ghost_particles) {
  auto reset_force = [](Particle &p) { p.rattle.correction.fill(0); };

  boost::for_each(particles, reset_force);
  boost::for_each(ghost_particles, reset_force);
}

/**
 * @brief Calculate the positional correction for the particles.
 *
 * @param ia_params Parameters
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool calculate_positional_correction(RigidBond const &ia_params,
                                            Particle &p1, Particle &p2) {
  auto const r_ij = box_geo.get_mi_vector(p1.r.p, p2.r.p);
  auto const r_ij2 = r_ij.norm2();

  if (std::abs(1.0 - r_ij2 / ia_params.d2) > ia_params.p_tol) {
    auto const r_ij_t =
        box_geo.get_mi_vector(p1.r.p_last_timestep, p2.r.p_last_timestep);
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

/**
 * @brief Compute the correction vectors using given kernel.
 *
 * @param cs cell structure
 * @param kernel kernel function
 * @return True if correction is necessary
 */
template <typename Kernel>
static bool compute_correction_vector(CellStructure &cs, Kernel kernel) {
  bool correction = false;
  cs.bond_loop([&correction, &kernel](Particle &p1, int bond_id,
                                      Utils::Span<Particle *> partners) {
    auto const &iaparams = *bonded_ia_params.at(bond_id);

    if (auto const *bond = boost::get<RigidBond>(&iaparams)) {
      auto const corrected = kernel(*bond, p1, *partners[0]);
      if (corrected)
        correction = true;
    }

    /* Rigid bonds cannot break */
    return false;
  });

  return correction;
}

/**
 * @brief Apply positional corrections
 *
 * @param particles particle range
 */
static void apply_positional_correction(const ParticleRange &particles) {
  boost::for_each(particles, [](Particle &p) {
    p.r.p += p.rattle.correction;
    p.m.v += p.rattle.correction;
  });
}

void correct_position_shake(CellStructure &cs) {
  cells_update_ghosts(Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES);

  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  int cnt;
  for (cnt = 0; cnt < SHAKE_MAX_ITERATIONS; ++cnt) {
    init_correction_vector(particles, ghost_particles);
    bool const repeat_ =
        compute_correction_vector(cs, calculate_positional_correction);
    bool const repeat =
        boost::mpi::all_reduce(comm_cart, repeat_, std::logical_or<bool>());

    // no correction is necessary, skip communication and bail out
    if (!repeat)
      break;

    cell_structure.ghosts_reduce_rattle_correction();

    apply_positional_correction(particles);
    cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);
  }
  if (cnt >= SHAKE_MAX_ITERATIONS) {
    runtimeErrorMsg() << "RATTLE failed to converge after " << cnt
                      << " iterations";
  }

  check_resort_particles();
}

/**
 * @brief Calculate the velocity correction for the particles.
 *
 * The position correction is accumulated in the forces
 * of the particles so that it can be reduced over the ghosts.
 *
 * @param ia_params Parameters
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool calculate_velocity_correction(RigidBond const &ia_params,
                                          Particle &p1, Particle &p2) {
  auto const v_ij = p1.m.v - p2.m.v;
  auto const r_ij = box_geo.get_mi_vector(p1.r.p, p2.r.p);

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

/**
 * @brief Apply velocity corrections
 *
 * @param particles particle range
 */
static void apply_velocity_correction(const ParticleRange &particles) {
  boost::for_each(particles, [](Particle &p) { p.m.v += p.rattle.correction; });
}

void correct_velocity_shake(CellStructure &cs) {
  cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);

  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  int cnt;
  for (cnt = 0; cnt < SHAKE_MAX_ITERATIONS; ++cnt) {
    init_correction_vector(particles, ghost_particles);
    bool const repeat_ =
        compute_correction_vector(cs, calculate_velocity_correction);
    bool const repeat =
        boost::mpi::all_reduce(comm_cart, repeat_, std::logical_or<bool>());

    // no correction is necessary, skip communication and bail out
    if (!repeat)
      break;

    cell_structure.ghosts_reduce_rattle_correction();

    apply_velocity_correction(particles);
    cs.ghosts_update(Cells::DATA_PART_MOMENTUM);
  }

  if (cnt >= SHAKE_MAX_ITERATIONS) {
    runtimeErrorMsg() << "VEL RATTLE failed to converge after " << cnt
                      << " iterations";
  }
}

#endif
