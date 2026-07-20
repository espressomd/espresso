/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#ifdef ESPRESSO_BOND_CONSTRAINT

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/range/algorithm.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <variant>
#include <vector>

/** Maximal number of iterations before the RATTLE algorithm bails out. */
static constexpr auto shake_max_iterations = 1000;

static void check_convergence(int cnt, char const *const name) {
  static constexpr char const *const msg = " failed to converge after ";
  if (cnt >= shake_max_iterations) {
    runtimeErrorMsg() << name << msg << cnt << " iterations";
  }
}

/**
 * @brief copy current position
 *
 * @param particles particle range
 * @param ghost_particles ghost particle range
 */
void save_old_position(const ParticleRange &particles,
                       const ParticleRange &ghost_particles) {
  auto save_pos = [](Particle &p) { p.pos_last_time_step() = p.pos(); };

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
  auto reset_force = [](Particle &p) { p.rattle_params().correction.fill(0); };

  boost::for_each(particles, reset_force);
  boost::for_each(ghost_particles, reset_force);
}

/**
 * @brief Calculate the positional correction for the particles.
 *
 * @param ia_params Parameters
 * @param box_geo Box geometry.
 * @param p1 First particle.
 * @param p2 Second particle.
 * @param bond_id Bonded interaction id of this specific bond.
 * @param rigid_bond_virial Per-bond-type RATTLE constraint virial
 *                          accumulator, indexed by @p bond_id.
 * @return True if there was a correction.
 */
static bool calculate_positional_correction(
    RigidBond const &ia_params, BoxGeometry const &box_geo, Particle &p1,
    Particle &p2, int bond_id,
    std::vector<Utils::Vector9d> &rigid_bond_virial) {
  auto const r_ij = box_geo.get_mi_vector(p1.pos(), p2.pos());
  auto const r_ij2 = r_ij.norm2();

  if (std::abs(1.0 - r_ij2 / ia_params.d2) > ia_params.p_tol) {
    auto const r_ij_t =
        box_geo.get_mi_vector(p1.pos_last_time_step(), p2.pos_last_time_step());
    auto const r_ij_dot = r_ij_t * r_ij;
    auto const G =
        0.50 * (ia_params.d2 - r_ij2) / r_ij_dot / (p1.mass() + p2.mass());

    auto const pos_corr = G * r_ij_t;
    p1.rattle_params().correction += pos_corr * p2.mass();
    p2.rattle_params().correction -= pos_corr * p1.mass();

    // Constraint force implied by this bond alone during this iteration:
    // the correction just applied to p1 is
    // @f$ \Delta r1 = pos_corr*m2 = (1/2)*a1*dt^2 @f$,
    // so @f$ F1 = m1*a1 = 2*m1*\Delta r1/dt^2 = 2*m1*m2*pos_corr/dt^2 @f$.
    // Division by dt^2 is deferred to @ref System::calculate_pressure(),
    // where the timestep is available. This uses r_ij_t, the bond vector
    // at the start of the MD step (fixed across all SHAKE iterations of this
    // step), so the contribution is exact for this bond alone, regardless of
    // how many other rigid bonds p1 or p2 participate in.
    rigid_bond_virial[static_cast<std::size_t>(bond_id)] += Utils::flatten(
        Utils::tensor_product(2.0 * p1.mass() * p2.mass() * pos_corr, r_ij_t));

    return true;
  }

  return false;
}

/**
 * @brief Compute the correction vectors using given kernel.
 *
 * @param cs cell structure
 * @param box_geo Box geometry
 * @param bonded_ias Bonded interactions
 * @param kernel kernel function
 * @return True if correction is necessary
 */
template <typename Kernel>
static bool compute_correction_vector(CellStructure &cs,
                                      BoxGeometry const &box_geo,
                                      BondedInteractionsMap const &bonded_ias,
                                      Kernel kernel) {
  bool correction = false;
  cs.bond_loop([&correction, &kernel, &box_geo, &bonded_ias](
                   Particle &p1, int bond_id, std::span<Particle *> partners) {
    auto const &iaparams = *bonded_ias.at(bond_id);

    if (auto const *bond = std::get_if<RigidBond>(&iaparams)) {
      auto const corrected = kernel(*bond, box_geo, p1, *partners[0], bond_id);
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
    p.pos() += p.rattle_params().correction;
    p.v() += p.rattle_params().correction;
  });
}

void correct_position_shake(CellStructure &cs, BoxGeometry const &box_geo,
                            BondedInteractionsMap &bonded_ias) {
  unsigned const flag = Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES;
  cs.update_ghosts_and_resort_particle(flag);

  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  // Reset the per-bond-type constraint virial for this timestep
  bonded_ias.rigid_bond_virial.assign(
      static_cast<std::size_t>(bonded_ias.get_next_key()),
      Utils::Vector9d::broadcast(0.));

  int cnt;
  for (cnt = 0; cnt < shake_max_iterations; ++cnt) {
    init_correction_vector(particles, ghost_particles);
    bool const repeat_ = compute_correction_vector(
        cs, box_geo, bonded_ias,
        [&bonded_ias](RigidBond const &bond, BoxGeometry const &box_geo_,
                      Particle &p1, Particle &p2, int bond_id) {
          return calculate_positional_correction(
              bond, box_geo_, p1, p2, bond_id, bonded_ias.rigid_bond_virial);
        });
    bool const repeat =
        boost::mpi::all_reduce(comm_cart, repeat_, std::logical_or<bool>());

    // no correction is necessary, skip communication and bail out
    if (!repeat)
      break;

    cs.ghosts_reduce_rattle_correction();

    apply_positional_correction(particles);
    cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);
  }
  check_convergence(cnt, "RATTLE");

  auto const resort_level =
      cs.check_resort_required() ? Cells::RESORT_LOCAL : Cells::RESORT_NONE;
  cs.set_resort_particles(resort_level);
}

/**
 * @brief Calculate the velocity correction for the particles.
 *
 * @param ia_params Parameters
 * @param box_geo Box geometry.
 * @param p1 First particle.
 * @param p2 Second particle.
 * @return True if there was a correction.
 */
static bool calculate_velocity_correction(RigidBond const &ia_params,
                                          BoxGeometry const &box_geo,
                                          Particle &p1, Particle &p2) {
  auto const v_ij = p1.v() - p2.v();
  auto const r_ij = box_geo.get_mi_vector(p1.pos(), p2.pos());

  auto const v_proj = v_ij * r_ij;
  if (std::abs(v_proj) > ia_params.v_tol) {
    auto const K = v_proj / ia_params.d2 / (p1.mass() + p2.mass());

    auto const vel_corr = K * r_ij;

    p1.rattle_params().correction -= vel_corr * p2.mass();
    p2.rattle_params().correction += vel_corr * p1.mass();

    return true;
  }

  return false;
}

/**
 * @brief Apply velocity corrections
 *
 * @param particles particle range
 */
static void apply_velocity_correction(ParticleRange const &particles) {
  boost::for_each(particles,
                  [](Particle &p) { p.v() += p.rattle_params().correction; });
}

void correct_velocity_shake(CellStructure &cs, BoxGeometry const &box_geo,
                            BondedInteractionsMap const &bonded_ias) {
  cs.ghosts_update(Cells::DATA_PART_POSITION | Cells::DATA_PART_MOMENTUM);

  auto particles = cs.local_particles();
  auto ghost_particles = cs.ghost_particles();

  int cnt;
  for (cnt = 0; cnt < shake_max_iterations; ++cnt) {
    init_correction_vector(particles, ghost_particles);
    bool const repeat_ = compute_correction_vector(
        cs, box_geo, bonded_ias,
        [](RigidBond const &bond, BoxGeometry const &box_geo_, Particle &p1,
           Particle &p2, int /* bond_id */) {
          return calculate_velocity_correction(bond, box_geo_, p1, p2);
        });
    bool const repeat =
        boost::mpi::all_reduce(comm_cart, repeat_, std::logical_or<bool>());

    // no correction is necessary, skip communication and bail out
    if (!repeat)
      break;

    cs.ghosts_reduce_rattle_correction();

    apply_velocity_correction(particles);
    cs.ghosts_update(Cells::DATA_PART_MOMENTUM);
  }
  check_convergence(cnt, "VEL RATTLE");
}

#endif
