/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_NPT

#include "PropagationMode.hpp"
#include "PropagationPredicate.hpp"
#include "npt.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/math/sqr.hpp>

struct PropagationPredicateNPT {
  int modes;
  PropagationPredicateNPT(int default_propagation) {
    modes = PropagationMode::TRANS_LANGEVIN_NPT;
    if (default_propagation & PropagationMode::TRANS_LANGEVIN_NPT) {
      modes |= PropagationMode::SYSTEM_DEFAULT;
    }
  }

  bool operator()(int prop) const { return (prop & modes); }
};

using ParticleRangeNPT = ParticleRangeFiltered<PropagationPredicateNPT>;

/**
 * @brief Special propagator for velocity Verlet NpT with the Andersen method.
 * Propagate the velocities and positions. Integration steps before force
 * calculation of the Velocity Verlet integrator:
 * \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t \cdot F(t)/m \f]
 * \f[ x(t+\Delta t) = x(t) + \Delta t \cdot v(t+0.5 \Delta t) \f]
 *
 * Propagate pressure, box_length (2 times) and positions, rescale
 * positions and velocities and check Verlet list criterion (only NpT).
 */
void velocity_verlet_npt_Andersen_step_1(ParticleRangeNPT const &particles,
                                         IsotropicNptThermostat const &npt_iso,
                                         double time_step,
                                         System::System &system);
/**
 * @brief Special propagator for velocity Verlet NpT with the Andersen method.
 */
void velocity_verlet_npt_MTK_step_1(ParticleRangeNPT const &particles,
                                    IsotropicNptThermostat const &npt_iso,
                                    double time_step, System::System &system);

/**
 * @brief Final integration step of the velocity Verlet NpT integrator
 * with the Andersen method.
 * Finalize instantaneous pressure calculation:
 * \f[ v(t+\Delta t) = v(t+0.5 \Delta t)
 *      + 0.5 \Delta t \cdot F(t+\Delta t)/m \f]
 */
void velocity_verlet_npt_Andersen_step_2(ParticleRangeNPT const &particles,
                                         double time_step,
                                         System::System &system);
/**
 * @brief Final integration step of the velocity Verlet NpT integrator
 * with the MTK method.
 */
void velocity_verlet_npt_MTK_step_2(ParticleRangeNPT const &particles,
                                    double time_step, System::System &system);

/**
 * @brief Propagate the particle's velocity.
 * @f$ v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) @f$
 */
inline void velocity_verlet_npt_propagate_vel_final(
    NptIsoParameters const &nptiso, InstantaneousPressure &npt_inst_pressure,
    ParticleRangeNPT const &particles, double time_step) {

  npt_inst_pressure.p_vel = {};
  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          npt_inst_pressure.p_vel[j] += Utils::sqr(p.v()[j]) * p.mass();
        }
        p.v()[j] += p.force()[j] * time_step / (2. * p.mass());
      }
    }
  }
}

/**
 * @brief Propagate the particle's velocity.
 * @f$ v(t+0.5*dt) = v(t) + 0.5*dt * a(t) @f$
 */
inline void velocity_verlet_npt_propagate_vel(
    NptIsoParameters const &nptiso, InstantaneousPressure &npt_inst_pressure,
    ParticleRangeNPT const &particles, double time_step) {
  npt_inst_pressure.p_vel = {};

  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        p.v()[j] += p.force()[j] * time_step / (2. * p.mass());
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          npt_inst_pressure.p_vel[j] += Utils::sqr(p.v()[j]) * p.mass();
        }
      }
    }
  }
}

#endif // ESPRESSO_NPT
