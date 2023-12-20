/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef NPT

#include "ParticleRange.hpp"
#include "PropagationMode.hpp"
#include "PropagationPredicate.hpp"

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

/** Special propagator for NpT isotropic.
 *  Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator:
 *  \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t \cdot F(t)/m \f]
 *  \f[ x(t+\Delta t) = x(t) + \Delta t \cdot v(t+0.5 \Delta t) \f]
 *
 *  Propagate pressure, box_length (2 times) and positions, rescale
 *  positions and velocities and check Verlet list criterion (only NpT).
 */
void velocity_verlet_npt_step_1(ParticleRangeNPT const &particles,
                                double time_step);

/** Final integration step of the Velocity Verlet+NpT integrator.
 *  Finalize instantaneous pressure calculation:
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t)
 *      + 0.5 \Delta t \cdot F(t+\Delta t)/m \f]
 */
void velocity_verlet_npt_step_2(ParticleRangeNPT const &particles,
                                double time_step);

#endif // NPT
