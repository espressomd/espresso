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

#include "Particle.hpp"
#include "particle_store/VectorReference.hpp"
#include "rotation.hpp"

#include <cstdint>

// The velocity-Verlet translation steps take per-row VectorReference bundles
// (one row-base pointer + stride, resolved ONCE per particle by the caller via
// p.v()/p.pos()/p.force()) instead of a 2D column view subscripted per
// component. This is the fastest access form: the inner axis loop is pure
// stride-1 pointer arithmetic (vel[j] == base[j] on LayoutRight) with the row
// base computed once, so the compiler does not re-derive row*stride per
// component and does not keep three heavyweight 2D-view handles live across
// the loop. Scalar fallbacks (mass 1.0, fixed byte 0u) are supplied by the
// caller under the same ifdefs the Particle accessors used.

/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
 *  v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+0.5 \Delta t) \f]
 */
inline void
velocity_verlet_propagator_1(VectorReference vel, VectorReference pos,
                             VectorReference const force, double const mass,
                             std::uint8_t const fixed, double time_step) {
  for (unsigned int j = 0; j < 3; j++) {
    if (not detail::get_nth_bit(fixed, j)) {
      /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
      vel[j] += 0.5 * time_step * force[j] / mass;

      /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
       * v(t+0.5*dt) */
      pos[j] += time_step * vel[j];
    }
  }
}

/** Final integration step of the Velocity Verlet integrator
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f]
 */
inline void velocity_verlet_propagator_2(VectorReference vel,
                                         VectorReference const force,
                                         double const mass,
                                         std::uint8_t const fixed,
                                         double time_step) {
  for (unsigned int j = 0; j < 3; j++) {
    if (not detail::get_nth_bit(fixed, j)) {
      /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
      vel[j] += 0.5 * time_step * force[j] / mass;
    }
  }
}

#ifdef ESPRESSO_ROTATION
// The VV-rotation steps are column kernels (rotation.hpp). The caller
// hoists the quaternion / angular-velocity / rinertia / torque / rotation view
// handles ONCE and passes them + the store row; the kernels guard the
// can_rotate() (rotation != 0) case per-row internally.
template <class QuatView, class OmegaView, class RinertiaView, class TorqueView,
          class RotationView>
inline void velocity_verlet_rotator_1(QuatView const &quat_view,
                                      OmegaView const &omega_view,
                                      RinertiaView const &rinertia_view,
                                      TorqueView const &torque_view,
                                      RotationView const &rotation_view,
                                      int const row, double time_step) {
  velocity_verlet_rotator_1_kernel(quat_view, omega_view, rinertia_view,
                                   torque_view, rotation_view, row, time_step);
}

template <class QuatView, class OmegaView, class RinertiaView, class TorqueView,
          class RotationView>
inline void velocity_verlet_rotator_2(QuatView const &quat_view,
                                      OmegaView const &omega_view,
                                      RinertiaView const &rinertia_view,
                                      TorqueView const &torque_view,
                                      RotationView const &rotation_view,
                                      int const row, double time_step) {
  velocity_verlet_rotator_2_kernel(quat_view, omega_view, rinertia_view,
                                   torque_view, rotation_view, row, time_step);
}
#endif // ESPRESSO_ROTATION
