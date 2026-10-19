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

#pragma once

#include "config/config.hpp"

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"
#include "particle_store/VectorReference.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"
#include "thermostats/brownian_inline.hpp"

#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <cstdint>

// The Brownian sub-kernels (thermostats/brownian_inline.hpp) and the two
// drivers below are templated on a "particle-like" accessor type. The column
// path passes a BrownianRowView built from hoisted *_view() column handles +
// a store row. The row-view exposes the SAME accessor names/semantics the
// drivers and sub-kernels use, but each access indexes the pre-resolved view
// handle by row instead of paying view_host() on a Particle proxy. All
// write-through references (pos/v/quat/omega/torque) alias the columns, so
// the multi-stage quaternion update (bd_drag_rot then bd_random_walk_rot
// reading it back) and every arithmetic step are bitwise identical to the
// Particle-view path.

/** @brief Non-owning column-view of one particle for the Brownian kernels. */
template <class PosView, class VelView, class ForceView, class TorqueView,
          class QuatView, class OmegaView, class RinertiaView, class IdView,
          class MassView, class RotationView, class ExtFlagView,
          class GammaView, class GammaRotView>
class BrownianRowView {
  PosView const &m_pos;
  VelView const &m_vel;
  ForceView const &m_force;
  TorqueView const &m_torque;
  QuatView const &m_quat;
  OmegaView const &m_omega;
  RinertiaView const &m_rinertia;
  IdView const &m_id;
  MassView const &m_mass;
  RotationView const &m_rotation;
  ExtFlagView const &m_ext_flag;
  GammaView const &m_gamma;
  GammaRotView const &m_gamma_rot;
  int m_row;

public:
  BrownianRowView(PosView const &pos, VelView const &vel,
                  ForceView const &force, TorqueView const &torque,
                  QuatView const &quat, OmegaView const &omega,
                  RinertiaView const &rinertia, IdView const &id,
                  MassView const &mass, RotationView const &rotation,
                  ExtFlagView const &ext_flag, GammaView const &gamma,
                  GammaRotView const &gamma_rot, int const row)
      : m_pos(pos), m_vel(vel), m_force(force), m_torque(torque), m_quat(quat),
        m_omega(omega), m_rinertia(rinertia), m_id(id), m_mass(mass),
        m_rotation(rotation), m_ext_flag(ext_flag), m_gamma(gamma),
        m_gamma_rot(gamma_rot), m_row(row) {}

  VectorReference pos() const {
    return VectorReference(const_cast<double *>(&m_pos(m_row, 0)),
                           m_pos.stride(1));
  }
  Utils::Vector3d force() const {
    return {m_force(m_row, 0), m_force(m_row, 1), m_force(m_row, 2)};
  }
  VectorReference v() const {
    return VectorReference(const_cast<double *>(&m_vel(m_row, 0)),
                           m_vel.stride(1));
  }
  int id() const { return m_id(m_row); }

#ifdef ESPRESSO_MASS
  double mass() const { return m_mass(m_row); }
#else
  double mass() const { return 1.0; }
#endif

#ifdef ESPRESSO_EXTERNAL_FORCES
  bool is_fixed_along(unsigned int const axis) const {
    return detail::get_nth_bit(m_ext_flag(m_row), axis);
  }
#else
  bool is_fixed_along(unsigned int const) const { return false; }
#endif

#ifdef ESPRESSO_ROTATION
  std::uint8_t rotation() const { return m_rotation(m_row); }
  bool can_rotate() const { return static_cast<bool>(m_rotation(m_row)); }
  bool can_rotate_around(unsigned int const axis) const {
    return detail::get_nth_bit(m_rotation(m_row), axis);
  }
  Utils::Quaternion<double> quat() const {
    Utils::Quaternion<double> q;
    for (unsigned int j = 0u; j < 4u; ++j)
      q[j] = m_quat(m_row, j);
    return q;
  }
  QuaternionReference quat() {
    return QuaternionReference(const_cast<double *>(&m_quat(m_row, 0)),
                               m_quat.stride(1));
  }
  Utils::Vector3d torque() const {
    return {m_torque(m_row, 0), m_torque(m_row, 1), m_torque(m_row, 2)};
  }
  VectorReference torque() {
    return VectorReference(const_cast<double *>(&m_torque(m_row, 0)),
                           m_torque.stride(1));
  }
  Utils::Vector3d omega() const {
    return {m_omega(m_row, 0), m_omega(m_row, 1), m_omega(m_row, 2)};
  }
  VectorReference omega() {
    return VectorReference(const_cast<double *>(&m_omega(m_row, 0)),
                           m_omega.stride(1));
  }
#endif

#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d rinertia() const {
    return {m_rinertia(m_row, 0), m_rinertia(m_row, 1), m_rinertia(m_row, 2)};
  }
#else
  Utils::Vector3d rinertia() const { return {1., 1., 1.}; }
#endif

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  ParticleStore::GammaValue gamma() const {
    return {m_gamma(m_row, 0), m_gamma(m_row, 1), m_gamma(m_row, 2)};
  }
  ParticleStore::GammaValue gamma_rot() const {
    return {m_gamma_rot(m_row, 0), m_gamma_rot(m_row, 1),
            m_gamma_rot(m_row, 2)};
  }
#else
  ParticleStore::GammaValue gamma() const { return m_gamma(m_row); }
  ParticleStore::GammaValue gamma_rot() const { return m_gamma_rot(m_row); }
#endif
#endif
};

/** @brief Deduce and build a @ref BrownianRowView from hoisted view handles. */
template <class PosView, class VelView, class ForceView, class TorqueView,
          class QuatView, class OmegaView, class RinertiaView, class IdView,
          class MassView, class RotationView, class ExtFlagView,
          class GammaView, class GammaRotView>
inline auto make_brownian_row_view(
    PosView const &pos, VelView const &vel, ForceView const &force,
    TorqueView const &torque, QuatView const &quat, OmegaView const &omega,
    RinertiaView const &rinertia, IdView const &id, MassView const &mass,
    RotationView const &rotation, ExtFlagView const &ext_flag,
    GammaView const &gamma, GammaRotView const &gamma_rot, int const row) {
  return BrownianRowView<PosView, VelView, ForceView, TorqueView, QuatView,
                         OmegaView, RinertiaView, IdView, MassView,
                         RotationView, ExtFlagView, GammaView, GammaRotView>(
      pos, vel, force, torque, quat, omega, rinertia, id, mass, rotation,
      ext_flag, gamma, gamma_rot, row);
}

// Templated drivers: work with a Particle view or a BrownianRowView. The
// torque conversion uses the value-parameter overload so the body compiles
// against either type.
template <class ParticleLike>
inline void brownian_dynamics_propagator(BrownianThermostat const &brownian,
                                         ParticleLike &&p, double time_step,
                                         double kT) {
  auto pos = p.pos();
  pos += bd_drag(brownian.gamma, p, time_step);
  p.v() = bd_drag_vel(brownian.gamma, p);
  pos += bd_random_walk(brownian, p, time_step, kT);
  p.v() += bd_random_walk_vel(brownian, p);
}

#ifdef ESPRESSO_ROTATION
template <class ParticleLike>
inline void brownian_dynamics_rotator(BrownianThermostat const &brownian,
                                      ParticleLike &&p, double time_step,
                                      double kT) {
  if (!p.can_rotate())
    return;
  // Convert the torque to the body frame and apply the axis fix (value overload
  // -- identical to convert_torque_to_body_frame_apply_fix(p) on a Particle).
  p.torque() = convert_torque_to_body_frame_apply_fix(
      Utils::Quaternion<double>(p.quat()), Utils::Vector3d(p.torque()),
      p.rotation());
  auto quat = p.quat();
  quat = bd_drag_rot(brownian.gamma_rotation, p, time_step);
  p.omega() = bd_drag_vel_rot(brownian.gamma_rotation, p);
  quat = bd_random_walk_rot(brownian, p, time_step, kT);
  p.omega() += bd_random_walk_vel_rot(brownian, p);
}
#endif // ESPRESSO_ROTATION
