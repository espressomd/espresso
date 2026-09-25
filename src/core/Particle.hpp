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

#include "BondList.hpp"
#include "PropagationMode.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/quaternion.hpp>

#include <boost/container/vector.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace detail {
constexpr inline bool get_nth_bit(uint8_t const bitfield,
                                  unsigned int const bit_idx) {
  return bitfield & (1u << bit_idx);
}
} // namespace detail

#ifdef ESPRESSO_ENGINE
/** Properties of a self-propelled particle. */
struct ParticleParametersSwimming {
  /** Imposed constant force. */
  double f_swim = 0.;
  /** Is the particle a swimmer. */
  bool swimming = false;
  /** Whether f_swim is applied to the particle or to the fluid. */
  bool is_engine_force_on_fluid = false;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & f_swim & swimming & is_engine_force_on_fluid;
  }
};
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
/** Properties for thermal Stoner-Wohlfarth magnetodynamics. */
struct ThermalStonerWohlfarthParameters {
  /**
   * Flag to distinguish virtual particles carrying the dipole moment in
   * the thermal Stoner-Wohlfarth model from other types of virtual sites.
   */
  bool is_enabled = false;
  /** angle between the director and dipole moment of a Stoner-Wohlfarth
   * particle */
  double phi0 = 0.;
  /** saturation magnetisation of a polarizable particle */
  double sat_mag = 0.;
  /**
   * @brief Inverse anisotropy field in reduced units.
   * anisotropy field = 2.*K1/(mu0 * Ms) in [A / m] where K1 is the magnetic
   * anisotropy constant [kg / (m s^2)].
   */
  double ani_fld_inv = 0.;
  /**
   * @brief Magnetic anisotropy energy (K1 * V) in units of energy.
   * Related to ani_param from Eq.3 in @cite mostarac25a
   * by: ani_param = ani_energy / kT
   */
  double ani_energy = 0.;
  /** Browns attempt frequency. Prefactor from Eq.9 in @cite mostarac25a.  */
  double tau0_inv = 0.;
  /** time units parameter for the kinetic Monte Carlo step */
  double dt_incr = 0.;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & is_enabled & phi0 & sat_mag & ani_fld_inv & ani_energy & tau0_inv &
        dt_incr;
  }
};
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 */
struct ParticleProperties {
  /** unique identifier for the particle. */
  int identity = -1;
  /** Molecule identifier. */
  int mol_id = 0;
  /** particle type, used for non-bonded interactions. */
  int type = 0;
  /** which propagation schemes should be applied to the particle **/
  int propagation = PropagationMode::SYSTEM_DEFAULT;

#ifdef ESPRESSO_ROTATION
  /** Bitfield for the particle axes of rotation.
   *  Values:
   *  - 0: no rotation
   *  - 1: allow rotation around the x axis
   *  - 2: allow rotation around the y axis
   *  - 4: allow rotation around the z axis
   *  By default, the particle cannot rotate.
   */
  uint8_t rotation = static_cast<uint8_t>(0b000u);
#else
  /** Bitfield for the particle axes of rotation. Particle cannot rotate. */
  static constexpr uint8_t rotation = static_cast<uint8_t>(0b000u);
#endif

#ifdef ESPRESSO_EXTERNAL_FORCES
  /** Flag for fixed particle coordinates.
   *  Values:
   *  - 0: no fixed coordinates
   *  - 1: fix translation along the x axis
   *  - 2: fix translation along the y axis
   *  - 4: fix translation along the z axis
   */
  uint8_t ext_flag = static_cast<uint8_t>(0b000u);
#else  // ESPRESSO_EXTERNAL_FORCES
  /** Bitfield for fixed particle coordinates. Coordinates cannot be fixed. */
  static constexpr uint8_t ext_flag = static_cast<uint8_t>(0b000u);
#endif // ESPRESSO_EXTERNAL_FORCES

  /** particle mass */
#ifdef ESPRESSO_MASS
  double mass = 1.0;
#else
  constexpr static double mass{1.0};
#endif

  /** rotational inertia */
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d rinertia = {1., 1., 1.};
#else
  static constexpr Utils::Vector3d rinertia = {1., 1., 1.};
#endif

  /** charge. */
#ifdef ESPRESSO_ELECTROSTATICS
  double q = 0.0;
#else
  constexpr static double q{0.0};
#endif

#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  Utils::Vector3d mu_E = {0., 0., 0.};
#endif

#ifdef ESPRESSO_DIPOLES
  /** dipole moment (absolute value) */
  double dipm = 0.;
#endif

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  /** total dipole field */
  Utils::Vector3d dip_fld = {0., 0., 0.};
#endif

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  /** The following properties define, with respect to which real particle a
   *  virtual site is placed and at what distance. The relative orientation of
   *  the vector pointing from real particle to virtual site with respect to the
   *  orientation of the real particle is stored in the virtual site's
   *  quaternion attribute.
   */
  struct VirtualSitesRelativeParameters {
    int to_particle_id = -1;
    double distance = 0.;
    /** Relative position of the virtual site. */
    Utils::Quaternion<double> rel_orientation =
        Utils::Quaternion<double>::identity();
    /** Orientation of the virtual particle in the body fixed frame. */
    Utils::Quaternion<double> quat = Utils::Quaternion<double>::identity();

    template <class Archive> void serialize(Archive &ar, long int) {
      ar & to_particle_id;
      ar & distance;
      ar & rel_orientation;
      ar & quat;
    }
  } vs_relative;
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE

#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
/** Friction coefficient for translation */
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double gamma = -1.;
#else
  Utils::Vector3d gamma = {-1., -1., -1.};
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#ifdef ESPRESSO_ROTATION
/** Friction coefficient for rotation */
#ifndef ESPRESSO_PARTICLE_ANISOTROPY
  double gamma_rot = -1.;
#else
  Utils::Vector3d gamma_rot = {-1., -1., -1.};
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE

#ifdef ESPRESSO_EXTERNAL_FORCES
  /** External force. */
  Utils::Vector3d ext_force = {0., 0., 0.};
#ifdef ESPRESSO_ROTATION
  /** External torque. */
  Utils::Vector3d ext_torque = {0., 0., 0.};
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_EXTERNAL_FORCES

#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming swim;
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  ThermalStonerWohlfarthParameters magnetodynamics;
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & identity;
    ar & mol_id;
    ar & type;
    ar & propagation;
#ifdef ESPRESSO_MASS
    ar & mass;
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    ar & rinertia;
#endif
#ifdef ESPRESSO_ROTATION
    ar & rotation;
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    ar & q;
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    ar & mu_E;
#endif
#ifdef ESPRESSO_DIPOLES
    ar & dipm;
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    ar & dip_fld;
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    ar & vs_relative;
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    ar & gamma;
#ifdef ESPRESSO_ROTATION
    ar & gamma_rot;
#endif
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & ext_flag;
    ar & ext_force;
#ifdef ESPRESSO_ROTATION
    ar & ext_torque;
#endif
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
    ar & swim;
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
    ar & magnetodynamics;
#endif
  }
};

/** Positional information on a particle. Information that is
 *  communicated to calculate interactions with ghost particles.
 */
struct ParticlePosition {
  /** periodically folded position. */
  Utils::Vector3d p = {0., 0., 0.};
  /** index of the simulation box image where the particle really sits. */
  Utils::Vector3i i = {0, 0, 0};

#ifdef ESPRESSO_ROTATION
  /** quaternion to define particle orientation */
  Utils::Quaternion<double> quat = Utils::Quaternion<double>::identity();
  /** unit director calculated from the quaternion */
  constexpr Utils::Vector3d calc_director() const {
    return Utils::convert_quaternion_to_director(quat);
  }
#endif

#ifdef ESPRESSO_BOND_CONSTRAINT
  /** particle position at the previous time step (RATTLE algorithm) */
  Utils::Vector3d p_last_timestep = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & p;
    ar & i;
#ifdef ESPRESSO_ROTATION
    ar & quat;
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    ar & p_last_timestep;
#endif
  }
};

/** Force information on a particle. Forces of ghost particles are
 *  collected and added up to the force of the original particle.
 */
struct ParticleForce {
  ParticleForce() = default;
  ParticleForce(ParticleForce const &) = default;
  ParticleForce &operator=(ParticleForce const &) = default;
  ParticleForce(const Utils::Vector3d &f) : f(f) {}
#ifdef ESPRESSO_ROTATION
  ParticleForce(const Utils::Vector3d &f, const Utils::Vector3d &torque)
      : f(f), torque(torque) {}
#endif

  friend ParticleForce operator+(ParticleForce const &lhs,
                                 ParticleForce const &rhs) {
    ParticleForce result = lhs;
    result += rhs;
    return result;
  }

  ParticleForce &operator+=(ParticleForce const &rhs) {
    f += rhs.f;
#ifdef ESPRESSO_ROTATION
    torque += rhs.torque;
#endif
    return *this;
  }

  /** force. */
  Utils::Vector3d f = {0., 0., 0.};

#ifdef ESPRESSO_ROTATION
  /** torque. */
  Utils::Vector3d torque = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & f;
#ifdef ESPRESSO_ROTATION
    ar & torque;
#endif
  }
};

/** Momentum information on a particle. Information not contained in
 *  communication of ghost particles so far, but a communication would
 *  be necessary for velocity-dependent potentials.
 */
struct ParticleMomentum {
  /** velocity. */
  Utils::Vector3d v = {0., 0., 0.};

#ifdef ESPRESSO_ROTATION
  /** angular velocity.
   *  ALWAYS IN PARTICLE FIXED, I.E., CO-ROTATING COORDINATE SYSTEM.
   */
  Utils::Vector3d omega = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & v;
#ifdef ESPRESSO_ROTATION
    ar & omega;
#endif
  }
};

/** Information on a particle that is needed only on the
 *  node the particle belongs to.
 */
struct ParticleLocal {
  /** is particle a ghost particle. */
  bool ghost = false;
  short int lees_edwards_flag = 0;
  /** position from the last Verlet list update. */
  Utils::Vector3d p_old = {0., 0., 0.};
  /** Accumulated applied Lees-Edwards offset. */
  double lees_edwards_offset = 0.;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & ghost;
    ar & lees_edwards_flag;
    ar & p_old;
    ar & lees_edwards_offset;
  }
};

#ifdef ESPRESSO_BOND_CONSTRAINT
struct ParticleRattle {
  /** position/velocity correction */
  Utils::Vector3d correction = {0., 0., 0.};

  friend ParticleRattle operator+(ParticleRattle const &lhs,
                                  ParticleRattle const &rhs) {
    return {lhs.correction + rhs.correction};
  }

  ParticleRattle &operator+=(ParticleRattle const &rhs) {
    return *this = *this + rhs;
  }

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & correction;
  }
};
#endif

/** Struct holding all information for one particle. */
struct Particle { // NOLINT(bugprone-exception-escape)
private:
  ParticleProperties p;
  ParticlePosition r;
  ParticleMomentum m;
  ParticleForce f;
  ParticleLocal l;
#ifdef ESPRESSO_BOND_CONSTRAINT
  ParticleRattle rattle;
#endif
  BondList bl;
#ifdef ESPRESSO_EXCLUSIONS
  /** list of particles, with which this particle has no non-bonded
   *  interactions
   */
  Utils::compact_vector<int> el;
#endif

public:
  constexpr auto const &id() const { return p.identity; }
  constexpr auto &id() { return p.identity; }
  constexpr auto const &mol_id() const { return p.mol_id; }
  constexpr auto &mol_id() { return p.mol_id; }
  constexpr auto const &type() const { return p.type; }
  constexpr auto &type() { return p.type; }

  constexpr auto const &propagation() const { return p.propagation; }
  constexpr auto &propagation() { return p.propagation; }

  constexpr bool operator==(Particle const &rhs) const {
    return id() == rhs.id();
  }

  constexpr bool operator!=(Particle const &rhs) const {
    return id() != rhs.id();
  }

  constexpr auto const &bonds() const { return bl; }
  constexpr auto &bonds() { return bl; }

  constexpr auto const &pos() const { return r.p; }
  constexpr auto &pos() { return r.p; }
  constexpr auto const &v() const { return m.v; }
  constexpr auto &v() { return m.v; }
  constexpr auto const &force() const { return f.f; }
  constexpr auto &force() { return f.f; }
  constexpr auto const &force_and_torque() const { return f; }
  constexpr auto &force_and_torque() { return f; }

  constexpr bool is_ghost() const { return l.ghost; }
  constexpr void set_ghost(bool const ghost_flag) { l.ghost = ghost_flag; }
  constexpr auto &pos_at_last_verlet_update() { return l.p_old; }
  constexpr auto const &pos_at_last_verlet_update() const { return l.p_old; }
  constexpr auto const &image_box() const { return r.i; }
  constexpr auto &image_box() { return r.i; }
  constexpr auto const &lees_edwards_offset() const {
    return l.lees_edwards_offset;
  }
  constexpr auto &lees_edwards_offset() { return l.lees_edwards_offset; }
  constexpr auto const &lees_edwards_flag() const {
    return l.lees_edwards_flag;
  }
  constexpr auto &lees_edwards_flag() { return l.lees_edwards_flag; }

  constexpr auto const &mass() const { return p.mass; }
#ifdef ESPRESSO_MASS
  constexpr auto &mass() { return p.mass; }
#endif
#ifdef ESPRESSO_ROTATION
  constexpr auto const &rotation() const { return p.rotation; }
  constexpr auto &rotation() { return p.rotation; }
  constexpr bool can_rotate() const { return static_cast<bool>(p.rotation); }
  constexpr bool can_rotate_around(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(p.rotation, axis);
  }
  constexpr void set_can_rotate_around(unsigned int const axis,
                                       bool const rot_flag) {
    assert(axis <= 2u);
    if (rot_flag) {
      p.rotation |= static_cast<uint8_t>(1u << axis);
    } else {
      p.rotation &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  constexpr void set_can_rotate_all_axes() {
    p.rotation = static_cast<uint8_t>(0b111u);
  }
  constexpr void set_cannot_rotate_all_axes() {
    p.rotation = static_cast<uint8_t>(0b000u);
  }
  constexpr auto const &quat() const { return r.quat; }
  constexpr auto &quat() { return r.quat; }
  constexpr auto const &torque() const { return f.torque; }
  constexpr auto &torque() { return f.torque; }
  constexpr auto const &omega() const { return m.omega; }
  constexpr auto &omega() { return m.omega; }
#ifdef ESPRESSO_EXTERNAL_FORCES
  constexpr auto const &ext_torque() const { return p.ext_torque; }
  constexpr auto &ext_torque() { return p.ext_torque; }
#endif // ESPRESSO_EXTERNAL_FORCES
  constexpr auto calc_director() const { return r.calc_director(); }
#else  // ESPRESSO_ROTATION
  constexpr auto can_rotate() const { return false; }
  constexpr auto can_rotate_around(unsigned int const) const { return false; }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_DIPOLES
  constexpr auto const &dipm() const { return p.dipm; }
  constexpr auto &dipm() { return p.dipm; }
  constexpr auto calc_dip() const { return calc_director() * dipm(); }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  constexpr auto const &stoner_wohlfarth_is_enabled() const {
    return p.magnetodynamics.is_enabled;
  }
  constexpr auto &stoner_wohlfarth_is_enabled() {
    return p.magnetodynamics.is_enabled;
  }
  constexpr auto const &stoner_wohlfarth_phi_0() const {
    return p.magnetodynamics.phi0;
  }
  constexpr auto &stoner_wohlfarth_phi_0() { return p.magnetodynamics.phi0; }
  constexpr auto const &saturation_magnetization() const {
    return p.magnetodynamics.sat_mag;
  }
  constexpr auto &saturation_magnetization() {
    return p.magnetodynamics.sat_mag;
  }
  constexpr auto const &magnetic_anisotropy_field_inv() const {
    return p.magnetodynamics.ani_fld_inv;
  }
  constexpr auto &magnetic_anisotropy_field_inv() {
    return p.magnetodynamics.ani_fld_inv;
  }
  constexpr auto const &magnetic_anisotropy_energy() const {
    return p.magnetodynamics.ani_energy;
  }
  constexpr auto &magnetic_anisotropy_energy() {
    return p.magnetodynamics.ani_energy;
  }
  constexpr auto const &stoner_wohlfarth_tau0_inv() const {
    return p.magnetodynamics.tau0_inv;
  }
  constexpr auto &stoner_wohlfarth_tau0_inv() {
    return p.magnetodynamics.tau0_inv;
  }
  constexpr auto const &stoner_wohlfarth_dt_incr() const {
    return p.magnetodynamics.dt_incr;
  }
  constexpr auto &stoner_wohlfarth_dt_incr() {
    return p.magnetodynamics.dt_incr;
  }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  constexpr auto const &dip_fld() const { return p.dip_fld; }
  constexpr auto &dip_fld() { return p.dip_fld; }
#endif
  constexpr auto const &rinertia() const { return p.rinertia; }
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  constexpr auto &rinertia() { return p.rinertia; }
#endif
  constexpr auto const &q() const { return p.q; }
#ifdef ESPRESSO_ELECTROSTATICS
  constexpr auto &q() { return p.q; }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  constexpr auto const &mu_E() const { return p.mu_E; }
  constexpr auto &mu_E() { return p.mu_E; }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES
  constexpr auto is_virtual() const {
    return (p.propagation & (PropagationMode::TRANS_VS_RELATIVE |
                             PropagationMode::TRANS_VS_CENTER_OF_MASS |
                             PropagationMode::ROT_VS_RELATIVE |
                             PropagationMode::ROT_VS_INDEPENDENT |
                             PropagationMode::TRANS_LB_TRACER)) != 0;
  }
#else
  constexpr auto is_virtual() const { return false; }
#endif // ESPRESSO_VIRTUAL_SITES
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  constexpr auto const &vs_relative() const { return p.vs_relative; }
  constexpr auto &vs_relative() { return p.vs_relative; }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  constexpr auto const &gamma() const { return p.gamma; }
  constexpr auto &gamma() { return p.gamma; }
#ifdef ESPRESSO_ROTATION
  constexpr auto const &gamma_rot() const { return p.gamma_rot; }
  constexpr auto &gamma_rot() { return p.gamma_rot; }
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
  constexpr auto const &fixed() const { return p.ext_flag; }
  constexpr auto &fixed() { return p.ext_flag; }
  constexpr bool has_fixed_coordinates() const {
    return static_cast<bool>(p.ext_flag);
  }
  constexpr bool is_fixed_along(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(p.ext_flag, axis);
  }
  constexpr void set_fixed_along(int const axis, bool const fixed_flag) {
    // set new flag
    if (fixed_flag) {
      p.ext_flag |= static_cast<uint8_t>(1u << axis);
    } else {
      p.ext_flag &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  constexpr auto const &ext_force() const { return p.ext_force; }
  constexpr auto &ext_force() { return p.ext_force; }
#else  // ESPRESSO_EXTERNAL_FORCES
  constexpr bool has_fixed_coordinates() const { return false; }
  constexpr bool is_fixed_along(unsigned int const) const { return false; }
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  constexpr auto const &swimming() const { return p.swim; }
  constexpr auto &swimming() { return p.swim; }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  constexpr auto const &pos_last_time_step() const { return r.p_last_timestep; }
  constexpr auto &pos_last_time_step() { return r.p_last_timestep; }
  constexpr auto const &rattle_params() const { return rattle; }
  constexpr auto &rattle_params() { return rattle; }
  constexpr auto const &rattle_correction() const { return rattle.correction; }
  constexpr auto &rattle_correction() { return rattle.correction; }
#endif

#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> &exclusions() { return el; }
  Utils::compact_vector<int> const &exclusions() const { return el; }
  bool has_exclusion(int pid) const {
    return std::ranges::find(el, pid) != el.end();
  }
#endif

private:
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & p;
    ar & r;
    ar & m;
    ar & f;
    ar & l;
    ar & bl;
#ifdef ESPRESSO_EXCLUSIONS
    ar & el;
#endif
  }
};

BOOST_CLASS_IMPLEMENTATION(Particle, object_serializable)
#ifdef ESPRESSO_ENGINE
BOOST_CLASS_IMPLEMENTATION(ParticleParametersSwimming, object_serializable)
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_CLASS_IMPLEMENTATION(ThermalStonerWohlfarthParameters,
                           object_serializable)
#endif
BOOST_CLASS_IMPLEMENTATION(ParticleProperties, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticlePosition, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleMomentum, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleForce, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleLocal, object_serializable)
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_CLASS_IMPLEMENTATION(ParticleRattle, object_serializable)
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
BOOST_CLASS_IMPLEMENTATION(decltype(ParticleProperties::vs_relative),
                           object_serializable)
#endif

#ifdef ESPRESSO_ENGINE
BOOST_IS_BITWISE_SERIALIZABLE(ParticleParametersSwimming)
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_IS_BITWISE_SERIALIZABLE(ThermalStonerWohlfarthParameters)
#endif
BOOST_IS_BITWISE_SERIALIZABLE(ParticleProperties)
BOOST_IS_BITWISE_SERIALIZABLE(ParticlePosition)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleMomentum)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleForce)
BOOST_IS_BITWISE_SERIALIZABLE(ParticleLocal)
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_IS_BITWISE_SERIALIZABLE(ParticleRattle)
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
BOOST_IS_BITWISE_SERIALIZABLE(decltype(ParticleProperties::vs_relative))
#endif
