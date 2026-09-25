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
#include "particle_store/ParticleParameters.hpp"
#include "particle_store/ParticleStore.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/quaternion.hpp>

#include <boost/serialization/is_bitwise_serializable.hpp>
#include <boost/serialization/level.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace detail {
constexpr inline bool get_nth_bit(uint8_t const bitfield,
                                  unsigned int const bit_idx) {
  return bitfield & (1u << bit_idx);
}
} // namespace detail

// ParticleParametersSwimming, ThermalStonerWohlfarthParameters and
// VirtualSitesRelativeParameters (the cold parameter PODs) are defined in
// particle_store/ParticleParameters.hpp so that both this header and
// ParticleStore can name the complete types.

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 *
 *  These parameter fields are not @ref Particle members; they live only in the
 *  @ref ParticleStore columns / host sidecars (single ownership).
 *  @ref ParticleProperties is retained purely as a type anchor: several call
 *  sites name its field types via @c decltype(ParticleProperties::identity)
 *  etc., the @c ParticleProperties::VirtualSitesRelativeParameters member
 *  typedef keeps the script-interface spelling compiling, and the standalone
 *  properties-serialization unit test exercises it. It is neither a member of
 *  @ref Particle nor part of @ref Particle serialization.
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
  /** @see ::VirtualSitesRelativeParameters (defined in
   *  particle_store/ParticleParameters.hpp). Kept as a member typedef so that
   *  the spelling @c ParticleProperties::VirtualSitesRelativeParameters used by
   *  the script interface keeps compiling. */
  using VirtualSitesRelativeParameters = ::VirtualSitesRelativeParameters;
  VirtualSitesRelativeParameters vs_relative;
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
 *
 *  The position, image box, quaternion and position-at-last-time-step fields
 *  live in the @ref ParticleStore columns (single ownership). This empty struct
 *  is retained as the ghost-transfer documentation anchor referenced in
 *  ghosts.hpp; it is neither a member of @ref Particle nor serialized.
 */
struct ParticlePosition {};

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

/** Momentum information on a particle. Information communicated to calculate
 *  velocity-dependent interactions with ghost particles.
 *
 *  The velocity and angular-velocity fields live in the @ref ParticleStore
 *  columns (single ownership). This empty struct is retained as the
 *  ghost-transfer documentation anchor referenced in ghosts.hpp; it is neither
 *  a member of @ref Particle nor serialized.
 */
struct ParticleMomentum {};

#ifdef ESPRESSO_BOND_CONSTRAINT
/** The position/velocity RATTLE correction lives in a @ref ParticleStore
 *  observable column (single ownership); it is neither a member of @ref
 *  Particle nor serialized. The RATTLE ghost wire archives one
 *  @ref Utils::Vector3d (24 B) directly via @ref Particle::rattle_correction().
 *
 *  This empty struct is retained purely as the doxygen anchor for
 *  @ref GHOSTTRANS_RATTLE in ghosts.hpp, parallel to @ref ParticleMomentum for
 *  GHOSTTRANS_MOMENTUM; removing it would dangle those cross-references.
 */
struct ParticleRattle {};
#endif

/** @brief Non-owning VIEW of one particle stored in a @ref ParticleStore.
 *
 *  @ref Particle is a two-word view {store pointer, store row}; ALL per-field
 *  data lives in the @ref ParticleStore columns / host sidecars (single
 *  ownership). A @c Particle never carries data: every accessor is an
 *  unconditional store read at @c m_store_row, guarded (in debug builds) by
 *  @c assert(m_particle_store != nullptr). Cross-rank migration and the
 *  head-node fetch cache transfer data per field (see
 *  particle_store/MigrationPack) sourced/sunk directly from store columns.
 *
 *  A default-constructed @c Particle is an UNBOUND view (null store, row -1);
 *  it must be bound via @ref attach_to_store (or produced by
 *  @ref ParticleStore::make_view) before any accessor is used. Copying a
 *  @c Particle copies the two handle words -- two views over the SAME row alias
 *  the same columns (they do not snapshot). The
 *  constexpr-when-a-feature-is-disabled accessors keep their static fallbacks
 *  (there is no column when the feature is off, so nothing to read).
 */
struct Particle { // NOLINT(bugprone-exception-escape)
private:
  /** The store this view reads/writes; null while unbound. Rank-local; never
   *  serialized. */
  ParticleStore *m_particle_store = nullptr;
  /** Row of this particle in @c m_particle_store; -1 while unbound. */
  int m_store_row = -1;

  /** Static fallbacks for the constexpr-when-disabled parameter accessors.
   *  When the feature is off there is no ParticleStore column; the accessor
   *  returns a reference to this immutable default so the read-only accessor
   *  keeps its constexpr semantics. */
#ifndef ESPRESSO_MASS
  static constexpr double mass_fallback{1.0};
#endif
#ifndef ESPRESSO_ROTATIONAL_INERTIA
  static constexpr Utils::Vector3d rinertia_fallback = {1., 1., 1.};
#endif
#ifndef ESPRESSO_ELECTROSTATICS
  static constexpr double q_fallback{0.0};
#endif

public:
  void attach_to_store(ParticleStore &store, int const row) {
    m_particle_store = &store;
    m_store_row = row;
  }
  auto store() const { return m_particle_store; }
  auto store_row() const { return m_store_row; }

  int const &id() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->id(m_store_row);
  }
  int &id() {
    assert(m_particle_store != nullptr);
    return m_particle_store->id(m_store_row);
  }
  int const &mol_id() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->mol_id(m_store_row);
  }
  int &mol_id() {
    assert(m_particle_store != nullptr);
    return m_particle_store->mol_id(m_store_row);
  }
  int const &type() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->type(m_store_row);
  }
  int &type() {
    assert(m_particle_store != nullptr);
    return m_particle_store->type(m_store_row);
  }

  int const &propagation() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->propagation(m_store_row);
  }
  int &propagation() {
    assert(m_particle_store != nullptr);
    return m_particle_store->propagation(m_store_row);
  }

  bool operator==(Particle const &rhs) const { return id() == rhs.id(); }

  bool operator!=(Particle const &rhs) const { return id() != rhs.id(); }

  BondList const &bonds() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->bonds_sidecar_reference(m_store_row);
  }
  BondList &bonds() {
    assert(m_particle_store != nullptr);
    return m_particle_store->bonds_sidecar_reference(m_store_row);
  }

  Utils::Vector3d pos() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_value(m_store_row);
  }
  VectorReference pos() {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_reference(m_store_row);
  }
  Utils::Vector3d v() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->velocity_value(m_store_row);
  }
  VectorReference v() {
    assert(m_particle_store != nullptr);
    return m_particle_store->velocity_reference(m_store_row);
  }
  auto force() {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_reference(m_store_row);
  }
  Utils::Vector3d force() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_value(m_store_row);
  }

  /** @brief Whether this is a ghost particle.
   *
   *  Ghost-ness is STRUCTURAL. The store lays out local rows @c [0, n_local)
   *  first and ghost rows @c [n_local, n_total) after, so a view is a ghost iff
   *  its row is in the ghost suffix. This is the sole source of truth; it
   *  requires an attached view (asserted). */
  bool is_ghost() const {
    assert(m_particle_store != nullptr);
    return static_cast<std::size_t>(m_store_row) >=
           m_particle_store->number_of_local_particles();
  }
  VectorReference pos_at_last_verlet_update() {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_at_last_verlet_update_reference(
        m_store_row);
  }
  Utils::Vector3d pos_at_last_verlet_update() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_at_last_verlet_update_value(m_store_row);
  }
  Utils::Vector3i image_box() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->image_box_value(m_store_row);
  }
  IntegerVectorReference image_box() {
    assert(m_particle_store != nullptr);
    return m_particle_store->image_box_reference(m_store_row);
  }
  double lees_edwards_offset() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->lees_edwards_offset(m_store_row);
  }
  double &lees_edwards_offset() {
    assert(m_particle_store != nullptr);
    return m_particle_store->lees_edwards_offset(m_store_row);
  }
  short lees_edwards_flag() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->lees_edwards_flag(m_store_row);
  }
  short &lees_edwards_flag() {
    assert(m_particle_store != nullptr);
    return m_particle_store->lees_edwards_flag(m_store_row);
  }

#ifdef ESPRESSO_MASS
  double const &mass() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->mass(m_store_row);
  }
  double &mass() {
    assert(m_particle_store != nullptr);
    return m_particle_store->mass(m_store_row);
  }
#else
  constexpr auto &mass() const { return mass_fallback; }
#endif
#ifdef ESPRESSO_ROTATION
  std::uint8_t const &rotation() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->rotation(m_store_row);
  }
  std::uint8_t &rotation() {
    assert(m_particle_store != nullptr);
    return m_particle_store->rotation(m_store_row);
  }
  bool can_rotate() const { return static_cast<bool>(rotation()); }
  bool can_rotate_around(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(rotation(), axis);
  }
  void set_can_rotate_around(unsigned int const axis, bool const rot_flag) {
    assert(axis <= 2u);
    if (rot_flag) {
      rotation() |= static_cast<uint8_t>(1u << axis);
    } else {
      rotation() &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  void set_can_rotate_all_axes() { rotation() = static_cast<uint8_t>(0b111u); }
  void set_cannot_rotate_all_axes() {
    rotation() = static_cast<uint8_t>(0b000u);
  }
  Utils::Quaternion<double> quat() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->quaternion_value(m_store_row);
  }
  QuaternionReference quat() {
    assert(m_particle_store != nullptr);
    return m_particle_store->quaternion_reference(m_store_row);
  }
  auto torque() {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_reference(m_store_row);
  }
  Utils::Vector3d torque() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_value(m_store_row);
  }
  Utils::Vector3d omega() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->angular_velocity_value(m_store_row);
  }
  VectorReference omega() {
    assert(m_particle_store != nullptr);
    return m_particle_store->angular_velocity_reference(m_store_row);
  }
#ifdef ESPRESSO_EXTERNAL_FORCES
  Utils::Vector3d ext_torque() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_torque_value(m_store_row);
  }
  VectorReference ext_torque() {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_torque_reference(m_store_row);
  }
#endif // ESPRESSO_EXTERNAL_FORCES
  auto calc_director() const {
    return Utils::convert_quaternion_to_director(quat());
  }
#else  // ESPRESSO_ROTATION
  constexpr auto can_rotate() const { return false; }
  constexpr auto can_rotate_around(unsigned int const) const { return false; }
#endif // ESPRESSO_ROTATION
#ifdef ESPRESSO_DIPOLES
  double const &dipm() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->dipm(m_store_row);
  }
  double &dipm() {
    assert(m_particle_store != nullptr);
    return m_particle_store->dipm(m_store_row);
  }
  auto calc_dip() const { return calc_director() * dipm(); }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  // The Stoner-Wohlfarth parameters live in the ParticleStore host sidecar (a
  // whole POD indexed by store row). The accessors return a reference into it;
  // the individual-field accessors then reference into that.
  ThermalStonerWohlfarthParameters &magnetodynamics() {
    assert(m_particle_store != nullptr);
    return m_particle_store->magnetodynamics(m_store_row);
  }
  ThermalStonerWohlfarthParameters const &magnetodynamics() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->magnetodynamics(m_store_row);
  }
  auto const &stoner_wohlfarth_is_enabled() const {
    return magnetodynamics().is_enabled;
  }
  auto &stoner_wohlfarth_is_enabled() { return magnetodynamics().is_enabled; }
  auto const &stoner_wohlfarth_phi_0() const { return magnetodynamics().phi0; }
  auto &stoner_wohlfarth_phi_0() { return magnetodynamics().phi0; }
  auto const &saturation_magnetization() const {
    return magnetodynamics().sat_mag;
  }
  auto &saturation_magnetization() { return magnetodynamics().sat_mag; }
  auto const &magnetic_anisotropy_field_inv() const {
    return magnetodynamics().ani_fld_inv;
  }
  auto &magnetic_anisotropy_field_inv() {
    return magnetodynamics().ani_fld_inv;
  }
  auto const &magnetic_anisotropy_energy() const {
    return magnetodynamics().ani_energy;
  }
  auto &magnetic_anisotropy_energy() { return magnetodynamics().ani_energy; }
  auto const &stoner_wohlfarth_tau0_inv() const {
    return magnetodynamics().tau0_inv;
  }
  auto &stoner_wohlfarth_tau0_inv() { return magnetodynamics().tau0_inv; }
  auto const &stoner_wohlfarth_dt_incr() const {
    return magnetodynamics().dt_incr;
  }
  auto &stoner_wohlfarth_dt_incr() { return magnetodynamics().dt_incr; }
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  Utils::Vector3d dip_fld() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->dip_fld_value(m_store_row);
  }
  VectorReference dip_fld() {
    assert(m_particle_store != nullptr);
    return m_particle_store->dip_fld_reference(m_store_row);
  }
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  Utils::Vector3d rinertia() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->rinertia_value(m_store_row);
  }
  VectorReference rinertia() {
    assert(m_particle_store != nullptr);
    return m_particle_store->rinertia_reference(m_store_row);
  }
#else
  constexpr auto &rinertia() const { return rinertia_fallback; }
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  double const &q() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->q(m_store_row);
  }
  double &q() {
    assert(m_particle_store != nullptr);
    return m_particle_store->q(m_store_row);
  }
#else
  constexpr auto &q() const { return q_fallback; }
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  Utils::Vector3d mu_E() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->mu_E_value(m_store_row);
  }
  VectorReference mu_E() {
    assert(m_particle_store != nullptr);
    return m_particle_store->mu_E_reference(m_store_row);
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES
  auto is_virtual() const {
    return (propagation() & (PropagationMode::TRANS_VS_RELATIVE |
                             PropagationMode::TRANS_VS_CENTER_OF_MASS |
                             PropagationMode::ROT_VS_RELATIVE |
                             PropagationMode::ROT_VS_INDEPENDENT |
                             PropagationMode::TRANS_LB_TRACER)) != 0;
  }
#else
  constexpr auto is_virtual() const { return false; }
#endif // ESPRESSO_VIRTUAL_SITES
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  VirtualSitesRelativeParameters &vs_relative() {
    assert(m_particle_store != nullptr);
    return m_particle_store->vs_relative(m_store_row);
  }
  VirtualSitesRelativeParameters const &vs_relative() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->vs_relative(m_store_row);
  }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  // gamma/gamma_rot: element reference (double&) when isotropic, or a
  // VectorReference when ESPRESSO_PARTICLE_ANISOTROPY selects per-axis
  // friction; const returns a value.
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d gamma() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_value(m_store_row);
  }
  VectorReference gamma() {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_reference(m_store_row);
  }
#else
  double gamma() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_value(m_store_row);
  }
  double &gamma() {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_reference(m_store_row);
  }
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#ifdef ESPRESSO_ROTATION
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  Utils::Vector3d gamma_rot() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_rot_value(m_store_row);
  }
  VectorReference gamma_rot() {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_rot_reference(m_store_row);
  }
#else
  double gamma_rot() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_rot_value(m_store_row);
  }
  double &gamma_rot() {
    assert(m_particle_store != nullptr);
    return m_particle_store->gamma_rot_reference(m_store_row);
  }
#endif // ESPRESSO_PARTICLE_ANISOTROPY
#endif // ESPRESSO_ROTATION
#endif // ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_EXTERNAL_FORCES
  std::uint8_t const &fixed() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_flag(m_store_row);
  }
  std::uint8_t &fixed() {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_flag(m_store_row);
  }
  bool has_fixed_coordinates() const { return static_cast<bool>(fixed()); }
  /** @brief Raw fixed-coordinate bitfield, by value.
   *  Hot integrator loops read this ONCE per particle and bit-test locally
   *  instead of calling @ref is_fixed_along (which re-reads the ParticleStore
   *  ext_flag column on every axis). Bit @c axis set means the coordinate is
   *  fixed. */
  std::uint8_t fixed_flags_byte() const { return fixed(); }
  bool is_fixed_along(unsigned int const axis) const {
    assert(axis <= 2u);
    return detail::get_nth_bit(fixed(), axis);
  }
  void set_fixed_along(int const axis, bool const fixed_flag) {
    // set new flag
    if (fixed_flag) {
      fixed() |= static_cast<uint8_t>(1u << axis);
    } else {
      fixed() &= static_cast<uint8_t>(~(1u << axis));
    }
  }
  Utils::Vector3d ext_force() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_force_value(m_store_row);
  }
  VectorReference ext_force() {
    assert(m_particle_store != nullptr);
    return m_particle_store->ext_force_reference(m_store_row);
  }
#else  // ESPRESSO_EXTERNAL_FORCES
  constexpr bool has_fixed_coordinates() const { return false; }
  constexpr std::uint8_t fixed_flags_byte() const {
    return static_cast<std::uint8_t>(0u);
  }
  constexpr bool is_fixed_along(unsigned int const) const { return false; }
#endif // ESPRESSO_EXTERNAL_FORCES
#ifdef ESPRESSO_ENGINE
  ParticleParametersSwimming &swimming() {
    assert(m_particle_store != nullptr);
    return m_particle_store->swimming(m_store_row);
  }
  ParticleParametersSwimming const &swimming() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->swimming(m_store_row);
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  Utils::Vector3d pos_last_time_step() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_last_time_step_value(m_store_row);
  }
  VectorReference pos_last_time_step() {
    assert(m_particle_store != nullptr);
    return m_particle_store->position_last_time_step_reference(m_store_row);
  }
  // RATTLE/SHAKE correction: a ParticleStore observable column (structurally
  // like force). It is per-iteration SHAKE scratch -- never persisted, never
  // migrated, never checkpointed -- so the accessor is attached-only (asserts a
  // store), exactly like force().
  VectorReference rattle_correction() {
    assert(m_particle_store != nullptr);
    return m_particle_store->rattle_correction_reference(m_store_row);
  }
  Utils::Vector3d rattle_correction() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->rattle_correction_value(m_store_row);
  }
#endif

#ifdef ESPRESSO_EXCLUSIONS
  Utils::compact_vector<int> &exclusions() {
    assert(m_particle_store != nullptr);
    return m_particle_store->exclusions_sidecar_reference(m_store_row);
  }
  Utils::compact_vector<int> const &exclusions() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->exclusions_sidecar_reference(m_store_row);
  }
  bool has_exclusion(int pid) const {
    return std::ranges::find(exclusions(), pid) != exclusions().end();
  }
#endif
};

// A Particle is a two-word view {ParticleStore* + int row}. sizeof therefore
// collapses to a pointer, an int, and the alignment padding between them
// (8 + 4 + 4 = 16 on LP64). The assertion documents and pins this ABI so a
// stray data member can never sneak back in.
static_assert(sizeof(Particle) <= 2u * sizeof(void *),
              "Particle must be a two-word view (store pointer + row)");

BOOST_CLASS_IMPLEMENTATION(ParticleProperties, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticlePosition, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleMomentum, object_serializable)
BOOST_CLASS_IMPLEMENTATION(ParticleForce, object_serializable)
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_CLASS_IMPLEMENTATION(ParticleRattle, object_serializable)
#endif
#ifdef ESPRESSO_ENGINE
BOOST_CLASS_IMPLEMENTATION(ParticleParametersSwimming, object_serializable)
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
BOOST_CLASS_IMPLEMENTATION(ThermalStonerWohlfarthParameters,
                           object_serializable)
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
#ifdef ESPRESSO_BOND_CONSTRAINT
BOOST_IS_BITWISE_SERIALIZABLE(ParticleRattle)
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
BOOST_IS_BITWISE_SERIALIZABLE(decltype(ParticleProperties::vs_relative))
#endif
