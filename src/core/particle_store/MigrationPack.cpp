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

#include "particle_store/MigrationPack.hpp"

#include "BondList.hpp"
#include "particle_store/ParticleParameters.hpp"

#include <utils/Vector.hpp>
#include <utils/compact_vector.hpp>
#include <utils/quaternion.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>

// Per-field migration pack. The covered field set is the canonical
// ParticleStore field list (see the sync note above ParticleStore::assign_row),
// minus the RATTLE correction (per-iteration scratch, recomputed
// post-migration) and minus the ghost flag (migrating particles are never
// ghosts). Any field added to assign_row/copy_row must be added to BOTH the
// writer (append_row_*) and the reader (read_row_*) below, and to
// fixed_size_per_row(); the maximal-population round-trip unit test enforces
// this.

namespace MigrationPack {
namespace {

// -- native-endian POD memcpy cursors (mirror the ghost wire practice) --------

class WriteCursor {
  char *m_at;

public:
  explicit WriteCursor(std::vector<char> &buffer) : m_at(buffer.data()) {}
  template <class T> void put(T const &value) {
    std::memcpy(m_at, &value, sizeof(T));
    m_at += sizeof(T);
  }
  void put_bytes(void const *src, std::size_t bytes) {
    // A zero-length ragged run (e.g. a particle with no bonds/exclusions) has a
    // null data() pointer; memcpy declares both operands nonnull, so calling it
    // with (dst, nullptr, 0) is UB even though it copies nothing (trips UBSan).
    if (bytes != 0u) {
      std::memcpy(m_at, src, bytes);
      m_at += bytes;
    }
  }
  char const *position() const { return m_at; }
};

class ReadCursor {
  char const *m_at;

public:
  explicit ReadCursor(std::vector<char> const &buffer) : m_at(buffer.data()) {}
  template <class T> T get() {
    T value;
    std::memcpy(&value, m_at, sizeof(T));
    m_at += sizeof(T);
    return value;
  }
  void get_bytes(void *dst, std::size_t bytes) {
    // See WriteCursor::put_bytes: a zero-length run has a null data() pointer,
    // and memcpy(dst, src, 0) with a null operand is UB (trips UBSan).
    if (bytes != 0u) {
      std::memcpy(dst, m_at, bytes);
      m_at += bytes;
    }
  }
  char const *position() const { return m_at; }
};

// The row-count header and every ragged run length travel as this type.
using LengthType = std::uint64_t;

// -- per-field helpers over the fixed-width vector/scalar types ---------------

void put_vector(WriteCursor &cursor, Utils::Vector3d const &v) {
  cursor.put(v[0]);
  cursor.put(v[1]);
  cursor.put(v[2]);
}
Utils::Vector3d get_vector(ReadCursor &cursor) {
  auto const x = cursor.get<double>();
  auto const y = cursor.get<double>();
  auto const z = cursor.get<double>();
  return {x, y, z};
}
void put_ivector(WriteCursor &cursor, Utils::Vector3i const &v) {
  cursor.put(v[0]);
  cursor.put(v[1]);
  cursor.put(v[2]);
}
Utils::Vector3i get_ivector(ReadCursor &cursor) {
  auto const x = cursor.get<int>();
  auto const y = cursor.get<int>();
  auto const z = cursor.get<int>();
  return {x, y, z};
}
#ifdef ESPRESSO_ROTATION
void put_quaternion(WriteCursor &cursor, Utils::Quaternion<double> const &q) {
  cursor.put(q[0]);
  cursor.put(q[1]);
  cursor.put(q[2]);
  cursor.put(q[3]);
}
Utils::Quaternion<double> get_quaternion(ReadCursor &cursor) {
  Utils::Quaternion<double> q;
  q[0] = cursor.get<double>();
  q[1] = cursor.get<double>();
  q[2] = cursor.get<double>();
  q[3] = cursor.get<double>();
  return q;
}
#endif
#if defined(ESPRESSO_THERMOSTAT_PER_PARTICLE) &&                               \
    defined(ESPRESSO_PARTICLE_ANISOTROPY)
void put_gamma(WriteCursor &cursor, Utils::Vector3d const &g) {
  put_vector(cursor, g);
}
#elif defined(ESPRESSO_THERMOSTAT_PER_PARTICLE)
void put_gamma(WriteCursor &cursor, double const g) { cursor.put(g); }
#endif

// -- fixed-width per-row size (constant for a given ifdef config) -------------

std::size_t fixed_size_per_row() {
  std::size_t size = 0u;
  constexpr auto d = sizeof(double);
  // POSITION leg.
  size += 3u * d;           // position
  size += 3u * sizeof(int); // image_box
#ifdef ESPRESSO_ROTATION
  size += 4u * d; // quaternion
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  size += 3u * d; // position_last_time_step
#endif
  size += 3u * d;        // position_at_last_verlet_update
  size += d;             // lees_edwards_offset
  size += sizeof(short); // lees_edwards_flag
  // MOMENTUM leg.
  size += 3u * d; // velocity
#ifdef ESPRESSO_ROTATION
  size += 3u * d; // angular_velocity
#endif
  // PROPRTS leg.
  size += 4u * sizeof(int); // id, mol_id, type, propagation
#ifdef ESPRESSO_ROTATION
  size += sizeof(std::uint8_t); // rotation
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  size += sizeof(std::uint8_t); // ext_flag
#endif
#ifdef ESPRESSO_MASS
  size += d; // mass
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  size += d; // q
#endif
#ifdef ESPRESSO_DIPOLES
  size += d; // dipm
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  size += 3u * d; // rinertia
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  size += 3u * d; // mu_E
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  size += 3u * d; // dip_fld
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  size += 3u * d; // ext_force
#ifdef ESPRESSO_ROTATION
  size += 3u * d; // ext_torque
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  size += 3u * d; // gamma
#ifdef ESPRESSO_ROTATION
  size += 3u * d; // gamma_rot
#endif
#else
  size += d; // gamma
#ifdef ESPRESSO_ROTATION
  size += d; // gamma_rot
#endif
#endif
#endif
#ifdef ESPRESSO_ENGINE
  size += sizeof(ParticleParametersSwimming);
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  size += sizeof(ThermalStonerWohlfarthParameters);
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  size += sizeof(VirtualSitesRelativeParameters);
#endif
  // FORCE leg.
  size += 3u * d; // force
#ifdef ESPRESSO_ROTATION
  size += 3u * d; // torque
#endif
  return size;
}

// -- fixed-width legs, grouped like the GHOSTTRANS groups ---------------------

void append_row_fixed(WriteCursor &cursor, ParticleStore const &store,
                      int const row) {
  // POSITION leg (+ the three migration-only fields).
  put_vector(cursor, store.position_value(row));
  put_ivector(cursor, store.image_box_value(row));
#ifdef ESPRESSO_ROTATION
  put_quaternion(cursor, store.quaternion_value(row));
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  put_vector(cursor, store.position_last_time_step_value(row));
#endif
  put_vector(cursor, store.position_at_last_verlet_update_value(row));
  cursor.put(store.lees_edwards_offset(row));
  cursor.put(store.lees_edwards_flag(row));

  // MOMENTUM leg.
  put_vector(cursor, store.velocity_value(row));
#ifdef ESPRESSO_ROTATION
  put_vector(cursor, store.angular_velocity_value(row));
#endif

  // PROPRTS leg.
  cursor.put(store.id(row));
  cursor.put(store.mol_id(row));
  cursor.put(store.type(row));
  cursor.put(store.propagation(row));
#ifdef ESPRESSO_ROTATION
  cursor.put(store.rotation(row));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  cursor.put(store.ext_flag(row));
#endif
#ifdef ESPRESSO_MASS
  cursor.put(store.mass(row));
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  cursor.put(store.q(row));
#endif
#ifdef ESPRESSO_DIPOLES
  cursor.put(store.dipm(row));
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  put_vector(cursor, store.rinertia_value(row));
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  put_vector(cursor, store.mu_E_value(row));
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  put_vector(cursor, store.dip_fld_value(row));
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  put_vector(cursor, store.ext_force_value(row));
#ifdef ESPRESSO_ROTATION
  put_vector(cursor, store.ext_torque_value(row));
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
  put_gamma(cursor, store.gamma_value(row));
#ifdef ESPRESSO_ROTATION
  put_gamma(cursor, store.gamma_rot_value(row));
#endif
#endif
#ifdef ESPRESSO_ENGINE
  {
    auto const swim = store.swimming(row);
    cursor.put_bytes(&swim, sizeof(swim));
  }
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  {
    auto const md = store.magnetodynamics(row);
    cursor.put_bytes(&md, sizeof(md));
  }
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  {
    auto const vsr = store.vs_relative(row);
    cursor.put_bytes(&vsr, sizeof(vsr));
  }
#endif

  // FORCE leg.
  put_vector(cursor, store.force_value(row));
#ifdef ESPRESSO_ROTATION
  put_vector(cursor, store.torque_value(row));
#endif
}

void read_row_fixed(ReadCursor &cursor, ParticleStore &store, int const row) {
  // POSITION leg.
  store.position_reference(row) = get_vector(cursor);
  store.image_box_reference(row) = get_ivector(cursor);
#ifdef ESPRESSO_ROTATION
  store.quaternion_reference(row) = get_quaternion(cursor);
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  store.position_last_time_step_reference(row) = get_vector(cursor);
#endif
  store.position_at_last_verlet_update_reference(row) = get_vector(cursor);
  store.lees_edwards_offset(row) = cursor.get<double>();
  store.lees_edwards_flag(row) = cursor.get<short>();

  // MOMENTUM leg.
  store.velocity_reference(row) = get_vector(cursor);
#ifdef ESPRESSO_ROTATION
  store.angular_velocity_reference(row) = get_vector(cursor);
#endif

  // PROPRTS leg.
  store.id(row) = cursor.get<int>();
  store.mol_id(row) = cursor.get<int>();
  store.type(row) = cursor.get<int>();
  store.propagation(row) = cursor.get<int>();
#ifdef ESPRESSO_ROTATION
  store.rotation(row) = cursor.get<std::uint8_t>();
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_flag(row) = cursor.get<std::uint8_t>();
#endif
#ifdef ESPRESSO_MASS
  store.mass(row) = cursor.get<double>();
#endif
#ifdef ESPRESSO_ELECTROSTATICS
  store.q(row) = cursor.get<double>();
#endif
#ifdef ESPRESSO_DIPOLES
  store.dipm(row) = cursor.get<double>();
#endif
#ifdef ESPRESSO_ROTATIONAL_INERTIA
  store.rinertia_reference(row) = get_vector(cursor);
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
  store.mu_E_reference(row) = get_vector(cursor);
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  store.dip_fld_reference(row) = get_vector(cursor);
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
  store.ext_force_reference(row) = get_vector(cursor);
#ifdef ESPRESSO_ROTATION
  store.ext_torque_reference(row) = get_vector(cursor);
#endif
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
#ifdef ESPRESSO_PARTICLE_ANISOTROPY
  store.gamma_reference(row) = get_vector(cursor);
#ifdef ESPRESSO_ROTATION
  store.gamma_rot_reference(row) = get_vector(cursor);
#endif
#else
  store.gamma_reference(row) = cursor.get<double>();
#ifdef ESPRESSO_ROTATION
  store.gamma_rot_reference(row) = cursor.get<double>();
#endif
#endif
#endif
#ifdef ESPRESSO_ENGINE
  cursor.get_bytes(&store.swimming(row), sizeof(ParticleParametersSwimming));
#endif
#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
  cursor.get_bytes(&store.magnetodynamics(row),
                   sizeof(ThermalStonerWohlfarthParameters));
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  cursor.get_bytes(&store.vs_relative(row),
                   sizeof(VirtualSitesRelativeParameters));
#endif

  // FORCE leg.
  store.force_reference(row) = get_vector(cursor);
#ifdef ESPRESSO_ROTATION
  store.torque_reference(row) = get_vector(cursor);
#endif
}

// -- ragged legs (length-prefixed, like the ghost bond buffer) ----------------

std::size_t bonds_run_size(ParticleStore const &store, int const row) {
  return sizeof(LengthType) +
         store.bonds_sidecar_reference(row).storage().size() * sizeof(int);
}
void append_row_bonds(WriteCursor &cursor, ParticleStore const &store,
                      int const row) {
  auto const &run = store.bonds_sidecar_reference(row).storage();
  cursor.put(static_cast<LengthType>(run.size()));
  cursor.put_bytes(run.data(), run.size() * sizeof(int));
}
void read_row_bonds(ReadCursor &cursor, ParticleStore &store, int const row) {
  auto const n = static_cast<std::size_t>(cursor.get<LengthType>());
  BondList::storage_type run;
  run.resize(static_cast<BondList::storage_type::size_type>(n));
  cursor.get_bytes(run.data(), n * sizeof(int));
  store.bonds_sidecar_reference(row).replace_storage(std::move(run));
}

#ifdef ESPRESSO_EXCLUSIONS
std::size_t exclusions_run_size(ParticleStore const &store, int const row) {
  return sizeof(LengthType) +
         store.exclusions_sidecar_reference(row).size() * sizeof(int);
}
void append_row_exclusions(WriteCursor &cursor, ParticleStore const &store,
                           int const row) {
  auto const &run = store.exclusions_sidecar_reference(row);
  cursor.put(static_cast<LengthType>(run.size()));
  cursor.put_bytes(run.data(), run.size() * sizeof(int));
}
void read_row_exclusions(ReadCursor &cursor, ParticleStore &store,
                         int const row) {
  auto const n = static_cast<std::size_t>(cursor.get<LengthType>());
  auto &run = store.exclusions_sidecar_reference(row);
  run.resize(static_cast<Utils::compact_vector<int>::size_type>(n));
  cursor.get_bytes(run.data(), n * sizeof(int));
}
#endif

} // namespace

std::size_t packed_size(ParticleStore const &store, std::span<int const> rows) {
  // Fixed part: the row-count header (u64) + the constant-per-row fixed legs.
  // The header carries no per-row id list (the id is re-read from the PROPRTS
  // leg on unpack), so the header is just the count.
  std::size_t size = sizeof(LengthType) + rows.size() * fixed_size_per_row();
  // Ragged actuals.
  for (auto const row : rows) {
    size += bonds_run_size(store, row);
  }
#ifdef ESPRESSO_EXCLUSIONS
  for (auto const row : rows) {
    size += exclusions_run_size(store, row);
  }
#endif
  return size;
}

void pack_rows(ParticleStore const &store, std::span<int const> rows,
               std::vector<char> &buffer) {
  auto const total = packed_size(store, rows);
  buffer.resize(total);
  WriteCursor cursor{buffer};

  // Header: row count only (the per-row id is carried in the PROPRTS leg and
  // re-read on unpack, so no separate id list is packed).
  cursor.put(static_cast<LengthType>(rows.size()));

  // Fixed-width legs, one contiguous block per row.
  for (auto const row : rows) {
    append_row_fixed(cursor, store, row);
  }
  // Ragged bond leg, then ragged exclusion leg.
  for (auto const row : rows) {
    append_row_bonds(cursor, store, row);
  }
#ifdef ESPRESSO_EXCLUSIONS
  for (auto const row : rows) {
    append_row_exclusions(cursor, store, row);
  }
#endif

  assert(static_cast<std::size_t>(cursor.position() - buffer.data()) == total);
  static_cast<void>(total);
}

std::size_t unpack_rows(ParticleStore &store, int const first_row,
                        std::vector<char> const &buffer) {
  ReadCursor cursor{buffer};

  // Header: row count. The receiver peeks this same u64 to reserve the target
  // rows before unpacking (see RegularDecomposition::exchange_neighbors); the
  // per-row id is re-read from the PROPRTS leg below, so no separate id list
  // travels in the header.
  auto const count = static_cast<std::size_t>(cursor.get<LengthType>());

  for (std::size_t i = 0u; i < count; ++i) {
    read_row_fixed(cursor, store, first_row + static_cast<int>(i));
  }
  for (std::size_t i = 0u; i < count; ++i) {
    read_row_bonds(cursor, store, first_row + static_cast<int>(i));
  }
#ifdef ESPRESSO_EXCLUSIONS
  for (std::size_t i = 0u; i < count; ++i) {
    read_row_exclusions(cursor, store, first_row + static_cast<int>(i));
  }
#endif

  assert(static_cast<std::size_t>(cursor.position() - buffer.data()) ==
         buffer.size());
  return count;
}

} // namespace MigrationPack
