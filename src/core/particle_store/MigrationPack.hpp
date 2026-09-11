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

/**
 * @file
 * @brief Per-field particle migration pack.
 *
 * Serializes a set of @ref ParticleStore rows into a flat byte buffer and back,
 * field-by-field, WITHOUT constructing a whole @ref Particle. This is the wire
 * for cross-rank particle migration (the global-resort migration and the
 * head-node fetch cache).
 *
 * Wire layout (native-endian POD @c memcpy, matching the ghost wire practice):
 *
 *   [ header : row-count (u64) ]   (the id is not listed per row in the header;
 *                                   it travels in the PROPRTS leg and the
 *                                   receiver peeks this count to reserve rows)
 *   [ fixed-width legs, grouped like the GHOSTTRANS groups, one contiguous
 *     block per group across all rows:
 *        POSITION  (pos, image_box, [quaternion], [pos_last_time_step],
 *                   pos_at_last_verlet_update, lees_edwards_offset,
 *                   lees_edwards_flag)   -- extends the ghost POSITION group
 *                   with the three migration-only fields
 *        MOMENTUM  (velocity, [angular_velocity])
 *        PROPRTS   (id, mol_id, type, propagation, [rotation], [ext_flag],
 *                   [mass], [q], [dipm], [rinertia], [mu_E], [dip_fld],
 *                   [ext_force], [ext_torque], [gamma], [gamma_rot],
 *                   [swimming], [magnetodynamics], [vs_relative])
 *        FORCE     (force, [torque]) ]
 *   [ ragged BONDS leg      : per row, run-length (u64) then that many ints ]
 *   [ ragged EXCLUSIONS leg : per row, run-length (u64) then that many ints ]
 *
 * The ghost flag (@c Particle::l.ghost) is deliberately OMITTED: a migrating
 * particle is never a ghost. @ref MigrationPack::unpack_rows asserts this
 * precondition holds.
 *
 * The RATTLE correction is omitted (as in @c Particle::serialize): it is a
 * per-iteration scratch recomputed on the owning rank post-migration.
 *
 * The covered field set is IDENTICAL to @ref ParticleStore::assign_row /
 * @ref ParticleStore::copy_row (see the sync note in ParticleStore.cpp), minus
 * the RATTLE correction. The maximal-population round-trip unit test enforces
 * completeness.
 */

#include "particle_store/ParticleStore.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace MigrationPack {

/**
 * @brief Exact byte size @ref pack_rows will produce for @p rows of @p store.
 *
 * A constant fixed part per row (the enabled fixed-width legs, determined at
 * compile time by the ifdef config) plus the row-count header (u64) and the
 * ragged bond/exclusion actuals summed over @p rows.
 */
std::size_t packed_size(ParticleStore const &store, std::span<int const> rows);

/**
 * @brief Pack @p rows of @p store into @p buffer (native-endian POD memcpy).
 *
 * @p buffer is resized to @ref packed_size and fully overwritten.
 */
void pack_rows(ParticleStore const &store, std::span<int const> rows,
               std::vector<char> &buffer);

/**
 * @brief Unpack a buffer produced by @ref pack_rows into consecutive store
 * rows.
 *
 * Writes the unpacked rows into @p store starting at @p first_row (so
 * @c first_row .. first_row + count - 1 must be valid, allocated rows -- the
 * caller sizes the store). Returns the number of rows unpacked (the row-count
 * header).
 */
std::size_t unpack_rows(ParticleStore &store, int first_row,
                        std::vector<char> const &buffer);

} // namespace MigrationPack
