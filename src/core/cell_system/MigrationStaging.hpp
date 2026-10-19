/*
 * Copyright (C) 2026 The ESPResSo project
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
 * @brief Handle a @ref ParticleDecomposition uses to route migrating particles
 * through the owning @ref CellStructure's staging @ref ParticleStore.
 *
 * A mis-celled non-local particle is COPIED (column-by-column) out of the live
 * store into a row of the @ref CellStructure staging store (@c stage_row); the
 * migration `displaced_parts` list holds staging-row INDICES; per direction
 * those rows are @ref MigrationPack "packed" into a flat byte buffer; received
 * buffers are unpacked into fresh reserved staging rows (@c reserve_rows); and
 * a received row that belongs on this node is staged into its home cell as a
 * reference to that staging row (@ref CellParticleStorage::insert_staged_row)
 * -- the next store rebuild copies it into a committed row, so the final store
 * row-assignment order matches the cell staging order.
 *
 * @ref CellStructure owns the staging store and installs this handle on the
 * decomposition (like @c set_commit_store); the function objects wrap the
 * @ref CellStructure staging helpers. @c store points at the address-stable
 * staging-store MEMBER of @ref CellStructure (its internals may be swapped out
 * when it grows, but the member address does not change), so packing reads the
 * current columns through it AND a staged row reference into a cell always
 * names the current staging-store address.
 */

#include <functional>

class ParticleStore;

struct MigrationStaging {
  /** The staging store to pack from / unpack into (address-stable). It is also
   *  the source store a received/moved staging row is staged into a cell FROM
   *  (@ref CellParticleStorage::insert_staged_row). */
  ParticleStore *store = nullptr;
  /** Copy a live store row into a fresh staging row; returns the staging row.
   */
  std::function<int(int live_row)> stage_row;
  /** Reserve @p count fresh staging rows; returns the first reserved index. */
  std::function<int(int count)> reserve_rows;
  /** Drop all staged rows (reset the row counter to zero). */
  std::function<void()> clear;

  explicit operator bool() const { return store != nullptr; }
};
