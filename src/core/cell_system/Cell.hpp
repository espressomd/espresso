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

#include "Particle.hpp"
#include "cell_system/CellRows.hpp"
#include "cell_system/RowParticleRange.hpp"
#include "particle_store/ParticleStore.hpp"

#include <boost/range/iterator_range.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

/** @brief A particle staged for insertion into a cell.
 *
 *  A cell stages a reference to a row in a SOURCE store from which the next
 *  store rebuild copies every field into a fresh row via
 *  @ref ParticleStore::copy_row. Two shapes:
 *  - @c source_store != nullptr : copy from @c source_store at @c source_row
 *    (a migrating particle held in the CellStructure staging store, or a
 *    freshly-built new particle in a creation staging store).
 *  - @c source_store == nullptr : a fresh, default-seeded particle (a new ghost
 *    from @ref CellParticleStorage::resize_ghost_storage); the rebuild seeds
 *    the row to defaults (@ref ParticleStore::seed_default_row).
 *  Ghost-ness is STRUCTURAL (the row's position in the store), so no flag is
 *  carried here: ghost rows always land in the ghost suffix by construction. */
struct StagedParticle {
  ParticleStore *source_store = nullptr;
  int source_row = -1;
};

template <class CellRef> class Neighbors {
  using storage_type = std::vector<CellRef>;

public:
  using value_type = storage_type::value_type;
  using iterator = storage_type::iterator;
  using const_iterator = storage_type::const_iterator;
  using cell_range = boost::iterator_range<iterator>;

private:
  void copy(const Neighbors &rhs) {
    m_neighbors = rhs.m_neighbors;
    m_red_black_divider =
        m_neighbors.begin() +
        std::distance(rhs.m_neighbors.begin(),
                      const_iterator(rhs.m_red_black_divider));
  }

public:
  Neighbors() = default;
  Neighbors(const Neighbors &rhs) { copy(rhs); }
  Neighbors &operator=(const Neighbors &rhs) {
    copy(rhs);
    return *this;
  }

  Neighbors(std::span<const CellRef> red_neighbors,
            std::span<const CellRef> black_neighbors) {
    m_neighbors.resize(red_neighbors.size() + black_neighbors.size());
    auto const res = std::ranges::copy(red_neighbors, m_neighbors.begin());
    m_red_black_divider = res.out;
    std::ranges::copy(black_neighbors, m_red_black_divider);
  }

  /**
   * @brief All neighbors.
   */
  cell_range all() { return {m_neighbors.begin(), m_neighbors.end()}; }
  /**
   * @brief Red partition of neighbors.
   *
   * An partition of the neighbors so that iterating over all
   * neighbors_red of all cells visits every pair exactly once.
   * Complement of neighbors_black.
   */
  cell_range red() { return {m_neighbors.begin(), m_red_black_divider}; }
  /**
   * @brief Black partition of neighbors.
   *
   * An partition of the neighbors so that iterating over all
   * neighbors_black of all cells visits every pair exactly once.
   * Complement of neighbors_red.
   */
  cell_range black() { return {m_red_black_divider, m_neighbors.end()}; }

private:
  /** Container with all the neighbors.
      Red neighbors are first, black second. */
  storage_type m_neighbors;
  /** Iterator pointing to the first black neighbor
      in the container. */
  iterator m_red_black_divider;
};

class Cell {
  using neighbors_type = Neighbors<Cell *>;

  /** Contiguous store-row range of this cell's committed particles.
   *  The cell's committed particles occupy rows @c [m_offset, m_offset+m_count)
   *  of the @ref ParticleStore, in cell-traversal order; the range is
   *  (re)established by the permutation rebuild
   *  (@ref CellStructure::ensure_particle_store_synchronized). Cells are laid
   *  out back-to-back in the store by the counting sort, so a cell's rows are a
   *  pure arithmetic range and @ref rows() hands out a @ref CellRowSpan over
   * it. Together with @c m_staged this is the cell's authoritative particle
   *  content: the cell owns a ROW RANGE into the store plus a small buffer of
   *  not-yet-committed row references. A row dropped mid-epoch is marked
   *  pending-removed on the store (not physically removed from the range); the
   *  live @ref CellRowSpan / @ref RowParticleRange skip it. */
  std::size_t m_offset = 0u;
  std::size_t m_count = 0u;

  /** Not-yet-committed staged particles. A cell mutation that ADDS a particle
   *  (@ref CellParticleStorage::insert_staged_row, the migration re-insert
   *  path, ghost-cell resize) appends a @ref StagedParticle here (a reference
   *  to a source-store row, or a fresh-default marker) and marks the store
   *  dirty; the next store rebuild copies a fresh row per staged entry (via
   *  @ref ParticleStore::copy_row / @ref ParticleStore::seed_default_row) and
   *  clears this buffer. Between the mutation and the rebuild the staged
   *  particles are part of the cell but have no committed row yet, so
   *  @ref particles() cannot yet see them -- callers that add particles must
   *  trigger the rebuild (mark-dirty + ensure_particle_store_synchronized)
   *  before iterating. */
  std::vector<StagedParticle> m_staged;

  /** The store the row indices point into. Wired by @ref CellStructure at
   *  construction / decomposition swap; used by @ref particles() to hand out
   *  views. */
  ParticleStore *m_store = nullptr;

public:
  /** @brief Bind the @ref ParticleStore this cell's rows index into. */
  void set_store(ParticleStore &store) { m_store = &store; }
  /** @brief The store this cell's rows index into. */
  ParticleStore &store() {
    assert(m_store != nullptr);
    return *m_store;
  }

  /** @brief Particles in this cell as a range of @ref Particle VIEWS over the
   *  committed store rows. NOT including staged (uncommitted) particles --
   * those become visible only after the next store rebuild; pending-removed
   * rows are skipped. The returned range's views alias the store and are
   * invalidated by a rebuild (see @ref RowParticleRange). */
  RowParticleRange particles() {
    assert(m_store != nullptr);
    return {m_offset, m_count, *m_store};
  }
  RowParticleRange particles() const {
    assert(m_store != nullptr);
    return {m_offset, m_count, *m_store};
  }

  /** @brief The cell's committed store-row range as a @ref CellRowSpan. The
   *  live surface (@c begin/end/size/empty) skips pending-removed rows;
   *  @c offset()/count() expose the raw range. */
  CellRowSpan rows() const { return {m_offset, m_count, m_store}; }

  /** @brief Raw store-row range of this cell's committed particles.
   *  Removal-agnostic; the resort loops iterate these directly and mark dropped
   *  rows pending-removed. */
  std::size_t offset() const { return m_offset; }
  std::size_t count() const { return m_count; }
  /** @brief Set the committed row range (written by the store rebuild). */
  void set_range(std::size_t offset, std::size_t count) {
    m_offset = offset;
    m_count = count;
  }

  /** @brief The staging area of not-yet-committed detached particles. */
  auto &staged() { return m_staged; }
  auto const &staged() const { return m_staged; }

  /** @brief Total live particle count = live committed rows + staged. Committed
   *  rows marked pending-removed are excluded (@ref rows().size() is the live
   *  count). Ghost-count (PARTNUM) communication uses this so that a cell
   *  resized as a ghost destination (which stages its ghosts, deferring the row
   *  commit) reports its new size immediately when read as a source for a
   *  downstream ghost layer within the same communicator pass. */
  std::size_t size() const { return rows().size() + m_staged.size(); }

  neighbors_type m_neighbors;

  /** Interaction pairs as @ref ParticleStore ROW indices. Pairs are recorded as
   *  store rows rather than pointers, mirroring @c
   * CellStructure::m_verlet_list. NOTE: this member has no reader or writer
   * anywhere in the code base (the live Verlet list is @c
   * CellStructure::m_verlet_list); it is a candidate for removal. */
  std::vector<std::pair<int, int>> m_verlet_list;

  /**
   * @brief All neighbors of the cell.
   */
  neighbors_type &neighbors() { return m_neighbors; }

  /**
   * @brief Interior/boundary classification.
   *
   * A local cell is *boundary* iff at least one of its neighbors is a ghost
   * cell.  Set by GhostComm::mark_boundary_cells() after neighbor setup.
   * Default false (interior).
   */
  bool m_is_boundary = false;

  bool is_boundary() const { return m_is_boundary; }
};
