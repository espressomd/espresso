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

#include "particle_store/ParticleStore.hpp"

#include <cstddef>
#include <iterator>

/**
 * @brief Lightweight view over a @ref Cell's contiguous store-row range.
 *
 * A @ref Cell's committed particles occupy a CONTIGUOUS range
 * @c [offset, offset+count) of the @ref ParticleStore, established by the
 * permutation rebuild (counting sort by cell). @ref Cell::rows() hands out a
 * @c CellRowSpan over that range, backed by pure arithmetic.
 *
 * MID-EPOCH REMOVAL: a range cannot swap-remove, so a row dropped between
 * rebuilds
 * (@ref CellParticleStorage::drop_row) is marked PENDING-REMOVED on the store
 * (@ref ParticleStore::mark_pending_removal). The LIVE surface of this span --
 * @ref begin / @ref end / @ref size / @ref empty -- skips such rows, so a
 * removed particle is invisible to iteration and id-resolution immediately, and
 * the next rebuild resolves it by not carrying the row into the new generation.
 * @ref offset / @ref count expose the RAW (removal-agnostic) range for the
 * removal-free contexts (the resort mutation loops, which iterate raw positions
 * and mark rows pending, and the hot physics loops, which never run mid-removal
 * so the raw range equals the live range there).
 *
 * Identity histories are removal-free, so the live and raw views coincide on
 * the identity path.
 */
class CellRowSpan {
  std::size_t m_offset;
  std::size_t m_count;
  ParticleStore const *m_store;

  // The next raw row at or after `raw` that is not pending-removed, clamped to
  // the range end (offset + count) when none remains.
  std::size_t skip_removed(std::size_t raw) const {
    auto const end = m_offset + m_count;
    while (raw < end and m_store != nullptr and
           m_store->is_pending_removal(static_cast<int>(raw))) {
      ++raw;
    }
    return raw;
  }

public:
  /** @brief Forward iterator over the LIVE rows (skips pending-removed). Yields
   *  the store row index (@c int) for each surviving row, in ascending order.
   */
  class const_iterator {
    std::size_t m_raw;
    CellRowSpan const *m_span;

  public:
    using value_type = int;
    using reference = int;
    using pointer = void;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    const_iterator() : m_raw(0u), m_span(nullptr) {}
    const_iterator(std::size_t raw, CellRowSpan const *span)
        : m_raw(raw), m_span(span) {}

    int operator*() const { return static_cast<int>(m_raw); }
    const_iterator &operator++() {
      m_raw = m_span->skip_removed(m_raw + 1u);
      return *this;
    }
    const_iterator operator++(int) {
      auto const copy = *this;
      ++(*this);
      return copy;
    }
    bool operator==(const_iterator const &rhs) const {
      return m_raw == rhs.m_raw;
    }
    bool operator!=(const_iterator const &rhs) const {
      return not(*this == rhs);
    }
  };

  CellRowSpan(std::size_t offset, std::size_t count, ParticleStore const *store)
      : m_offset(offset), m_count(count), m_store(store) {}

  /** @brief Raw range start (removal-agnostic). */
  std::size_t offset() const { return m_offset; }
  /** @brief Raw committed row count (removal-agnostic). */
  std::size_t count() const { return m_count; }

  const_iterator begin() const { return {skip_removed(m_offset), this}; }
  const_iterator end() const { return {m_offset + m_count, this}; }

  /** @brief Number of LIVE rows (raw count minus pending-removed). */
  std::size_t size() const {
    if (m_store == nullptr) {
      return m_count;
    }
    std::size_t live = 0u;
    for (std::size_t raw = m_offset; raw < m_offset + m_count; ++raw) {
      if (not m_store->is_pending_removal(static_cast<int>(raw))) {
        ++live;
      }
    }
    return live;
  }

  bool empty() const { return begin() == end(); }
};
