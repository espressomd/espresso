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

#include <boost/iterator/iterator_facade.hpp>

#include <cassert>
#include <iterator>
#include <optional>
#include <utility>

namespace detail {
/* Detect the particle iterator type for a given cell iterator type. */
template <class CellIterator>
using particle_iterator_t =
    decltype((*std::declval<CellIterator>())->particles().begin());
/* Detect the particle type for a given cell iterator type. */
template <class CellIterator>
using particle_t =
    std::iterator_traits<particle_iterator_t<CellIterator>>::value_type;
} // namespace detail

template <typename BidirectionalIterator>
struct ParticleIterator
    : public boost::iterator_facade<ParticleIterator<BidirectionalIterator>,
                                    detail::particle_t<BidirectionalIterator>,
                                    boost::forward_traversal_tag> {
private:
  using base_type =
      boost::iterator_facade<ParticleIterator<BidirectionalIterator>,
                             detail::particle_t<BidirectionalIterator>,
                             boost::forward_traversal_tag>;
  using particle_iterator = detail::particle_iterator_t<BidirectionalIterator>;

  BidirectionalIterator m_cell, m_end;
  particle_iterator m_part;

public:
  ParticleIterator(BidirectionalIterator cell, BidirectionalIterator end)
      : m_cell(cell), m_end(end) {
    m_part = (m_cell != m_end) ? (*m_cell)->particles().begin()
                               : particle_iterator();
  }

  ParticleIterator(BidirectionalIterator end)
      : m_cell(end), m_end(end), m_part() {}

private:
  friend base_type::difference_type distance(ParticleIterator const &begin,
                                             ParticleIterator const &end) {
    if (begin == end)
      return 0;

    /* Remaining parts in this cell */
    auto dist = std::distance(begin.m_part, (*begin.m_cell)->particles().end());
    /* Now add the size of all cells between the next
       one and the last one */
    auto it = std::next(begin.m_cell);

    while (it != end.m_cell) {
      dist += static_cast<long>((*it)->particles().size());
      ++it;
    }

    return dist;
  }

  friend class boost::iterator_core_access;

  void increment() {
    assert(m_cell != m_end);

    ++m_part;
    /* If we are at the end of the particle range of the current cell,
     * we have to go to the next cell with particles. */
    if (m_part == (*m_cell)->particles().end()) {
      /* Find next cell with particles, without running over the end. */
      do {
        ++m_cell;
      } while ((m_cell != m_end) && ((*m_cell)->particles().empty()));

      /* If there is a cell, start go to its beginning. */
      m_part = (m_cell != m_end) ? (*m_cell)->particles().begin()
                                 : particle_iterator();
    }
  }

  bool equal(ParticleIterator const &rhs) const {
    return (m_cell == (rhs.m_cell)) && (m_part == rhs.m_part);
  }

  // True iff the inner particle iterator exposes the store-handle accessors of
  // @ref RowParticleIterator -- i.e. the production Cell::particles() iterator.
  // Generic/mock cell iterators (unit tests) do not, and take the plain
  // fallback below.
  static constexpr bool rebindable_part =
      requires(particle_iterator const &it) {
        it.current_store();
        it.current_row();
      };

  // Rebind this iterator's OWN reused view to the current (store, row) and
  // hand back a reference to it. Going through m_part's own dereference would
  // materialise/rebind the inner RowParticleIterator's cache instead, and that
  // inner iterator is REPLACED wholesale at every cell boundary (increment()
  // reassigns m_part), so its cache would be reconstructed once per cell --
  // costly when cells hold few particles. Keeping the view here means it
  // survives cell boundaries and is only ever rebound (two handle-field
  // writes). For non-rebindable (mock) cell iterators this falls back to the
  // plain by-reference dereference.
  auto &dereference() const {
    if constexpr (rebindable_part) {
      auto &self = const_cast<ParticleIterator &>(*this);
      if (not self.m_view) {
        self.m_view.emplace();
      }
      self.m_view->attach_to_store(*m_part.current_store(),
                                   m_part.current_row());
      return *self.m_view;
    } else {
      return *m_part;
    }
  }

  /** This iterator's reused view; rebound per dereference, never copied per
   *  element (see @c dereference), for the production rebindable path. Unused
   *  (and default-empty) on the generic fallback path. Mutable so the const
   *  dereference (boost iterator_facade contract) can rebind the
   *  logically-transient cache. */
  mutable std::optional<detail::particle_t<BidirectionalIterator>> m_view;
};
