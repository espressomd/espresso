/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CORE_PARTICLE_ITERATOR_HPP
#define CORE_PARTICLE_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <iterator>

namespace detail {
/* Detect the particle iterator type for a given cell iterator type. */
template <class CellIterator>
using particle_iterator_t =
    decltype((*std::declval<CellIterator>())->particles().begin());
/* Detect the particle type for a given cell iterator type. */
template <class CellIterator>
using particle_t = typename std::iterator_traits<
    particle_iterator_t<CellIterator>>::value_type;
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
  friend typename base_type::difference_type
  distance(ParticleIterator const &begin, ParticleIterator const &end) {
    if (begin == end)
      return 0;

    /* Remaining parts in this cell */
    auto dist = std::distance(begin.m_part, (*begin.m_cell)->particles().end());
    /* Now add the size of all cells between the next
       one and the last one */
    auto it = std::next(begin.m_cell);

    while (it != end.m_cell) {
      dist += (*it)->particles().size();
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

  Particle &dereference() const { return *m_part; }
};

#endif
