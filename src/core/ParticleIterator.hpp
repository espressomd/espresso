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
  int m_part_id;

public:
  ParticleIterator(BidirectionalIterator cell, BidirectionalIterator end,
                   int part_id)
      : m_cell(cell), m_end(end), m_part_id(part_id) {
    /* Jump to first actual particle */
    if ((m_cell != m_end) && (*m_cell)->particles().empty()) {
      increment();
    }
  }

private:
  friend typename base_type::difference_type
  distance(ParticleIterator const &begin, ParticleIterator const &end) {
    if (begin == end)
      return 0;

    /* Remaining parts in this cell */
    auto dist = ((*begin.m_cell)->particles().size() - begin.m_part_id);
    /* Now add the size of all cells between the next
       one and the last one */
    auto it = std::next(begin.m_cell);

    while (it != end.m_cell) {
      dist += (*it)->particles().size();
      ++it;
    }

    /* Remaining in the last cell */
    dist += end.m_part_id;

    return dist;
  }

  friend class boost::iterator_core_access;

  void increment() {
    /* If we are not at the end of the cells,
       there actually are particles in this cells
       and if we are not at the last particle in the
       cell we can just increment the particle id.
    */
    if ((m_cell != m_end) && (not(*m_cell)->particles().empty()) &&
        (m_part_id < ((*m_cell)->particles().size() - 1))) {
      /* Next part in same cell */
      ++m_part_id;
    } else {
      m_part_id = 0;

      if (m_cell != m_end)
        ++m_cell;

      /* Find next cell with particles */
      while ((m_cell != m_end) && ((*m_cell)->particles().empty())) {
        ++m_cell;
      }
    }
  }

  bool equal(ParticleIterator const &rhs) const {
    return (m_cell == (rhs.m_cell)) && (m_part_id == rhs.m_part_id);
  }

  Particle &dereference() const { return (*m_cell)->particles()[m_part_id]; }
};

#endif
