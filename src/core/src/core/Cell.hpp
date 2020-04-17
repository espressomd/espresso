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
#ifndef CORE_CELL_HPP
#define CORE_CELL_HPP

#include "Particle.hpp"
#include "ParticleList.hpp"

#include <utils/Span.hpp>

#include <boost/range/iterator_range.hpp>

#include <functional>
#include <vector>

template <class CellRef> class Neighbors {
  using storage_type = std::vector<CellRef>;

public:
  using value_type = typename storage_type::value_type;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;
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

  Neighbors(Utils::Span<const CellRef> red_neighbors,
            Utils::Span<const CellRef> black_neighbors) {
    m_neighbors.resize(red_neighbors.size() + black_neighbors.size());
    m_red_black_divider = std::copy(red_neighbors.begin(), red_neighbors.end(),
                                    m_neighbors.begin());
    std::copy(black_neighbors.begin(), black_neighbors.end(),
              m_red_black_divider);
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

  ParticleList m_particles;

public:
  /** Particles */
  auto &particles() { return m_particles; }
  auto const &particles() const { return m_particles; }

  neighbors_type m_neighbors;

  /** Interaction pairs */
  std::vector<std::pair<Particle *, Particle *>> m_verlet_list;

  /**
   * @brief All neighbors of the cell.
   */
  neighbors_type &neighbors() { return m_neighbors; }
};

#endif
