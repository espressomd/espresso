/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef ESPRESSO_SRC_CORE_PARTICLE_RANGE_HPP
#define ESPRESSO_SRC_CORE_PARTICLE_RANGE_HPP

#include "Particle.hpp"
#include "ParticleIterator.hpp"
#include "PropagationModes.hpp"
#include "cell_system/Cell.hpp"

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <cassert>

template <PropagationMode criterion> struct PropagationPredicate {
  bool operator()(Particle const &p) {
    // return p.p.propagation == criterion;
    return true;
  };
};

template <PropagationMode criterion>
class ParticleRangeFiltered
    : public boost::iterator_range<boost::iterators::filter_iterator<
          PropagationPredicate<criterion>, ParticleIterator<Cell **>>> {
  using base_type = boost::iterator_range<boost::iterators::filter_iterator<
      PropagationPredicate<criterion>, ParticleIterator<Cell **>>>;

public:
  using base_type::base_type;
  auto  size() const{ return std::distance(this->begin(), this->end()); } ;
};

using CellParticleIterator = ParticleIterator<Cell **>;

/**
 * @brief A range of particles.
 *
 * This is a boost::iterator_range with the addition that the
 * size of the range is cached.
 */
class ParticleRange : public boost::iterator_range<CellParticleIterator> {
  using base_type = boost::iterator_range<CellParticleIterator>;

public:
  using base_type::base_type;

  base_type::size_type size() const {
    if (m_size < 0) {
      m_size = distance(begin(), end());
    }

    return assert(m_size >= 0), static_cast<base_type::size_type>(m_size);
  }

  template <PropagationMode criterion>
  ParticleRangeFiltered<criterion> filter() const {
    return {boost::make_filter_iterator<PropagationPredicate<criterion>>(
                begin(), end()),
            boost::make_filter_iterator<PropagationPredicate<criterion>>(
                end(), end())};
  }

private:
  base_type::difference_type mutable m_size = -1;
};

using ParticleRangeDefault =
    ParticleRangeFiltered<PropagationMode::TRANS_SYSTEM_DEFAULT>;
using ParticleRangeLangevin =
    ParticleRangeFiltered<PropagationMode::TRANS_LANGEVIN>;
using ParticleRangeStokesian =
    ParticleRangeFiltered<PropagationMode::TRANS_STOKESIAN>;

#endif
