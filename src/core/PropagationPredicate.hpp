/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "CellParticleIterator.hpp"
#include "ParticleIterator.hpp"

#include <boost/iterator/filter_iterator.hpp>

#include <iterator>

template <typename Predicate> struct PropagationPredicate {
  Predicate predicate;

  PropagationPredicate(Predicate pred) : predicate(pred) {}

  bool operator()(Particle const &p) { return predicate(p.propagation()); };
};

template <typename Predicate>
class ParticleRangeFiltered
    : public boost::iterator_range<boost::iterators::filter_iterator<
          PropagationPredicate<Predicate>, CellParticleIterator>> {
  using base_type = boost::iterator_range<boost::iterators::filter_iterator<
      PropagationPredicate<Predicate>, CellParticleIterator>>;

public:
  using base_type::base_type;
  auto size() const { return std::distance(this->begin(), this->end()); };
};
