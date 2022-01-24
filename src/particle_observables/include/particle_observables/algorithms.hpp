/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef SRC_PARTICLE_OBSERVABLES_ALGORITHMS_HPP
#define SRC_PARTICLE_OBSERVABLES_ALGORITHMS_HPP

/**
 * @file
 *
 * Generic algorithms for the calculation of particle
 * property derived observables.
 */

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

namespace ParticleObservables {
namespace detail {
struct One {
  template <class Particle> auto operator()(Particle const &p) const {
    return 1;
  }
};

template <class ValueOp, class WeightOp> struct WeightedSum {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    using particle_type = typename ParticleRange::value_type;
    using value_op_type = decltype(ValueOp{}(std::declval<particle_type>()));
    using weight_op_type = decltype(WeightOp{}(std::declval<particle_type>()));
    auto func = [](auto sum, auto const &p) {
      auto const w = WeightOp{}(p);
      return std::make_pair(sum.first + ValueOp{}(p)*w, sum.second + w);
    };

    return std::accumulate(std::begin(particles), std::end(particles),
                           std::pair<value_op_type, weight_op_type>(), func);
  }
};
} // namespace detail

template <class ValueOp, class WeightOp> struct WeightedSum {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    return detail::WeightedSum<ValueOp, WeightOp>()(particles).first;
  }
};

template <class ValueOp> struct Sum {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    return detail::WeightedSum<ValueOp, detail::One>()(particles).first;
  }
};

template <class ValueOp, class WeightOp> struct WeightedAverage {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    auto const ws = detail::WeightedSum<ValueOp, WeightOp>()(particles);
    return (ws.second) ? ws.first / ws.second : ws.first;
  }
};

template <class ValueOp> struct Average {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    return WeightedAverage<ValueOp, detail::One>()(particles);
  }
};

template <class ValueOp> struct Map {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) const {
    using particle_type = typename ParticleRange::value_type;
    using value_op_type = decltype(ValueOp{}(std::declval<particle_type>()));
    std::vector<value_op_type> res;
    std::transform(std::begin(particles), std::end(particles),
                   std::back_inserter(res),
                   [](auto const &p) { return ValueOp{}(p); });
    return res;
  }
};
} // namespace ParticleObservables
#endif // SRC_PARTICLE_OBSERVABLES_ALGORITHMS_HPP
