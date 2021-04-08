/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include <particle_observables/observable.hpp>

#include "Observable.hpp"
#include "Particle.hpp"
#include "ParticleTraits.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/flatten.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <cstddef>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace Observables {

using ParticleReferenceRange =
    Utils::Span<std::reference_wrapper<const Particle>>;

/** %Particle-based observable.
 *
 *  Base class for observables extracting raw data from particle subsets and
 *  returning either the data or a statistic derived from it.
 */
class PidObservable : virtual public Observable {
  /** Identifiers of particles measured by this observable */
  std::vector<int> m_ids;

  virtual std::vector<double>
  evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &traits) const = 0;

public:
  explicit PidObservable(std::vector<int> ids) : m_ids(std::move(ids)) {}
  std::vector<double> operator()() const final;
  std::vector<int> const &ids() const { return m_ids; }
};

namespace detail {
/**
 * Recursive implementation for finding the shape of a given `std::vector` of
 * types. A vector of extents is constructed starting at
 * the template specialization for `std::vector<T>`.
 */
template <class T> struct shape_impl;

template <> struct shape_impl<double> {
  static std::vector<size_t> eval(size_t /* n_part */) { return {1}; }
};
template <class _, size_t N> struct shape_impl<Utils::Vector<_, N>> {
  static std::vector<size_t> eval(size_t /* n_part */) { return {N}; }
};
template <class T> struct shape_impl<std::vector<T>> {
  static std::vector<size_t> eval(size_t n_part) {
    std::vector<size_t> ret{n_part};
    boost::copy(shape_impl<T>::eval(n_part), std::back_inserter(ret));

    return ret;
  }
};
} // namespace detail

/**
 * This class implements an interface to the `particle_observables` library that
 * implements necessary algorithms needed for observables that are based on
 * single particle properties.
 * @tparam ObsType An observables composed of an algorithm from
 * src/particle_observables/include/particle_observables/algorithms.hpp and two
 * particle properties.
 *
 *  Example usage:
 *  @code{.cpp}
 *  using namespace ParticleObservables;
 *  using CenterOfMass = ParticleObservable<WeightedAverage<Position, Mass>>;
 *  @endcode
 */
template <class ObsType> class ParticleObservable : public PidObservable {
public:
  using PidObservable::PidObservable;
  std::vector<size_t> shape() const override {
    using std::declval;

    return detail::shape_impl<decltype(declval<ObsType>()(
        declval<ParticleReferenceRange>()))>::eval(ids().size());
  }

  std::vector<double>
  evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    std::vector<double> res;
    Utils::flatten(ObsType{}(particles), std::back_inserter(res));
    return res;
  }
};

} // namespace Observables
#endif
