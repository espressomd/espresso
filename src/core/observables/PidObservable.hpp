/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#include <particle_observables/observable.hpp>

#include "Observable.hpp"
#include "Particle.hpp"
#include "ParticleTraits.hpp"

#include <utils/Vector.hpp>
#include <utils/flatten.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace Observables {

using ParticleReferenceRange =
    std::vector<std::reference_wrapper<Particle const>>;

/** Particle-based observable.
 *
 *  Base class for observables extracting raw data from particle subsets and
 *  returning either the data or a statistic derived from it.
 */
class PidObservable : virtual public Observable {
  /** Identifiers of particles measured by this observable */
  std::vector<int> m_ids;

  virtual std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const = 0;

public:
  explicit PidObservable(std::vector<int> ids) : m_ids(std::move(ids)) {}
  std::vector<double>
  operator()(boost::mpi::communicator const &comm) const final;
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
  static std::vector<std::size_t> eval(std::size_t /* n_part */) { return {1}; }
};
template <class _, std::size_t N> struct shape_impl<Utils::Vector<_, N>> {
  static std::vector<std::size_t> eval(std::size_t /* n_part */) { return {N}; }
};
template <class T> struct shape_impl<std::vector<T>> {
  static std::vector<std::size_t> eval(std::size_t n_part) {
    std::vector<std::size_t> ret{n_part};
    boost::copy(shape_impl<T>::eval(n_part), std::back_inserter(ret));

    return ret;
  }
};
template <class T, class U> struct shape_impl<std::pair<T, U>> {
  static std::vector<std::size_t> eval(std::size_t n_part) {
    return shape_impl<T>::eval(n_part);
  }
};

inline auto get_argsort(boost::mpi::communicator const &comm,
                        std::vector<int> const &local_pids,
                        std::vector<int> const &sorted_pids) {
  std::vector<unsigned int> argsort{};

  std::vector<std::vector<int>> global_pids;
  boost::mpi::gather(comm, local_pids, global_pids, 0);
  if (comm.rank() == 0) {
    auto const n_part = sorted_pids.size();
    std::vector<int> unsorted_pids;
    unsorted_pids.reserve(n_part);
    for (auto const &vec : global_pids) {
      for (auto const pid : vec) {
        unsorted_pids.emplace_back(pid);
      }
    }
    // get vector of indices that sorts the data vectors
    std::vector<unsigned int> iota(n_part);
    std::iota(iota.begin(), iota.end(), 0u);
    argsort.reserve(n_part);
    auto const pid_begin = std::begin(unsorted_pids);
    auto const pid_end = std::end(unsorted_pids);
    for (auto const pid : sorted_pids) {
      auto const pid_pos = std::find(pid_begin, pid_end, pid);
      auto const i =
          static_cast<std::size_t>(std::distance(pid_begin, pid_pos));
      argsort.emplace_back(iota[i]);
    }
  }
  return argsort;
}

/** Get the positions of all the particles in the system in the order
 * specified by \a sorted_pids. Only the head node returns a vector
 * of positions, all other nodes return an empty vector.
 * This requires several MPI communications to construct.
 */
inline auto
get_all_particle_positions(boost::mpi::communicator const &comm,
                           ParticleReferenceRange const &local_particles,
                           std::vector<int> const &sorted_pids,
                           ParticleObservables::traits<Particle> const &traits,
                           bool use_folded_positions = false) {
  using pos_type = decltype(traits.position(std::declval<Particle>()));
  std::vector<pos_type> local_positions{};
  std::vector<int> local_pids{};
  local_positions.reserve(local_particles.size());
  local_pids.reserve(local_particles.size());

  for (auto const &particle : local_particles) {
    if (use_folded_positions) {
      local_positions.emplace_back(traits.position_folded(particle));
    } else {
      local_positions.emplace_back(traits.position(particle));
    }
    local_pids.emplace_back(traits.id(particle));
  }

  auto const argsort = detail::get_argsort(comm, local_pids, sorted_pids);

  std::vector<std::vector<pos_type>> global_positions{};
  global_positions.reserve(static_cast<std::size_t>(comm.size()));
  boost::mpi::gather(comm, local_positions, global_positions, 0);

  if (comm.rank() != 0) {
    return std::vector<pos_type>();
  }

  std::vector<pos_type> global_positions_flattened{};
  global_positions_flattened.reserve(sorted_pids.size());
  for (auto const &vec : global_positions) {
    for (auto const &pos : vec) {
      global_positions_flattened.emplace_back(pos);
    }
  }
  assert(global_positions_flattened.size() == sorted_pids.size());

  std::vector<pos_type> positions_sorted{};
  positions_sorted.reserve(sorted_pids.size());
  for (auto const i : argsort) {
    positions_sorted.emplace_back(global_positions_flattened[i]);
  }

  return positions_sorted;
}
} // namespace detail

/**
 * This class implements an interface to the `particle_observables` library that
 * implements necessary algorithms needed for observables that are based on
 * single particle properties.
 * @tparam ObsType An observables composed of an algorithm from
 * @ref src/particle_observables/include/particle_observables/algorithms.hpp
 * and two particle properties.
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
  std::vector<std::size_t> shape() const override {
    using std::declval;

    return detail::shape_impl<decltype(declval<ObsType>()(
        declval<ParticleReferenceRange const &>()))>::eval(ids().size());
  }

  template <typename T> struct is_map : std::false_type {};
  template <typename T>
  struct is_map<ParticleObservables::Map<T>> : std::true_type {};

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           ParticleObservables::traits<Particle> const &) const override {
    if constexpr (is_map<ObsType>::value) {
      std::vector<double> local_traits{};
      local_traits.reserve(local_particles.size());
      Utils::flatten(ObsType{}(local_particles),
                     std::back_inserter(local_traits));
      std::vector<int> local_pids{};
      local_pids.reserve(local_particles.size());
      Utils::flatten(ParticleObservables::Identities{}(local_particles),
                     std::back_inserter(local_pids));

      std::vector<std::vector<double>> global_traits{};
      boost::mpi::gather(comm, local_traits, global_traits, 0);

      auto const argsort = detail::get_argsort(comm, local_pids, ids());

      if (comm.rank() != 0) {
        return {};
      }

      // get total size of the global traits vector
      auto const size = std::accumulate(
          global_traits.begin(), global_traits.end(), 0u,
          [](auto const acc, auto const &vec) { return acc + vec.size(); });

      auto const n_dims = size / ids().size();

      std::vector<double> global_traits_flattened{};
      global_traits_flattened.reserve(size);

      for (auto const &vec : global_traits) {
        for (auto const val : vec) {
          global_traits_flattened.emplace_back(val);
        }
      }

      std::vector<double> output{};
      output.reserve(size);

      for (auto const i : argsort) {
        for (std::size_t j = 0; j < n_dims; ++j) {
          output.emplace_back(global_traits_flattened[i * n_dims + j]);
        }
      }
      return output;
    } else {
      auto observable = ObsType{};
      auto const local_result = observable(local_particles);
      std::remove_const_t<decltype(local_result)> result{};
      boost::mpi::reduce(comm, local_result, result, observable, 0);

      return result.first;
    }
  }
};

} // namespace Observables
