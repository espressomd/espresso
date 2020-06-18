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
#ifndef ESPRESSO_OBSERVABLE_STAT_HPP
#define ESPRESSO_OBSERVABLE_STAT_HPP

#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/numeric.hpp>

#include <utils/Span.hpp>
#include <utils/index.hpp>

#include <algorithm>
#include <utility>
#include <vector>

/** Observable for the scalar pressure, pressure tensor and energy. */
class Observable_stat {
  /** Array for observables on each node. */
  std::vector<double> data;
  /** Number of doubles per data item */
  size_t m_chunk_size;

  /** Whether this observable is a pressure or energy observable */
  bool m_pressure_obs;

  /** Calculate the maximal number of non-bonded interaction pairs in the
   *  system.
   */
  static size_t max_non_bonded_pairs() {
    extern int max_seen_particle_type;
    return static_cast<size_t>(
        (max_seen_particle_type * (max_seen_particle_type + 1)) / 2);
  }

  /** Get contribution from a non-bonded interaction */
  Utils::Span<double> non_bonded_contribution(Utils::Span<double> base_pointer,
                                              int type1, int type2) const {
    extern int max_seen_particle_type;
    int offset = Utils::upper_triangular(
        std::min(type1, type2), std::max(type1, type2), max_seen_particle_type);
    return {base_pointer.begin() + offset * m_chunk_size, m_chunk_size};
  }

public:
  explicit Observable_stat(size_t chunk_size, bool pressure_obs = true)
      : m_chunk_size(chunk_size), m_pressure_obs(pressure_obs) {
    resize();
  }

  /** Accumulate values.
   *  @param acc    Initial value for the accumulator.
   *  @param column Which column to sum up (only relevant for multi-dimensional
   *                observables).
   */
  double accumulate(double acc = 0.0, size_t column = 0) {
    assert(column < m_chunk_size);
    if (m_chunk_size == 1)
      return boost::accumulate(data, acc);

    for (auto it = data.begin() + column; it < data.end(); it += m_chunk_size)
      acc += *it;
    return acc;
  }

  /** Rescale values */
  void rescale(double volume) {
    auto const factor = 1. / volume;
    for (auto &e : data) {
      e *= factor;
    }
  }

  /** Reduce contributions from all MPI ranks. */
  void reduce(boost::mpi::communicator const &comm) {
    BOOST_MPI_CHECK_RESULT(
        MPI_Reduce, ((comm.rank() == 0) ? MPI_IN_PLACE : data.data(),
                     data.data(), data.size(), MPI_DOUBLE, MPI_SUM, 0, comm));
  }

  /** Contribution from linear and angular kinetic energy (accumulated). */
  Utils::Span<double> kinetic;
  /** Contribution(s) from bonded interactions. */
  Utils::Span<double> bonded;
  /** Contribution(s) from non-bonded interactions. */
  Utils::Span<double> non_bonded;
  /** Contribution(s) from Coulomb interactions. */
  Utils::Span<double> coulomb;
  /** Contribution(s) from dipolar interactions. */
  Utils::Span<double> dipolar;
  /** Contribution(s) from virtual sites (accumulated). */
  Utils::Span<double> virtual_sites;
  /** Contribution from external fields (accumulated). */
  Utils::Span<double> external_fields;
  /** Contribution(s) from non-bonded intramolecular interactions. */
  Utils::Span<double> non_bonded_intra;
  /** Contribution(s) from non-bonded intermolecular interactions. */
  Utils::Span<double> non_bonded_inter;

  /** Resize the observable */
  void resize();

  /** Get contribution from a bonded interaction */
  Utils::Span<double> bonded_contribution(int bond_id) {
    return Utils::Span<double>(bonded.data() + m_chunk_size * bond_id,
                               m_chunk_size);
  }

  /** Get contribution from a non-bonded interaction */
  Utils::Span<double> non_bonded_contribution(int type1, int type2) const {
    return non_bonded_contribution(non_bonded, type1, type2);
  }

  /** Get contribution from a non-bonded intramolecular interaction */
  Utils::Span<double> non_bonded_intra_contribution(int type1,
                                                    int type2) const {
    return non_bonded_contribution(non_bonded_intra, type1, type2);
  }

  /** Get contribution from a non-bonded intermolecular interaction */
  Utils::Span<double> non_bonded_inter_contribution(int type1,
                                                    int type2) const {
    return non_bonded_contribution(non_bonded_inter, type1, type2);
  }
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
