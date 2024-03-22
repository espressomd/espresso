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
#ifndef ESPRESSO_OBSERVABLE_STAT_HPP
#define ESPRESSO_OBSERVABLE_STAT_HPP

#include <boost/range/algorithm/transform.hpp>
#include <boost/range/numeric.hpp>

#include <utils/Span.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

/** Observable for the pressure and energy. */
class Observable_stat {
  /** Array for observables on each node. */
  std::vector<double> m_data;
  /** Number of doubles per data item */
  std::size_t m_chunk_size;

  /** Get contribution from a non-bonded interaction */
  Utils::Span<double>
  get_non_bonded_contribution(Utils::Span<double> base_pointer, int type1,
                              int type2) const;

public:
  explicit Observable_stat(std::size_t chunk_size);

  auto get_chunk_size() const { return m_chunk_size; }

  /** Accumulate values.
   *  @param acc    Initial value for the accumulator.
   *  @param column Which column to sum up (only relevant for multi-dimensional
   *                observables).
   */
  double accumulate(double acc = 0.0, std::size_t column = 0) const {
    assert(column < m_chunk_size);
    if (m_chunk_size == 1)
      return boost::accumulate(m_data, acc);

    for (auto it = m_data.begin() + static_cast<std::ptrdiff_t>(column);
         it < m_data.end(); it += static_cast<std::ptrdiff_t>(m_chunk_size))
      acc += *it;
    return acc;
  }

  /** Rescale values */
  void rescale(double volume) {
    auto const fac = 1. / volume;
    boost::transform(m_data, m_data.begin(), [fac](auto e) { return e * fac; });
  }

  /** Contribution from linear and angular kinetic energy (accumulated). */
  Utils::Span<double> kinetic;
  /** Contribution(s) from bonded interactions. */
  Utils::Span<double> bonded;
  /** Contribution(s) from Coulomb interactions. */
  Utils::Span<double> coulomb;
  /** Contribution(s) from dipolar interactions. */
  Utils::Span<double> dipolar;
  /** Contribution from virtual sites (accumulated). */
  Utils::Span<double> virtual_sites;
  /** Contribution from external fields (accumulated). */
  Utils::Span<double> external_fields;
  /** Contribution(s) from non-bonded intramolecular interactions. */
  Utils::Span<double> non_bonded_intra;
  /** Contribution(s) from non-bonded intermolecular interactions. */
  Utils::Span<double> non_bonded_inter;

  /** Get contribution from a bonded interaction */
  Utils::Span<double> bonded_contribution(int bond_id) const {
    auto const offset = m_chunk_size * static_cast<std::size_t>(bond_id);
    return {bonded.data() + offset, m_chunk_size};
  }

  void add_non_bonded_contribution(int type1, int type2, int molid1, int molid2,
                                   Utils::Span<const double> data) {
    assert(data.size() == m_chunk_size);
    auto const span = (molid1 == molid2) ? non_bonded_intra : non_bonded_inter;
    auto const dest = get_non_bonded_contribution(span, type1, type2);

    boost::transform(dest, data, dest.begin(), std::plus<>{});
  }

  void add_non_bonded_contribution(int type1, int type2, int molid1, int molid2,
                                   double data) {
    add_non_bonded_contribution(type1, type2, molid1, molid2, {&data, 1});
  }

  /** Get contribution from a non-bonded intramolecular interaction */
  Utils::Span<double> non_bonded_intra_contribution(int type1,
                                                    int type2) const {
    return get_non_bonded_contribution(non_bonded_intra, type1, type2);
  }

  /** Get contribution from a non-bonded intermolecular interaction */
  Utils::Span<double> non_bonded_inter_contribution(int type1,
                                                    int type2) const {
    return get_non_bonded_contribution(non_bonded_inter, type1, type2);
  }

  /** MPI reduction. */
  void mpi_reduce();
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
