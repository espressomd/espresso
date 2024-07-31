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

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <numeric>
#include <span>
#include <vector>

/** Observable for the pressure and energy. */
class Observable_stat {
  /** Array for observables on each node. */
  std::vector<double> m_data;
  /** Number of doubles per data item */
  std::size_t m_chunk_size;

  std::size_t get_non_bonded_offset(int type1, int type2) const;

  /** Get contribution from a non-bonded interaction */
  auto get_non_bonded_contribution(std::span<double> view, int type1,
                                   int type2) const {
    auto const offset = get_non_bonded_offset(type1, type2);
    return view.subspan(offset * m_chunk_size, m_chunk_size);
  }

public:
  Observable_stat(std::size_t chunk_size, std::size_t n_bonded, int max_type);

  auto get_chunk_size() const { return m_chunk_size; }

  /** Accumulate values.
   *  @param acc    Initial value for the accumulator.
   *  @param column Which column to sum up (only relevant for multi-dimensional
   *                observables).
   */
  double accumulate(double acc = 0.0, std::size_t column = 0ul) const {
    assert(column < m_chunk_size);
    if (m_chunk_size == 1ul)
      return std::accumulate(m_data.begin(), m_data.end(), acc);

    for (auto it = m_data.begin() + static_cast<std::ptrdiff_t>(column);
         it < m_data.end(); it += static_cast<std::ptrdiff_t>(m_chunk_size))
      acc += *it;
    return acc;
  }

  /** Rescale values */
  void rescale(double volume) {
    std::ranges::transform(m_data, m_data.begin(),
                           std::bind_front(std::multiplies{}, 1. / volume));
  }

  /** Contribution from linear and angular kinetic energy (accumulated). */
  std::span<double> kinetic;
  /** Contribution(s) from bonded interactions. */
  std::span<double> bonded;
  /** Contribution(s) from Coulomb interactions. */
  std::span<double> coulomb;
  /** Contribution(s) from dipolar interactions. */
  std::span<double> dipolar;
  /** Contribution from virtual sites (accumulated). */
  std::span<double> virtual_sites;
  /** Contribution from external fields (accumulated). */
  std::span<double> external_fields;
  /** Contribution(s) from non-bonded intramolecular interactions. */
  std::span<double> non_bonded_intra;
  /** Contribution(s) from non-bonded intermolecular interactions. */
  std::span<double> non_bonded_inter;

  /** Get contribution from a bonded interaction */
  std::span<double> bonded_contribution(int bond_id) const {
    auto const offset = m_chunk_size * static_cast<std::size_t>(bond_id);
    return {bonded.data() + offset, m_chunk_size};
  }

  void add_non_bonded_contribution(int type1, int type2, int molid1, int molid2,
                                   std::span<const double> data) {
    assert(data.size() == m_chunk_size);
    auto const view = (molid1 == molid2) ? non_bonded_intra : non_bonded_inter;
    auto const dest = get_non_bonded_contribution(view, type1, type2);

    std::ranges::transform(dest, data, dest.begin(), std::plus{});
  }

  void add_non_bonded_contribution(int type1, int type2, int molid1, int molid2,
                                   double data) {
    add_non_bonded_contribution(type1, type2, molid1, molid2, {&data, 1});
  }

  /** Get contribution from a non-bonded intramolecular interaction */
  auto non_bonded_intra_contribution(int type1, int type2) const {
    return get_non_bonded_contribution(non_bonded_intra, type1, type2);
  }

  /** Get contribution from a non-bonded intermolecular interaction */
  auto non_bonded_inter_contribution(int type1, int type2) const {
    return get_non_bonded_contribution(non_bonded_inter, type1, type2);
  }

  /** MPI reduction. */
  void mpi_reduce();
};
