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

#include <boost/range/numeric.hpp>

#include <utils/Span.hpp>
#include <utils/index.hpp>

#include <algorithm>
#include <utility>
#include <vector>

/** Observable for system statistics.
 *  Store unidimensional (energy, scalar pressure) and multi-dimensional
 *  (pressure tensor) properties of the system and provide accumulation and
 *  reduction functionality.
 */
class Observable_stat_base {
protected:
  /** Array for observables on each node. */
  std::vector<double> data;
  /** Number of doubles per data item */
  size_t m_chunk_size;

public:
  /** @param chunk_size Dimensionality of the data
   */
  explicit Observable_stat_base(size_t chunk_size) : m_chunk_size(chunk_size) {}

  /** Resize the observable */
  virtual void resize() = 0;

  /** Resize and zero out the observable */
  void resize_and_clear() {
    resize();
    std::fill(data.begin(), data.end(), 0);
  }

  /** Gather the contributions from the current MPI rank.
   *  @param[out] out Destination of the reduction.
   */
  void reduce(double *out) const;

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

protected:
  /** Get contribution from a non-bonded interaction */
  Utils::Span<double> non_bonded_contribution(double *base_pointer, int type1,
                                              int type2) const {
    extern int max_seen_particle_type;
    if (type1 > type2) {
      using std::swap;
      swap(type1, type2);
    }
    int offset = Utils::upper_triangular(type1, type2, max_seen_particle_type);
    return Utils::Span<double>(base_pointer + offset * m_chunk_size,
                               m_chunk_size);
  }

  /** Calculate the maximal number of non-bonded interaction pairs in the
   *  system.
   */
  size_t max_non_bonded_pairs() {
    extern int max_seen_particle_type;
    return static_cast<size_t>(
        (max_seen_particle_type * (max_seen_particle_type + 1)) / 2);
  }
};

/** Observable for the scalar pressure, pressure tensor and energy. */
class Observable_stat : public Observable_stat_base {
private:
  /** Whether this observable is a pressure or energy observable */
  bool m_pressure_obs;

public:
  explicit Observable_stat(size_t chunk_size, bool pressure_obs = true)
      : Observable_stat_base(chunk_size), m_pressure_obs(pressure_obs) {
    resize_and_clear();
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

  /** Resize the observable */
  void resize() final;

  /** Get contribution from a bonded interaction */
  Utils::Span<double> bonded_contribution(int bond_id) {
    return Utils::Span<double>(bonded.data() + m_chunk_size * bond_id,
                               m_chunk_size);
  }

  /** Get contribution from a non-bonded interaction */
  Utils::Span<double> non_bonded_contribution(int type1, int type2) const {
    return Observable_stat_base::non_bonded_contribution(non_bonded.data(),
                                                         type1, type2);
  }
};

/** Structure used only in the pressure calculation to distinguish
 *  non-bonded intra- and inter- molecular contributions.
 */
class Observable_stat_non_bonded : public Observable_stat_base {
public:
  explicit Observable_stat_non_bonded(size_t chunk_size)
      : Observable_stat_base(chunk_size) {
    resize_and_clear();
  }

  /** Contribution(s) from non-bonded intramolecular interactions. */
  Utils::Span<double> non_bonded_intra;
  /** Contribution(s) from non-bonded intermolecular interactions. */
  Utils::Span<double> non_bonded_inter;

  /** Resize the observable */
  void resize() final;

  /** Get contribution from a non-bonded intramolecular interaction */
  Utils::Span<double> non_bonded_intra_contribution(int type1,
                                                    int type2) const {
    return non_bonded_contribution(non_bonded_intra.data(), type1, type2);
  }

  /** Get contribution from a non-bonded intermolecular interaction */
  Utils::Span<double> non_bonded_inter_contribution(int type1,
                                                    int type2) const {
    return non_bonded_contribution(non_bonded_inter.data(), type1, type2);
  }
};

class Observable_stat_wrapper : public Observable_stat {
public:
  /** Observed statistic for the current MPI rank. */
  Observable_stat local;
  /** Flag to signal if the observable measures instantaneous pressure, i.e.
   *  the pressure with velocity compensation (half a time step), instead of
   *  the conventional pressure. Only relevant for NpT simulations.
   */
  bool v_comp;

  explicit Observable_stat_wrapper(size_t chunk_size, bool pressure_obs = true)
      : Observable_stat{chunk_size, pressure_obs}, local{chunk_size,
                                                         pressure_obs},
        v_comp(false) {}

  /** Gather the contributions from all MPI ranks. */
  void reduce() { local.reduce(data.data()); }
};

class Observable_stat_non_bonded_wrapper : public Observable_stat_non_bonded {
public:
  /** Observed statistic for the current MPI rank. */
  Observable_stat_non_bonded local;

  explicit Observable_stat_non_bonded_wrapper(size_t chunk_size)
      : Observable_stat_non_bonded{chunk_size}, local{chunk_size} {}

  /** Gather the contributions from all MPI ranks. */
  void reduce() { local.reduce(data.data()); }
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
