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

#include <utility>
#include <vector>

/** Cache for system statistics.
 *  Store unidimensional (energy, scalar pressure) and multi-dimensional
 *  (stress tensors) properties of the system and provide accumulation and
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

  /** Reinitialize the observable */
  virtual void realloc_and_clear() = 0;

  /** Gather the contributions from all nodes */
  void reduce(Observable_stat_base *output) const;

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

protected:
  /** Reinitialize the observable */
  void realloc_and_clear(size_t n_items) {
    data.resize(m_chunk_size * n_items);
    for (auto &value : data)
      value = 0;
  }

  /** Get contribution from a non-bonded interaction */
  Utils::Span<double> non_bonded_contribution(double *base_pointer, int type1,
                                              int type2) const {
    extern int max_seen_particle_type;
    if (type1 > type2) {
      using std::swap;
      swap(type1, type2);
    }
    int offset = ((2 * max_seen_particle_type - 1 - type1) * type1) / 2 + type2;
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

private:
  /** Register this observable. */
  virtual void register_obs() = 0;
};

/** Structure used to cache the results of the scalar pressure, stress tensor
 *  and energy calculations.
 */
class Observable_stat : public Observable_stat_base {
private:
  /** Callback for the Coulomb interaction */
  size_t (*m_get_n_coulomb)();
  /** Callback for the dipolar interaction */
  size_t (*m_get_n_dipolar)();

public:
  explicit Observable_stat(size_t chunk_size, size_t (*get_n_coulomb)(),
                           size_t (*get_n_dipolar)())
      : Observable_stat_base(chunk_size), m_get_n_coulomb(get_n_coulomb),
        m_get_n_dipolar(get_n_dipolar), is_initialized(false), v_comp(false) {
    register_obs();
    realloc_and_clear();
  }
  /** Flag to signal if the observable is initialized */
  bool is_initialized;
  /** Flag to signal if the observable measures instantaneous pressure, i.e.
   *  the pressure with velocity compensation (half a time step), instead of
   *  the conventional pressure. Only relevant for NpT simulations.
   */
  bool v_comp;
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

  /** Reinitialize the observable */
  void realloc_and_clear() final;

  /** Rescale values */
  void rescale(double volume, double time_step) {
    auto const factor1 = 1. / (volume * time_step * time_step);
    auto const factor2 = 1. / volume;
    for (auto it = data.begin(); it != data.begin() + m_chunk_size; ++it)
      *it *= factor1;
    for (auto it = data.begin() + m_chunk_size; it != data.end(); ++it)
      *it *= factor2;
  }

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

private:
  /** Register this observable. */
  void register_obs() final;
};

/** Structure used only in the pressure and stress tensor calculation to
 *  distinguish non-bonded intra- and inter- molecular contributions.
 */
class Observable_stat_non_bonded : public Observable_stat_base {
public:
  explicit Observable_stat_non_bonded(size_t chunk_size)
      : Observable_stat_base(chunk_size) {
    register_obs();
    realloc_and_clear();
  }
  /** Contribution(s) from non-bonded intramolecular interactions. */
  Utils::Span<double> non_bonded_intra;
  /** Contribution(s) from non-bonded intermolecular interactions. */
  Utils::Span<double> non_bonded_inter;

  /** Reinitialize the observable */
  void realloc_and_clear() final;

  /** Rescale values */
  void rescale(double volume) {
    auto const factor = 1. / volume;
    for (auto &value : data)
      value *= factor;
  }

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

private:
  /** Register this observable. */
  void register_obs() final;
};

/** Invalidate observables.
 *  This function is called whenever the system has changed in such a way
 *  that the cached observables are no longer accurate or that the size of
 *  the cache is no longer suitable (e.g. after addition of a new actor or
 *  interaction).
 */
void invalidate_obs();

/** Resize all observables.
 *  This function is called whenever the number of Coulomb actors, dipolar
 *  actors, bonded interactions or non-nonded interactions changed. This is
 *  necessary because the data array must be resized accordingly.
 */
void realloc_and_clear_all_obs();

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
