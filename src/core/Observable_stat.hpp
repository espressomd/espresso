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

struct Observable_stat_base {
  /** Get the first field in the observable */
  Utils::Span<double> first_field() {
    return Utils::Span<double>(data.data(), m_chunk_size);
  }

  /** Gather the contributions from all nodes */
  void reduce(Observable_stat_base *output) const;

  /** Accumulate values */
  double accumulate(double acc = 0.0, size_t from = 0) {
    if (m_chunk_size == 1 and from == 0)
      return boost::accumulate(data, acc);

    for (auto it = data.begin() + from; it != data.end(); ++it)
      acc += *it;
    return acc;
  }

  /** Accumulate values along one axis */
  double accumulate_along_dim(double acc = 0.0, size_t from = 0) {
    if (m_chunk_size == 1 and from == 0 and m_chunk_size == 1)
      return accumulate(acc);

    for (auto it = data.begin() + from; it < data.end(); it += m_chunk_size)
      acc += *it;
    return acc;
  }

protected:
  /** Array for observables on each node. */
  std::vector<double> data;
  /** number of doubles per data item */
  size_t m_chunk_size;
  /** Get contribution from a non-bonded interaction */
  double *nonbonded_ia(double *base_pointer, int type1, int type2) const {
    extern int max_seen_particle_type;
    if (type1 > type2) {
      using std::swap;
      swap(type1, type2);
    }
    return base_pointer +
           m_chunk_size *
               (((2 * max_seen_particle_type - 1 - type1) * type1) / 2 + type2);
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

/** Structure used for pressure, stress tensor and energy calculations. */
struct Observable_stat : Observable_stat_base {
  /** Flag to signal if the observable is initialized */
  bool is_initialized;
  /** Flag to signal if the observable measures instantaneous pressure, i.e.
   *  the pressure with velocity compensation (half a time step), instead of
   *  the conventional pressure. Only relevant for NpT simulations.
   */
  bool v_comp;
  /** number of Coulomb interactions */
  size_t n_coulomb;
  /** number of dipolar interactions */
  size_t n_dipolar;
  /** Number of virtual sites relative (rigid body) contributions */
  size_t n_virtual_sites;
  /** Number of external field contributions */
  const static size_t n_external_field = 1;

  /** start of bonded interactions. Right after the special ones */
  double *bonded;
  /** start of observables for non-bonded interactions. */
  double *non_bonded;
  /** start of observables for Coulomb interaction. */
  double *coulomb;
  /** start of observables for Coulomb interaction. */
  double *dipolar;
  /** Start of observables for virtual sites relative (rigid bodies) */
  double *virtual_sites;
  /** Start of observables for external fields */
  double *external_fields;

  /** Reinitialize the observable */
  void realloc_and_clear(size_t n_coulomb, size_t n_dipolar, size_t n_vs,
                         size_t c_size);

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
  double *bonded_ia(int bond_id) { return bonded + (m_chunk_size * bond_id); }

  /** Get contribution from a non-bonded interaction */
  double *nonbonded_ia(int type1, int type2) const {
    return Observable_stat_base::nonbonded_ia(non_bonded, type1, type2);
  }
};

/** Structure used only in the pressure and stress tensor calculation to
 *  distinguish non-bonded intra- and inter- molecular contributions.
 */
struct Observable_stat_non_bonded : Observable_stat_base {
private:
  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

public:
  /** Reinitialize the observable */
  void realloc_and_clear(size_t c_size);

  /** Rescale values */
  void rescale(double volume) {
    auto const factor = 1. / volume;
    for (auto &value : data)
      value *= factor;
  }

  /** Get contribution from a non-bonded intramolecular interaction */
  double *nonbonded_intra_ia(int type1, int type2) const {
    return nonbonded_ia(non_bonded_intra, type1, type2);
  }

  /** Get contribution from a non-bonded intermolecular interaction */
  double *nonbonded_inter_ia(int type1, int type2) const {
    return nonbonded_ia(non_bonded_inter, type1, type2);
  }
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
