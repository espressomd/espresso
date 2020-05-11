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

#include <vector>

struct Observable_stat {
  /** Status flag for observable calculation. For 'analyze energy': 0
   *  re-initialize observable struct, else everything is fine,
   *  calculation can start. For 'analyze pressure' and 'analyze p_inst':
   *  0 or !(1+v_comp) re-initialize, else all OK.
   */
  int init_status;

  /** Array for observables on each node. */
  std::vector<double> data;

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

  /** number of doubles per data item */
  size_t chunk_size;

  void realloc_and_clear(size_t n_coulomb, size_t n_dipolar, size_t n_vs,
                         size_t c_size);

  void rescale(double volume, double time_step) {
    for (size_t i = 0; i < chunk_size; ++i)
      data[i] /= (volume * time_step * time_step);
    for (size_t i = chunk_size; i < data.size(); ++i)
      data[i] /= volume;
  }
};

/** Structure used only in the pressure and stress tensor calculation to
   distinguish
    non-bonded intra- and inter- molecular contributions. */
struct Observable_stat_non_bonded {
  /** Array for observables on each node. */
  std::vector<double> data_nb;

  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

  /** number of doubles per data item */
  size_t chunk_size_nb;

  void realloc_and_clear_non_bonded(size_t c_size);

  void rescale(double volume) {
    for (auto &value : data_nb)
      value /= volume;
  }
};

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
