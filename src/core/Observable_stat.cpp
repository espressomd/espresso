/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "Observable_stat.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/mpi/communicator.hpp>

extern boost::mpi::communicator comm_cart;

/** Calculate the maximal number of non-bonded interaction pairs in the system.
 */
size_t max_non_bonded_pairs() {
  return static_cast<size_t>(
      (max_seen_particle_type * (max_seen_particle_type + 1)) / 2);
}

void Observable_stat::realloc_and_clear(size_t n_coulomb, size_t n_dipolar,
                                        size_t n_vs, size_t c_size) {
  // Number of doubles per interaction (pressure=1, stress tensor=9,...)
  chunk_size = c_size;

  auto const n_bonded = bonded_ia_params.size();
  auto const n_non_bonded = max_non_bonded_pairs();

  // Number of doubles to store pressure in
  size_t const total = chunk_size * (1 + n_bonded + n_non_bonded + n_coulomb +
                                     n_dipolar + n_vs + n_external_field);

  // Allocate mem for the double list
  data.resize(total);

  // Number of chunks for different interaction types
  this->n_coulomb = n_coulomb;
  this->n_dipolar = n_dipolar;
  this->n_virtual_sites = n_vs;
  // Pointers to the start of different contributions
  bonded = data.data() + chunk_size;
  non_bonded = bonded + chunk_size * n_bonded;
  coulomb = non_bonded + chunk_size * n_non_bonded;
  dipolar = coulomb + chunk_size * n_coulomb;
  virtual_sites = dipolar + chunk_size * n_dipolar;
  external_fields = virtual_sites + chunk_size * n_vs;

  // Set all observables to zero
  for (int i = 0; i < total; i++)
    data[i] = 0.0;

  init_status = 0;
}

void Observable_stat_non_bonded::realloc_and_clear_non_bonded(size_t c_size) {
  chunk_size_nb = c_size;
  auto const n_non_bonded = max_non_bonded_pairs();
  size_t const total = chunk_size_nb * 2 * n_non_bonded;

  data_nb.resize(total);
  non_bonded_intra = data_nb.data();
  non_bonded_inter = non_bonded_intra + chunk_size_nb * n_non_bonded;

  for (int i = 0; i < total; i++)
    data_nb[i] = 0.0;
}

void Observable_stat::reduce(double *array) const {
  MPI_Reduce(data.data(), array, data.size(), MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);
}

void Observable_stat_non_bonded::reduce(double *array) const {
  MPI_Reduce(data_nb.data(), array, data_nb.size(), MPI_DOUBLE, MPI_SUM, 0,
             comm_cart);
}
