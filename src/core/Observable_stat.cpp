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

#include "config.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/mpi/communicator.hpp>

extern boost::mpi::communicator comm_cart;

/** Tracker of observables */
std::vector<Observable_stat *> &registered_obs_stat() {
  static std::vector<Observable_stat *> s_registered_obs_stat;
  return s_registered_obs_stat;
}

/** Tracker of non-bonded observables */
std::vector<Observable_stat_non_bonded *> &registered_obs_stat_non_bonded() {
  static std::vector<Observable_stat_non_bonded *>
      s_registered_obs_stat_non_bonded;
  return s_registered_obs_stat_non_bonded;
}

void Observable_stat::register_obs() { registered_obs_stat().push_back(this); }

void Observable_stat_non_bonded::register_obs() {
  registered_obs_stat_non_bonded().push_back(this);
}

void invalidate_obs() {
  for (auto *obs : registered_obs_stat())
    obs->is_initialized = false;
}

void realloc_and_clear_all_obs() {
  for (auto *obs : registered_obs_stat())
    obs->realloc_and_clear();
  for (auto *obs : registered_obs_stat_non_bonded())
    obs->realloc_and_clear();
}

void Observable_stat::realloc_and_clear() {
  // number of chunks for different interaction types
  auto const n_coulomb =
      m_pressure_obs ? Coulomb::pressure_n() : Coulomb::energy_n();
  auto const n_dipolar =
      m_pressure_obs ? Dipole::pressure_n() : Dipole::energy_n();
  size_t n_vs = 0;
#ifdef VIRTUAL_SITES
  n_vs = 1;
#endif
  auto const n_bonded = bonded_ia_params.size();
  auto const n_non_bonded = max_non_bonded_pairs();
  constexpr size_t n_ext_fields = 1; // energies from all fields: accumulated
  constexpr size_t n_kinetic = 1; // linear+angular kinetic energy: accumulated

  // reallocate and reset memory
  auto const total = n_kinetic + n_bonded + n_non_bonded + n_coulomb +
                     n_dipolar + n_vs + n_ext_fields;
  Observable_stat_base::realloc_and_clear(total);
  is_initialized = false;

  // spans for the different contributions
  kinetic = Utils::Span<double>(data.data(), m_chunk_size);
  bonded = Utils::Span<double>(kinetic.end(), n_bonded * m_chunk_size);
  non_bonded = Utils::Span<double>(bonded.end(), n_non_bonded * m_chunk_size);
  coulomb = Utils::Span<double>(non_bonded.end(), n_coulomb * m_chunk_size);
  dipolar = Utils::Span<double>(coulomb.end(), n_dipolar * m_chunk_size);
  virtual_sites = Utils::Span<double>(dipolar.end(), n_vs * m_chunk_size);
  external_fields =
      Utils::Span<double>(virtual_sites.end(), n_ext_fields * m_chunk_size);
}

void Observable_stat_non_bonded::realloc_and_clear() {
  // number of chunks for different interaction types
  auto const n_non_bonded = max_non_bonded_pairs();
  // reallocate and reset memory
  Observable_stat_base::realloc_and_clear(2 * n_non_bonded);
  // spans for the different contributions
  auto const span_size = n_non_bonded * m_chunk_size;
  non_bonded_intra = Utils::Span<double>(data.data(), span_size);
  non_bonded_inter = Utils::Span<double>(non_bonded_intra.end(), span_size);
}

void Observable_stat_base::reduce(double *out) const {
  MPI_Reduce(data.data(), out, data.size(), MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}
