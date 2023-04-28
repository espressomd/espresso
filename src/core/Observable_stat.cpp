/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/Span.hpp>
#include <utils/index.hpp>

#include <boost/mpi/collectives/reduce.hpp>

#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

static auto max_non_bonded_pairs() {
  auto const n = static_cast<std::size_t>(max_seen_particle_type);
  return (n * (n + 1)) / 2;
}

Observable_stat::Observable_stat(std::size_t chunk_size)
    : m_chunk_size(chunk_size) {
  // number of chunks for different interaction types
  constexpr std::size_t n_coulomb = 2;
  constexpr std::size_t n_dipolar = 2;
#ifdef VIRTUAL_SITES
  constexpr std::size_t n_vs = 1;
#else
  constexpr std::size_t n_vs = 0;
#endif
  auto const n_bonded =
      static_cast<std::size_t>(bonded_ia_params.get_next_key());
  auto const n_non_bonded = max_non_bonded_pairs();
  constexpr std::size_t n_ext_fields = 1; // reduction over all fields
  constexpr std::size_t n_kinetic = 1; // linear+angular kinetic contributions

  auto const n_elements = n_kinetic + n_bonded + 2 * n_non_bonded + n_coulomb +
                          n_dipolar + n_vs + n_ext_fields;
  m_data = std::vector<double>(m_chunk_size * n_elements);

  // spans for the different contributions
  kinetic = Utils::Span<double>(m_data.data(), m_chunk_size);
  bonded = Utils::Span<double>(kinetic.end(), n_bonded * m_chunk_size);
  coulomb = Utils::Span<double>(bonded.end(), n_coulomb * m_chunk_size);
  dipolar = Utils::Span<double>(coulomb.end(), n_dipolar * m_chunk_size);
  virtual_sites = Utils::Span<double>(dipolar.end(), n_vs * m_chunk_size);
  external_fields =
      Utils::Span<double>(virtual_sites.end(), n_ext_fields * m_chunk_size);
  non_bonded_intra =
      Utils::Span<double>(external_fields.end(), n_non_bonded * m_chunk_size);
  non_bonded_inter =
      Utils::Span<double>(non_bonded_intra.end(), n_non_bonded * m_chunk_size);
  assert(non_bonded_inter.end() == (m_data.data() + m_data.size()));
}

Utils::Span<double>
Observable_stat::non_bonded_contribution(Utils::Span<double> base_pointer,
                                         int type1, int type2) const {
  auto const offset = static_cast<std::size_t>(Utils::upper_triangular(
      std::min(type1, type2), std::max(type1, type2), max_seen_particle_type));
  return {base_pointer.begin() + offset * m_chunk_size, m_chunk_size};
}

void Observable_stat::mpi_reduce() {
  if (comm_cart.rank() == 0) {
    std::vector<double> temp(m_data);
    boost::mpi::reduce(comm_cart, temp, m_data, std::plus<>{}, 0);
  } else {
    boost::mpi::reduce(comm_cart, m_data, std::plus<>{}, 0);
  }
}
