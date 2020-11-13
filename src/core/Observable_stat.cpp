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

#include <utils/Span.hpp>

#include <cstddef>

Observable_stat::Observable_stat(size_t chunk_size) : m_chunk_size(chunk_size) {
  // number of chunks for different interaction types
  auto constexpr n_coulomb = 2;
  auto constexpr n_dipolar = 2;
#ifdef VIRTUAL_SITES
  auto constexpr n_vs = 1;
#else
  auto constexpr n_vs = 0;
#endif
  auto const n_bonded = bonded_ia_params.size();
  auto const n_non_bonded = max_non_bonded_pairs();
  constexpr size_t n_ext_fields = 1; // energies from all fields: accumulated
  constexpr size_t n_kinetic = 1; // linear+angular kinetic energy: accumulated

  // resize vector
  auto const total = n_kinetic + n_bonded + 2 * n_non_bonded + n_coulomb +
                     n_dipolar + n_vs + n_ext_fields;
  m_data = std::vector<double>(m_chunk_size * total);

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
