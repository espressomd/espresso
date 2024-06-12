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

#include "communication.hpp"

#include <utils/index.hpp>

#include <boost/mpi/collectives/reduce.hpp>

#include <cassert>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

Observable_stat::Observable_stat(std::size_t chunk_size, std::size_t n_bonded,
                                 int max_type)
    : m_data{}, m_chunk_size{chunk_size} {
  // number of chunks for different interaction types
  constexpr std::size_t n_coulomb = 2;
  constexpr std::size_t n_dipolar = 2;
#ifdef VIRTUAL_SITES
  constexpr std::size_t n_vs = 1;
#else
  constexpr std::size_t n_vs = 0;
#endif
  auto const n_non_bonded = get_non_bonded_offset(max_type, max_type) + 1ul;
  constexpr std::size_t n_ext_fields = 1; // reduction over all fields
  constexpr std::size_t n_kinetic = 1; // linear+angular kinetic contributions

  auto const n_elements = n_kinetic + n_bonded + 2ul * n_non_bonded +
                          n_coulomb + n_dipolar + n_vs + n_ext_fields;
  m_data = std::vector<double>(m_chunk_size * n_elements);

  // spans for the different contributions
  kinetic = std::span<double>(m_data.data(), m_chunk_size);
  bonded = std::span<double>(kinetic.end(), n_bonded * m_chunk_size);
  coulomb = std::span<double>(bonded.end(), n_coulomb * m_chunk_size);
  dipolar = std::span<double>(coulomb.end(), n_dipolar * m_chunk_size);
  virtual_sites = std::span<double>(dipolar.end(), n_vs * m_chunk_size);
  external_fields =
      std::span<double>(virtual_sites.end(), n_ext_fields * m_chunk_size);
  non_bonded_intra =
      std::span<double>(external_fields.end(), n_non_bonded * m_chunk_size);
  non_bonded_inter =
      std::span<double>(non_bonded_intra.end(), n_non_bonded * m_chunk_size);
  assert(&*non_bonded_inter.end() == (m_data.data() + m_data.size()));
}

std::size_t Observable_stat::get_non_bonded_offset(int type1, int type2) const {
  return static_cast<std::size_t>(
      Utils::lower_triangular(std::max(type1, type2), std::min(type1, type2)));
}

void Observable_stat::mpi_reduce() {
  if (comm_cart.rank() == 0) {
    std::vector<double> temp(m_data);
    boost::mpi::reduce(comm_cart, temp, m_data, std::plus<>{}, 0);
  } else {
    boost::mpi::reduce(comm_cart, m_data, std::plus<>{}, 0);
  }
}
