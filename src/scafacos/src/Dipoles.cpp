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

#include "scafacos/Dipoles.hpp"

#include "utils.hpp"

#include <fcs.h>

#include <cassert>
#include <string>
#include <utility>
#include <vector>

#ifdef FCS_ENABLE_DIPOLES

namespace Scafacos {

Dipoles::Dipoles(MPI_Comm comm, std::string method, std::string parameters)
    : Scafacos{comm, std::move(method), std::move(parameters)} {}

void Dipoles::set_runtime_parameters(double const *const box_l,
                                     int const *const periodicity,
                                     int const total_particles) {
  // magnetostatics: ScaFaCoS calculates near field
  auto const near_field_flag = fcs_int{1};
  Scafacos::set_runtime_parameters(box_l, periodicity, total_particles,
                                   near_field_flag);
  handle_error(fcs_set_total_dipole_particles(m_handle, total_particles));
}

void Dipoles::run(std::vector<double> &dipoles, std::vector<double> &positions,
                  std::vector<double> &fields,
                  std::vector<double> &potentials) {

  assert(dipoles.size() % 3ul == 0ul);

  fields.resize(2ul * dipoles.size());
  potentials.resize(dipoles.size());

  auto const n_part = static_cast<int>(dipoles.size() / 3ul);
  handle_error(fcs_set_dipole_particles(m_handle, n_part, positions.data(),
                                        dipoles.data(), fields.data(),
                                        potentials.data()));
  handle_error(fcs_run(m_handle, 0, nullptr, nullptr, nullptr, nullptr));
}

} // namespace Scafacos

#endif // FCS_ENABLE_DIPOLES
