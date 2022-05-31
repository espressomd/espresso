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

#include "scafacos/Coulomb.hpp"

#include "utils.hpp"

#include <fcs.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Scafacos {

Coulomb::Coulomb(MPI_Comm comm, std::string method, std::string parameters)
    : Scafacos{comm, std::move(method), std::move(parameters)} {
  fcs_int near_field_delegation;
  fcs_get_near_field_delegation(m_handle, &near_field_delegation);
  m_method_can_delegate_near_field = static_cast<bool>(near_field_delegation);
  m_delegate_near_field = m_method_can_delegate_near_field;
}

void Coulomb::set_runtime_parameters(double const *const box_l,
                                     int const *const periodicity,
                                     int const total_particles) {
  auto const near_field_flag = get_near_field_flag();
  Scafacos::set_runtime_parameters(box_l, periodicity, total_particles,
                                   near_field_flag);
}

void Coulomb::set_near_field_delegation(bool delegate) {
  if (delegate != m_delegate_near_field) {
    if (delegate and not m_method_can_delegate_near_field) {
      throw std::runtime_error("Method '" + get_method() +
                               "' cannot delegate short-range calculation");
    }
    m_delegate_near_field = delegate;
    auto const near_field_flag = get_near_field_flag();
    handle_error(fcs_set_near_field_flag(m_handle, near_field_flag));
  }
}

double Coulomb::r_cut() const {
  if (m_delegate_near_field) {
    fcs_float r_cut;
    fcs_get_r_cut(m_handle, &r_cut);
    return r_cut;
  }
  return 0.0;
}

void Coulomb::set_r_cut(double r_cut) {
  if (m_delegate_near_field) {
    fcs_set_r_cut(m_handle, r_cut);
  }
}

void Coulomb::run(std::vector<double> &charges, std::vector<double> &positions,
                  std::vector<double> &fields,
                  std::vector<double> &potentials) {

  auto const n_part = charges.size();
  fields.resize(3ul * n_part);
  potentials.resize(n_part);

  tune(charges, positions);
  auto const size = static_cast<int>(n_part);
  handle_error(fcs_run(m_handle, size, positions.data(), charges.data(),
                       fields.data(), potentials.data()));
}

void Coulomb::tune(std::vector<double> &charges,
                   std::vector<double> &positions) {
  auto const n_part = charges.size();
  auto const size = static_cast<int>(n_part);
  handle_error(fcs_tune(m_handle, size, positions.data(), charges.data()));
}

} // namespace Scafacos
