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

#include "scafacos/Scafacos.hpp"

#include "utils.hpp"

#include <string>
#include <utility>
#include <vector>

namespace Scafacos {

std::vector<std::string> Scafacos::available_methods() {
  std::vector<std::string> methods;

#ifdef FCS_ENABLE_DIRECT
  methods.emplace_back("direct");
#endif
#ifdef FCS_ENABLE_EWALD
  methods.emplace_back("ewald");
#endif
#ifdef FCS_ENABLE_FMM
  methods.emplace_back("fmm");
#endif
#ifdef FCS_ENABLE_MEMD
  methods.emplace_back("memd");
#endif
#ifdef FCS_ENABLE_MMM1D
  methods.emplace_back("mmm1d");
#endif
#ifdef FCS_ENABLE_MMM2D
  methods.emplace_back("mmm2d");
#endif
#ifdef FCS_ENABLE_P2NFFT
  methods.emplace_back("p2nfft");
#endif
#ifdef FCS_ENABLE_P3M
  methods.emplace_back("p3m");
#endif
#ifdef FCS_ENABLE_PEPC
  methods.emplace_back("pepc");
#endif
#ifdef FCS_ENABLE_PP3MG
  methods.emplace_back("pp3mg");
#endif
#ifdef FCS_ENABLE_VMG
  methods.emplace_back("vmg");
#endif
#ifdef FCS_ENABLE_WOLF
  methods.emplace_back("wolf");
#endif

  return methods;
}

Scafacos::Scafacos(MPI_Comm comm, std::string method, std::string parameters)
    : m_method_name{std::move(method)}, m_parameters{std::move(parameters)} {

  handle_error(fcs_init(&m_handle, m_method_name.c_str(), comm));

  fcs_set_resort(m_handle, 0);

  handle_error(fcs_parser(m_handle, m_parameters.c_str(), 0));
}

Scafacos::~Scafacos() { fcs_destroy(m_handle); }

void Scafacos::set_runtime_parameters(double const *box_l,
                                      int const *periodicity,
                                      int total_particles,
                                      int near_field_flag) {
  // define rectangular box
  double boxa[3] = {box_l[0], 0., 0.};
  double boxb[3] = {0., box_l[1], 0.};
  double boxc[3] = {0., 0., box_l[2]};
  double off[3] = {0., 0., 0.};
  handle_error(fcs_set_common(m_handle, near_field_flag, boxa, boxb, boxc, off,
                              periodicity, total_particles));
}

} // namespace Scafacos
