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

#include "Scafacos.hpp"
#include <cassert>

namespace Scafacos {

#define handle_error(stmt)                                                     \
  {                                                                            \
    const FCSResult res = stmt;                                                \
    if (res)                                                                   \
      throw std::runtime_error(fcs_result_get_message(res));                   \
  }

std::string Scafacos::get_parameters() { return m_last_parameters; }
std::string Scafacos::get_method() { return method; }

std::list<std::string> Scafacos::available_methods() {
  std::list<std::string> methods;

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

Scafacos::Scafacos(std::string method, MPI_Comm comm,
                   const std::string &parameters)
    : method(std::move(method)) {
  handle_error(fcs_init(&handle, this->method.c_str(), comm));

  int near_flag;
  fcs_get_near_field_delegation(handle, &near_flag);
  has_near = near_flag != 0;

  fcs_set_resort(handle, 0);

  parse_parameters(parameters);
}

Scafacos::~Scafacos() { fcs_destroy(handle); }

void Scafacos::parse_parameters(const std::string &s) {
  m_last_parameters = s;
  handle_error(fcs_parser(handle, s.c_str(), 0));
}

double Scafacos::r_cut() const {
  if (has_near) {
    fcs_float r_cut;
    fcs_get_r_cut(handle, &r_cut);
    return r_cut;
  }
  return 0.0;
}

void Scafacos::set_r_cut(double r_cut) {
  if (has_near) {
    fcs_set_r_cut(handle, r_cut);
  }
}

void Scafacos::run(std::vector<double> &charges, std::vector<double> &positions,
                   std::vector<double> &fields,
                   std::vector<double> &potentials) {
  const int local_n_part = charges.size();

  fields.resize(3 * local_n_part);
  potentials.resize(local_n_part);

  handle_error(fcs_tune(handle, local_n_part, &(positions[0]), &(charges[0])));

  handle_error(fcs_run(handle, local_n_part, &(positions[0]), &(charges[0]),
                       &(fields[0]), &(potentials[0])));
}

#ifdef FCS_ENABLE_DIPOLES
void Scafacos::run_dipolar(std::vector<double> &dipoles,
                           std::vector<double> &positions,
                           std::vector<double> &fields,
                           std::vector<double> &potentials) {
  assert(dipoles.size() % 3 == 0);
  auto const local_n_part = dipoles.size() / 3;

  fields.resize(6 * local_n_part);
  potentials.resize(3 * local_n_part);

  handle_error(fcs_set_dipole_particles(handle, static_cast<int>(local_n_part),
                                        &(positions[0]), &(dipoles[0]),
                                        &(fields[0]), &(potentials[0])));
  handle_error(fcs_run(handle, 0, nullptr, nullptr, nullptr, nullptr));
}
#endif

void Scafacos::tune(std::vector<double> &charges,
                    std::vector<double> &positions) {
  handle_error(
      fcs_tune(handle, charges.size(), &(positions[0]), &(charges[0])));
}

void Scafacos::set_common_parameters(const double *box_l,
                                     const int *periodicity,
                                     int total_particles) {
  double boxa[3] = {0., 0., 0.}, boxb[3] = {0., 0., 0.}, boxc[3] = {0., 0., 0.},
         off[3] = {0., 0., 0.};
  boxa[0] = box_l[0];
  boxb[1] = box_l[1];
  boxc[2] = box_l[2];
  // Does scafacos calculate the near field part
  // For charges, if the method supports it, Es calculates near field
  int sr = 0;
  if (!dipolar() && has_near) {
    sr = 0;
  } else {
    // Scafacos does near field calc
    sr = 1;
  }
  handle_error(fcs_set_common(handle, sr, boxa, boxb, boxc, off, periodicity,
                              total_particles));
#ifdef FCS_ENABLE_DIPOLES
  if (m_dipolar)
    handle_error(fcs_set_total_dipole_particles(handle, total_particles));
#endif
}

void Scafacos::set_dipolar(bool d) {
#ifndef FCS_ENABLE_DIPOLES
  if (d) {
    throw std::runtime_error("Scafacos library does not have dipole support.");
  }
#endif
  m_dipolar = d;
}

} // namespace Scafacos
