/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <algorithm>
#include <chrono>
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

  std::cout << "RUN: allocating copies, size = " << n_part << std::endl;
  std::vector<double> positions2(3ul * n_part);
  // positions2.reserve(positions.size());
  std::vector<double> charges2(n_part);
  // charges2.reserve(charges.size());
  std::vector<double> fields2(3ul * n_part);
  // fields2.reserve(fields.size());
  std::vector<double> potentials2(n_part);
  // potentials2.reserve(potentials.size());

  std::cout << "RUN: copying data" << std::endl;
  std::copy(positions.begin(), positions.end(), positions2.begin());
  std::copy(charges.begin(), charges.end(), charges2.begin());
  std::copy(fields.begin(), fields.end(), fields2.begin());
  std::copy(potentials.begin(), potentials.end(), potentials2.begin());

  std::cout << "RUN: MPI parameters" << std::endl;
  int rank;
  MPI_Comm_rank(exData.globalComm(), &rank);
  int nranks;
  MPI_Comm_size(exData.globalComm(), &nranks);
  std::cout << "RUN: assigning sizes for partitions" << std::endl;
  exObj.setParameter(ExscalicosP3MShortRangeProcess(), exData, 0, nranks - 1);
  exObj.setParameter(ExscalicosP3MLongRangeProcess(), exData, 0, nranks - 1);

  // std::cout << rank << " running ExScaLiCoS::tune now - using " << nranks
  //           << " ranks : " << std::endl;
  if (setupExscalicos != 0) {
    // auto rMax = r_cut();
    // exObj.setRMax(rMax);
    exObj.tune(exData, positions2, charges2);
    setupExscalicos = 0;
  }

  std::cout << "RUN: running ScaFaCoS" << std::endl;
  std::cout << rank << " running fcs_run now: " << std::endl;
  // std::cout << "RUN: checking charges / positions content:" << std::endl;
  // for (int i = 0; i < std::min(10, (int)charges.size()); ++i)
  //   std::cout << rank << " " << i << " " << positions[3 * i + 0] << " "
  //             << positions[3 * i + 1] << " " << positions[3 * i + 2] << " "
  //             << charges[i] << std::endl;

  auto const size = static_cast<int>(n_part);

  handle_error(fcs_run(m_handle, size, positions.data(), charges.data(),
                       fields.data(), potentials.data()));
  MPI_Barrier(exData.globalComm());
  const auto startScafacos{std::chrono::steady_clock::now()};
  handle_error(fcs_run(m_handle, size, positions.data(), charges.data(),
                       fields.data(), potentials.data()));
  const auto endScafacos{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> timeScafacos{endScafacos - startScafacos};
  MPI_Barrier(exData.globalComm());

  // std::cout << "RUN: checking charges2 / positions2 content:" << std::endl;
  // for (int i = 0; i < std::min(10, (int)charges2.size()); ++i)
  //   std::cout << rank << " " << i << " " << positions2[3 * i + 0] << " "
  //             << positions2[3 * i + 1] << " " << positions2[3 * i + 2] << "
  //             "
  //             << charges2[i] << std::endl;

  std::cout << "RUN: running ExScaLiCoS" << std::endl;
  std::cout << rank << " running ExScaLiCoS::run now: " << std::endl;
  MPI_Barrier(exData.globalComm());
  const auto startExscalicos{std::chrono::steady_clock::now()};
  exObj.run(exData, positions2, charges2, potentials2, fields2);
  const auto endExscalicos{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> timeExscalicos{endExscalicos -
                                                     startExscalicos};
  MPI_Barrier(exData.globalComm());

  double tScafacos = timeScafacos.count();
  double tExscalicos = timeExscalicos.count();

  MPI_Allreduce(&tScafacos, MPI_IN_PLACE, 1, MPI_DOUBLE, MPI_MAX,
                exData.globalComm());
  MPI_Allreduce(&tExscalicos, MPI_IN_PLACE, 1, MPI_DOUBLE, MPI_MAX,
                exData.globalComm());

  std::vector<double> diff = {0.0, 0.0, 0.0, 0.0, 0.0};
  for (auto i = 0; i < charges2.size(); ++i) {
    diff[0] +=
        (fields[3 * i] - fields2[3 * i]) * (fields[3 * i] - fields2[3 * i]);
    diff[1] += (fields[3 * i + 1] - fields2[3 * i + 1]) *
               (fields[3 * i + 1] - fields2[3 * i + 1]);
    diff[2] += (fields[3 * i + 2] - fields2[3 * i + 2]) *
               (fields[3 * i + 2] - fields2[3 * i + 2]);
    diff[3] +=
        (potentials[i] - potentials2[i]) * (potentials[i] - potentials2[i]);
  }
  diff[4] = (double)(charges.size());

  MPI_Allreduce(diff.data(), MPI_IN_PLACE, 5, MPI_DOUBLE, MPI_SUM,
                exData.globalComm());

  diff[0] /= diff[4];
  diff[1] /= diff[4];
  diff[2] /= diff[4];
  diff[3] /= diff[4];

  diff[0] = std::sqrt(diff[0]);
  diff[1] = std::sqrt(diff[1]);
  diff[2] = std::sqrt(diff[2]);
  diff[3] = std::sqrt(diff[3]);

  if (rank == 0) {
    std::cout << "Comparision ScaFaCoS - ExScaLiCoS:\n deltaP = " << diff[3]
              << '\n'
              << "deltaf = " << diff[0] << " " << diff[1] << " " << diff[2]
              << '\n'
              << "tScaFaCoS: " << tScafacos << "\n"
              << "tExScaLiCoS: " << tExscalicos << "\n"
              << std::endl;
  }

  /*
  if (rank == 0) {
    std::cout << "Comparison results (fields / potentials): " << std::endl;
    for (int i = 0; i < potentials.size(); ++i) {
      std::cout << " ( " << potentials[i] << " , " << potentials2[i] << " ) "
                << " ( " << fields[3 * i + 0] << " , " << fields2[3 * i + 0]
                << " ) "
                << " ( " << fields[3 * i + 1] << " , " << fields2[3 * i + 1]
                << " ) "
                << " ( " << fields[3 * i + 2] << " , " << fields2[3 * i + 2]
                << " ) " << std::endl;
    }
  }
  */

  fields = fields2;
  potentials = potentials2;
}

void Coulomb::tune(std::vector<double> &charges,
                   std::vector<double> &positions) {
  auto const n_part = charges.size();
  auto const size = static_cast<int>(n_part);

  handle_error(fcs_tune(m_handle, size, positions.data(), charges.data()));

  std::cout << "RUN: allocating copies, size = " << n_part << std::endl;
  std::vector<double> positions2(3ul * n_part);
  // positions2.reserve(positions.size());
  std::vector<double> charges2(n_part);
  // charges2.reserve(charges.size());

  std::cout << "RUN: copying data" << std::endl;
  std::copy(positions.begin(), positions.end(), positions2.begin());
  std::copy(charges.begin(), charges.end(), charges2.begin());

  std::cout << "RUN: MPI parameters" << std::endl;
  int rank;
  MPI_Comm_rank(exData.globalComm(), &rank);
  int nranks;
  MPI_Comm_size(exData.globalComm(), &nranks);
  std::cout << "RUN: assigning sizes for partitions" << std::endl;
  exObj.setParameter(ExscalicosP3MShortRangeProcess(), exData, 0, nranks - 1);
  exObj.setParameter(ExscalicosP3MLongRangeProcess(), exData, 0, nranks - 1);

  // exObj.tune(exData, charges2, positions2);

  std::cout << "FINISHED tuning" << std::endl;
}

} // namespace Scafacos
