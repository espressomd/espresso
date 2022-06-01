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

#ifndef ESPRESSO_SRC_SCAFACOS_DIPOLES_HPP
#define ESPRESSO_SRC_SCAFACOS_DIPOLES_HPP

#include "scafacos/Scafacos.hpp"

#include <fcs.h>

#include <mpi.h>

#include <string>
#include <vector>

#ifdef FCS_ENABLE_DIPOLES

namespace Scafacos {

/** @brief Abstraction of a dipolar method from the ScaFaCoS library. */
struct Dipoles : public Scafacos {
  Dipoles(MPI_Comm comm, std::string method, std::string parameters);

  /** @brief Set box geometry and number of particles. */
  void set_runtime_parameters(double const *box_l, int const *periodicity,
                              int total_particles);

  /** @brief Calculate the fields and potentials for dipoles. */
  void run(std::vector<double> &dipoles, std::vector<double> &positions,
           std::vector<double> &fields, std::vector<double> &potentials);
};

} // namespace Scafacos

#endif // FCS_ENABLE_DIPOLES

#endif
