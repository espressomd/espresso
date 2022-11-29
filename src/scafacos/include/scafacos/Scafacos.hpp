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

#ifndef ESPRESSO_SRC_SCAFACOS_SCAFACOS_HPP
#define ESPRESSO_SRC_SCAFACOS_SCAFACOS_HPP

#include <fcs.h>

#include <mpi.h>

#include <string>
#include <type_traits>
#include <vector>

namespace Scafacos {

/** @brief Abstraction of a method from the ScaFaCoS library. */
struct Scafacos {
  Scafacos(MPI_Comm comm, std::string method, std::string parameters);
  ~Scafacos();
  /** Get the parameters from the library */
  std::string get_parameters() const { return m_parameters; }
  /** Get active method name */
  std::string get_method() const { return m_method_name; }

  /** Set box geometry, number of particles and calculation type. */
  void set_runtime_parameters(double const *box_l, int const *periodicity,
                              int total_particles, int near_field_flag);

  /** Get a list of methods supported by the library */
  static std::vector<std::string> available_methods();

protected:
  /** Handle from the library */
  FCS m_handle;

private:
  /** The method name */
  std::string m_method_name;
  /** The parameters set */
  std::string m_parameters;
};

static_assert(std::is_same_v<fcs_int, int>,
              "ScaFaCoS must be compiled with fcs_int = int");
static_assert(std::is_same_v<fcs_float, double>,
              "ScaFaCoS must be compiled with fcs_float = double");

} // namespace Scafacos
#endif
