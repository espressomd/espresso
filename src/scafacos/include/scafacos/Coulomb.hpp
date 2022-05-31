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

#ifndef ESPRESSO_SRC_SCAFACOS_COULOMB_HPP
#define ESPRESSO_SRC_SCAFACOS_COULOMB_HPP

#include "scafacos/Scafacos.hpp"

#include <fcs.h>

#include <mpi.h>

#include <string>
#include <vector>

namespace Scafacos {

/** @brief Abstraction of a Coulomb method from the ScaFaCoS library. */
struct Coulomb : public Scafacos {
  Coulomb(MPI_Comm comm, std::string method, std::string parameters);

  /** @brief Set box geometry and number of particles. */
  void set_runtime_parameters(double const *box_l, int const *periodicity,
                              int total_particles);

  /** @brief Calculate the fields and potentials for charges. */
  void run(std::vector<double> &charges, std::vector<double> &positions,
           std::vector<double> &fields, std::vector<double> &potentials);

  /**
   * @brief Delegate the short-range calculation.
   * By default, ESPResSo calculates the short-range forces and energies
   * if the chosen ScaFaCoS method support delegation. This decision can
   * be overriden to obtain the result from a full ScaFaCos calculation.
   * @param delegate  Delegate short-range calculation to ESPResSo if true,
   *                  or to ScaFaCoS if false.
   */
  void set_near_field_delegation(bool delegate);

  bool get_near_field_delegation() const { return m_delegate_near_field; }

  /** @brief Calculate short-range pair force. */
  double pair_force(double dist) const {
    if (m_delegate_near_field) {
      fcs_float field;
      fcs_compute_near_field(m_handle, dist, &field);
      return field;
    }

    return 0.0;
  }

  /** @brief Calculate short-range pair energy. */
  double pair_energy(double dist) const {
    if (m_delegate_near_field) {
      fcs_float potential;
      fcs_compute_near_potential(m_handle, dist, &potential);
      return potential;
    }

    return 0.0;
  }

  /** @brief Tune parameters */
  void tune(std::vector<double> &charges, std::vector<double> &positions);
  /** @brief Get short-range cutoff (0.0 if not supported by the method). */
  double r_cut() const;
  /** @brief Set short-range cutoff */
  void set_r_cut(double r_cut);

private:
  auto get_near_field_flag() const {
    return static_cast<fcs_int>((m_delegate_near_field) ? 0 : 1);
  }

  /** Whether the method supports delegating the short-range calculation. */
  bool m_method_can_delegate_near_field;
  /** Whether to delegate short-range calculation to ESPResSo. */
  bool m_delegate_near_field;
};

} // namespace Scafacos
#endif
