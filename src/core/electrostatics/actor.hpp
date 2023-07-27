/*
 * Copyright (C) 2022 The ESPResSo project
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
#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_ACTOR_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_ACTOR_HPP

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/solver.hpp"

#include <stdexcept>

namespace Coulomb {

template <typename Class> class Actor {
public:
  static auto constexpr charge_neutrality_tolerance_default = 2e-12;
  /**
   * @brief Electrostatics prefactor.
   */
  double prefactor = 0.;
  /**
   * @brief Relative tolerance for the charge excess during neutrality checks.
   * To deactivate neutrality checks, set this value to -1.
   */
  double charge_neutrality_tolerance = charge_neutrality_tolerance_default;

  void set_prefactor(double new_prefactor) {
    if (new_prefactor <= 0.) {
      throw std::domain_error("Parameter 'prefactor' must be > 0");
    }
    prefactor = new_prefactor;
  }

  void sanity_checks_charge_neutrality() const {
    if (charge_neutrality_tolerance != -1.) {
      check_charge_neutrality(charge_neutrality_tolerance);
    }
  }
};

} // namespace Coulomb

#endif // ELECTROSTATICS
#endif
