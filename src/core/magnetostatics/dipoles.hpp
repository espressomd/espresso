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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLES_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLES_HPP

#include "config/config.hpp"

#ifdef DIPOLES

#include "actor/traits.hpp"

#include "magnetostatics/solver.hpp"

#include "magnetostatics/barnes_hut_gpu.hpp"
#include "magnetostatics/dipolar_direct_sum.hpp"
#include "magnetostatics/dipolar_direct_sum_gpu.hpp"
#include "magnetostatics/dlc.hpp"
#include "magnetostatics/dp3m.hpp"
#include "magnetostatics/scafacos.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>

/** Get the magnetostatics prefactor. */
struct GetDipolesPrefactor {
  template <typename T>
  double operator()(std::shared_ptr<T> const &actor) const {
    return actor->prefactor;
  }
};

/** Run actor sanity checks. */
struct DipolesSanityChecks {
  template <typename T> void operator()(std::shared_ptr<T> const &actor) const {
    actor->sanity_checks();
  }
};

namespace Dipoles {
namespace traits {

/** @brief Whether an actor is a solver. */
template <typename T>
using is_solver = std::is_convertible<std::shared_ptr<T>, MagnetostaticsActor>;

/** @brief The dipolar method supports dipole fields calculation. */
template <class T> struct has_dipole_fields : std::false_type {};
#ifdef DIPOLE_FIELD_TRACKING
template <> struct has_dipole_fields<DipolarDirectSum> : std::true_type {};
#endif // DIPOLE_FIELD_TRACKING

} // namespace traits
} // namespace Dipoles

#endif // DIPOLES
#endif
