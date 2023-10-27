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

#pragma once

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
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <variant>

namespace Dipoles {

using MagnetostaticsActor =
    std::variant<std::shared_ptr<DipolarDirectSum>,
#ifdef DIPOLAR_DIRECT_SUM
                 std::shared_ptr<DipolarDirectSumGpu>,
#endif
#ifdef DIPOLAR_BARNES_HUT
                 std::shared_ptr<DipolarBarnesHutGpu>,
#endif
#ifdef DP3M
                 std::shared_ptr<DipolarP3M>,
#endif
#ifdef SCAFACOS_DIPOLES
                 std::shared_ptr<DipolarScafacos>,
#endif
                 std::shared_ptr<DipolarLayerCorrection>>;

struct Solver::Implementation {
  /// @brief Main electrostatics solver.
  std::optional<MagnetostaticsActor> solver;
  Implementation() : solver{} {}
};

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
