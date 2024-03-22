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
#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_HPP

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "actor/traits.hpp"

#include "electrostatics/solver.hpp"

#include "electrostatics/debye_hueckel.hpp"
#include "electrostatics/elc.hpp"
#include "electrostatics/icc.hpp"
#include "electrostatics/mmm1d.hpp"
#include "electrostatics/mmm1d_gpu.hpp"
#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m_gpu.hpp"
#include "electrostatics/reaction_field.hpp"
#include "electrostatics/scafacos.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <type_traits>
#include <variant>

namespace Coulomb {

using ElectrostaticsActor =
    std::variant<std::shared_ptr<DebyeHueckel>,
#ifdef P3M
                 std::shared_ptr<CoulombP3M>,
#ifdef CUDA
                 std::shared_ptr<CoulombP3MGPU>,
#endif // CUDA
                 std::shared_ptr<ElectrostaticLayerCorrection>,
#endif // P3M
                 std::shared_ptr<CoulombMMM1D>,
#ifdef MMM1D_GPU
                 std::shared_ptr<CoulombMMM1DGpu>,
#endif // MMM1D_GPU
#ifdef SCAFACOS
                 std::shared_ptr<CoulombScafacos>,
#endif // SCAFACOS
                 std::shared_ptr<ReactionField>>;

using ElectrostaticsExtension = std::variant<std::shared_ptr<ICCStar>>;

struct Solver::Implementation {
  /// @brief Main electrostatics solver.
  std::optional<ElectrostaticsActor> solver;
  /// @brief Extension that modifies the solver behavior.
  std::optional<ElectrostaticsExtension> extension;
  Implementation() : solver{}, extension{} {}
};

namespace traits {

#ifdef P3M
/** @brief Whether an actor can be adapted by ELC. */
template <typename T>
using elc_adaptable =
    std::is_convertible<std::shared_ptr<T>,
                        ElectrostaticLayerCorrection::BaseSolver>;
#endif // P3M

/** @brief Whether an actor is a solver. */
template <typename T>
using is_solver = std::is_convertible<std::shared_ptr<T>, ElectrostaticsActor>;
/** @brief Whether an actor is an extension. */
template <typename T>
using is_extension =
    std::is_convertible<std::shared_ptr<T>, ElectrostaticsExtension>;

/** @brief The electrostatic method supports pressure calculation. */
template <class T> struct has_pressure : std::true_type {};
#ifdef P3M
template <>
struct has_pressure<ElectrostaticLayerCorrection> : std::false_type {};
#endif // P3M
#ifdef MMM1D_GPU
template <> struct has_pressure<CoulombMMM1DGpu> : std::false_type {};
#endif // MMM1D_GPU
#ifdef SCAFACOS
template <> struct has_pressure<CoulombScafacos> : std::false_type {};
#endif // SCAFACOS
template <> struct has_pressure<CoulombMMM1D> : std::false_type {};

} // namespace traits
} // namespace Coulomb
#endif // ELECTROSTATICS
#endif
