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

#include "electrostatics/coulomb.hpp"
#include "electrostatics/solver.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>
#include <utils/demangle.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/matrix.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <tuple>
#include <variant>

namespace Coulomb {

struct ShortRangeForceKernel {

  using kernel_type = Solver::ShortRangeForceKernel;
  using result_type = std::optional<kernel_type>;

#ifdef ESPRESSO_ELECTROSTATICS
  template <typename T>
  result_type operator()(std::shared_ptr<T> const &ptr) const {
    auto const &actor = *ptr;
    return kernel_type{
        [&actor](double q1q2, Utils::Vector3d const &d, double dist) {
          return actor.pair_force(q1q2, d, dist);
        }};
  }

#ifdef ESPRESSO_P3M
  auto
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &ptr) const {
    return std::visit(*this, ptr->base_solver);
  }
#endif // ESPRESSO_P3M
#endif // ESPRESSO_ELECTROSTATICS
};

struct ShortRangeForceCorrectionsKernel {

  using kernel_type = Solver::ShortRangeForceCorrectionsKernel;
  using result_type = std::optional<kernel_type>;

  template <typename T>
  result_type operator()(std::shared_ptr<T> const &) const {
    return {};
  }

#ifdef ESPRESSO_P3M
  result_type
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &ptr) const {
    auto const &actor = *ptr;
    return kernel_type{[&actor](Utils::Vector3d const &pos1,
                                Utils::Vector3d const &pos2,
                                Utils::Vector3d &p1f_asym,
                                Utils::Vector3d &p2f_asym, double q1q2) {
      actor.add_pair_force_corrections(pos1, pos2, p1f_asym, p2f_asym, q1q2);
    }};
  }
#endif // ESPRESSO_P3M
};

struct ShortRangePressureKernel {

  using kernel_type = Solver::ShortRangePressureKernel;
  using result_type = std::optional<kernel_type>;

#ifdef ESPRESSO_ELECTROSTATICS
  template <typename T>
  result_type operator()(std::shared_ptr<T> const &ptr) const {
    if constexpr (traits::has_pressure<T>::value) {
      return kernel_type{
          [&actor = *ptr](double q1q2, Utils::Vector3d const &d, double dist) {
            return Utils::tensor_product(actor.pair_force(q1q2, d, dist), d);
          }};
    }
    return {};
  }
#endif // ESPRESSO_ELECTROSTATICS
};

struct ShortRangeEnergyKernel {

  using kernel_type = Solver::ShortRangeEnergyKernel;
  using result_type = std::optional<kernel_type>;

#ifdef ESPRESSO_ELECTROSTATICS
  template <typename T>
  result_type operator()(std::shared_ptr<T> const &ptr) const {
    auto const &actor = *ptr;
    return kernel_type{
        [&actor](Utils::Vector3d const &, Utils::Vector3d const &, double q1q2,
                 Utils::Vector3d const &,
                 double dist) { return actor.pair_energy(q1q2, dist); }};
  }
#ifdef ESPRESSO_P3M
  result_type
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &ptr) const {
    auto const &actor = *ptr;
    auto const energy_kernel = std::visit(*this, actor.base_solver);
    return kernel_type{
        [&actor, energy_kernel](Utils::Vector3d const &pos1,
                                Utils::Vector3d const &pos2, double q1q2,
                                Utils::Vector3d const &d, double dist) {
          auto energy = 0.;
          if (energy_kernel) {
            energy = (*energy_kernel)(pos1, pos2, q1q2, d, dist);
          }
          return energy + actor.pair_energy_correction(pos1, pos2, q1q2);
        }};
  }
#endif // ESPRESSO_P3M
#ifdef ESPRESSO_GSL
  result_type operator()(std::shared_ptr<CoulombMMM1D> const &actor) const {
    return kernel_type{
        [&actor](Utils::Vector3d const &, Utils::Vector3d const &, double q1q2,
                 Utils::Vector3d const &d,
                 double dist) { return actor->pair_energy(q1q2, d, dist); }};
  }
#endif // ESPRESSO_GSL
#endif // ESPRESSO_ELECTROSTATICS
};

inline std::optional<Solver::ShortRangeForceKernel>
Solver::pair_force_kernel() const {
#ifdef ESPRESSO_ELECTROSTATICS
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Coulomb::ShortRangeForceKernel();
    return std::visit(visitor, *solver);
  }
#endif // ESPRESSO_ELECTROSTATICS
  return std::nullopt;
}

inline std::optional<Solver::ShortRangeForceCorrectionsKernel>
Solver::pair_force_elc_kernel() const {
#ifdef ESPRESSO_ELECTROSTATICS
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Coulomb::ShortRangeForceCorrectionsKernel();
    return std::visit(visitor, *solver);
  }
#endif // ESPRESSO_ELECTROSTATICS
  return std::nullopt;
}

inline std::optional<Solver::ShortRangePressureKernel>
Solver::pair_pressure_kernel() const {
#ifdef ESPRESSO_ELECTROSTATICS
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Coulomb::ShortRangePressureKernel();
    return std::visit(visitor, *solver);
  }
#endif // ESPRESSO_ELECTROSTATICS
  return std::nullopt;
}

inline std::optional<Solver::ShortRangeEnergyKernel>
Solver::pair_energy_kernel() const {
#ifdef ESPRESSO_ELECTROSTATICS
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Coulomb::ShortRangeEnergyKernel();
    return std::visit(visitor, *solver);
  }
#endif // ESPRESSO_ELECTROSTATICS
  return std::nullopt;
}

} // namespace Coulomb
