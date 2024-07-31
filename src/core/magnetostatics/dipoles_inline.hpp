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

#include "magnetostatics/dipoles.hpp"
#include "magnetostatics/dp3m.hpp"
#include "magnetostatics/solver.hpp"

#include "Particle.hpp"

#include "actor/traits.hpp"
#include "actor/visitors.hpp"

#include <utils/Vector.hpp>

#include <functional>
#include <optional>
#include <variant>

namespace Dipoles {

struct ShortRangeForceKernel {

  using kernel_type = Solver::ShortRangeForceKernel;
  using result_type = std::optional<kernel_type>;

#ifdef DIPOLES
  template <typename T>
  result_type operator()(std::shared_ptr<T> const &) const {
    return {};
  }

#ifdef DP3M
  result_type operator()(std::shared_ptr<DipolarP3M> const &ptr) const {
    auto const &actor = *ptr;
    return kernel_type{[&actor](Particle const &p1, Particle const &p2,
                                Utils::Vector3d const &d, double dist,
                                double dist2) {
      return actor.pair_force(p1, p2, d, dist2, dist);
    }};
  }
#endif // DP3M

  result_type
  operator()(std::shared_ptr<DipolarLayerCorrection> const &ptr) const {
    return std::visit(*this, ptr->base_solver);
  }
#endif // DIPOLES
};

struct ShortRangeEnergyKernel {

  using kernel_type = Solver::ShortRangeEnergyKernel;
  using result_type = std::optional<kernel_type>;

#ifdef DIPOLES
  template <typename T>
  result_type operator()(std::shared_ptr<T> const &) const {
    return {};
  }

#ifdef DP3M
  result_type operator()(std::shared_ptr<DipolarP3M> const &ptr) const {
    auto const &actor = *ptr;
    return kernel_type{[&actor](Particle const &p1, Particle const &p2,
                                Utils::Vector3d const &d, double dist,
                                double dist2) {
      return actor.pair_energy(p1, p2, d, dist2, dist);
    }};
  }
#endif // DP3M

  result_type
  operator()(std::shared_ptr<DipolarLayerCorrection> const &ptr) const {
    return std::visit(*this, ptr->base_solver);
  }
#endif // DIPOLES
};

inline std::optional<Solver::ShortRangeForceKernel>
Solver::pair_force_kernel() const {
#ifdef DIPOLES
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Dipoles::ShortRangeForceKernel();
    return std::visit(visitor, *solver);
  }
#endif // DIPOLES
  return std::nullopt;
}

inline std::optional<Solver::ShortRangeEnergyKernel>
Solver::pair_energy_kernel() const {
#ifdef DIPOLES
  if (auto &solver = impl->solver; solver.has_value()) {
    auto const visitor = Dipoles::ShortRangeEnergyKernel();
    return std::visit(visitor, *solver);
  }
#endif // DIPOLES
  return std::nullopt;
}

} // namespace Dipoles
