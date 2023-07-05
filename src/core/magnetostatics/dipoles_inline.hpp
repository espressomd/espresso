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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLES_INLINE_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_DIPOLES_INLINE_HPP

#include "config/config.hpp"

#include "magnetostatics/dipoles.hpp"
#include "magnetostatics/dp3m.hpp"
#include "magnetostatics/solver.hpp"

#include "Particle.hpp"

#include "actor/traits.hpp"
#include "actor/visitors.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <functional>

namespace Dipoles {

struct ShortRangeForceKernel
    : public boost::static_visitor<
          boost::optional<Solver::ShortRangeForceKernel>> {

  using kernel_type = result_type::value_type;

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
    return boost::apply_visitor(*this, ptr->base_solver);
  }
#endif // DIPOLES
};

struct ShortRangeEnergyKernel
    : public boost::static_visitor<
          boost::optional<Solver::ShortRangeEnergyKernel>> {

  using kernel_type = result_type::value_type;

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
    return boost::apply_visitor(*this, ptr->base_solver);
  }
#endif // DIPOLES
};

inline boost::optional<Solver::ShortRangeForceKernel>
Solver::pair_force_kernel() const {
#ifdef DIPOLES
  if (solver) {
    auto const visitor = Dipoles::ShortRangeForceKernel();
    return boost::apply_visitor(visitor, *solver);
  }
#endif // DIPOLES
  return {};
}

inline boost::optional<Solver::ShortRangeEnergyKernel>
Solver::pair_energy_kernel() const {
#ifdef DIPOLES
  if (solver) {
    auto const visitor = Dipoles::ShortRangeEnergyKernel();
    return boost::apply_visitor(visitor, *solver);
  }
#endif // DIPOLES
  return {};
}

} // namespace Dipoles

#endif
