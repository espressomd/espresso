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
#ifndef OBSERVABLES_LB_FLUID_STRESS_HPP
#define OBSERVABLES_LB_FLUID_STRESS_HPP

#include "Observable.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include <utils/math/sqr.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
class LBFluidPressureTensor : public Observable {
public:
  std::vector<std::size_t> shape() const override { return {3, 3}; }
  std::vector<double> operator()() const override {
    auto const unit_conversion =
        1. / (LB::get_agrid() * Utils::sqr(LB::get_tau()));
    auto const tensor = LB::get_pressure_tensor() * unit_conversion;
    return tensor.as_vector();
  }
};

} // Namespace Observables

#endif
