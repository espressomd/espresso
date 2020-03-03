/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <vector>

namespace Observables {
class LBFluidStress : public Observable {
public:
  std::vector<size_t> shape() const override { return {3, 3}; }
  std::vector<double> operator()() const override {

    auto const unit_conversion =
        1. / (lb_lbfluid_get_agrid() * pow(lb_lbfluid_get_tau(), 2));
    auto const lower_triangle = lb_lbfluid_get_stress() * unit_conversion;
    return {lower_triangle[0], lower_triangle[1], lower_triangle[3],
            lower_triangle[1], lower_triangle[2], lower_triangle[4],
            lower_triangle[3], lower_triangle[4], lower_triangle[5]};
  }
};

} // Namespace Observables

#endif
