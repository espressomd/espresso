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
#ifndef OBSERVABLES_PRESSURE_HPP
#define OBSERVABLES_PRESSURE_HPP

#include "Observable.hpp"
#include "pressure.hpp"
#include <cstddef>
#include <vector>

namespace Observables {

class Pressure : public Observable {
public:
  std::vector<std::size_t> shape() const override { return {1}; }
  std::vector<double> operator()() const override {
    auto const ptensor = observable_compute_pressure_tensor();
    std::vector<double> res{1};
    res[0] = (ptensor[0] + ptensor[4] + ptensor[8]) / 3.;
    return res;
  }
};

} // Namespace Observables

#endif
