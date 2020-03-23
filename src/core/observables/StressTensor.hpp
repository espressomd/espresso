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
#ifndef OBSERVABLES_STRESSTENSOR_HPP

#include "Observable.hpp"
#include "Particle.hpp"
#include "pressure.hpp"
#include <vector>

namespace Observables {

class StressTensor : public Observable {
public:
  std::vector<size_t> shape() const override { return {3, 3}; }
  std::vector<double> operator()() const override {
    std::vector<double> res(n_values());
    observable_compute_stress_tensor(1, res.data());
    return res;
  }
};

} // Namespace Observables

#endif
