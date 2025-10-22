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

#include "Observable.hpp"
#include "Observable_stat.hpp"
#include "system/System.hpp"

#include <cstddef>
#include <vector>

namespace Observables {

class PressureTensor : public Observable {
public:
  std::vector<std::size_t> shape() const override { return {3u, 3u}; }
  std::vector<double>
  operator()(boost::mpi::communicator const &) const override {
    auto const obs = System::get_system().calculate_pressure();

    std::vector<double> result;
    result.reserve(9);
    for (std::size_t i = 0u; i < 9u; ++i) {
      result.emplace_back(obs->accumulate(0., i));
    }
    return result;
  }
};

} // Namespace Observables
