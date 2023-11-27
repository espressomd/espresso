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
#include "system/System.hpp"

#include <utils/math/sqr.hpp>
#include <utils/serialization/array.hpp>

#include <boost/mpi/collectives/reduce.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace Observables {
class LBFluidPressureTensor : public Observable {
public:
  std::vector<std::size_t> shape() const override { return {3, 3}; }
  std::vector<double>
  operator()(boost::mpi::communicator const &comm) const override {
    auto const &lb = System::get_system().lb;
    auto const pressure_conv = 1. / (lb.get_agrid() * Utils::sqr(lb.get_tau()));
    auto const local_tensor = lb.get_pressure_tensor() * pressure_conv;
    std::remove_const_t<decltype(local_tensor)> tensor;
    boost::mpi::reduce(comm, local_tensor, tensor, std::plus<>(), 0);
    return tensor.as_vector();
  }
};

} // Namespace Observables

#endif
