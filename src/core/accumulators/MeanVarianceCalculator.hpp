/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#include "AccumulatorBase.hpp"
#include "observables/Observable.hpp"

#include <utils/Accumulator.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Accumulators {

class MeanVarianceCalculator : public AccumulatorBase {
public:
  // The accumulator struct has to be initialized with the correct vector size,
  // therefore the order of init is important.
  MeanVarianceCalculator(::System::System const *system, int delta_N,
                         std::shared_ptr<Observables::Observable> obs)
      : AccumulatorBase(system, delta_N), m_obs(obs), m_acc(obs->n_values()) {}

  void update(boost::mpi::communicator const &comm) override;
  std::vector<double> mean();
  std::vector<double> variance();
  std::vector<double> std_error();
  std::string get_internal_state() const final;
  void set_internal_state(std::string const &) final;
  std::vector<std::size_t> shape() const override { return m_obs->shape(); }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  Utils::Accumulator m_acc;
};

} // namespace Accumulators
