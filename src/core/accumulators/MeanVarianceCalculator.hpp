/*
 * Copyright (C) 2016-2019 The ESPResSo project
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
#ifndef _ACCUMULATORS_ACCUMULATOR_H
#define _ACCUMULATORS_ACCUMULATOR_H

#include "AccumulatorBase.hpp"
#include "observables/Observable.hpp"
#include <utils/Accumulator.hpp>

namespace Accumulators {

class MeanVarianceCalculator : public AccumulatorBase {
public:
  // The accumulator struct has to be initialized with the correct vector size,
  // therefore the order of init is important.
  MeanVarianceCalculator(std::shared_ptr<Observables::Observable> const &obs,
                         int delta_N)
      : AccumulatorBase(delta_N), m_obs(obs), m_acc(obs->n_values()) {}

  void update() override;
  std::vector<double> get_mean();
  std::vector<double> get_variance();
  /* Partial serialization of state that is not accessible
     via the interface. */
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);
  std::vector<std::size_t> shape() const { return m_obs->shape(); }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  ::Utils::Accumulator m_acc;
};

} // namespace Accumulators

#endif
