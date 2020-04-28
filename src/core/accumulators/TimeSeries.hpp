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
#ifndef CORE_ACCUMULATORS_TIMESERIES_HPP
#define CORE_ACCUMULATORS_TIMESERIES_HPP

#include "AccumulatorBase.hpp"
#include "observables/Observable.hpp"

#include <memory>

namespace Accumulators {

/**
 * @brief Record values of an observable.
 *
 * This is a very simple accumulator that stores
 * the current value of an observable every time
 * it is updated.
 *
 */
class TimeSeries : public AccumulatorBase {
public:
  TimeSeries(std::shared_ptr<Observables::Observable> obs, int delta_N)
      : AccumulatorBase(delta_N), m_obs(std::move(obs)) {}

  void update() override;
  std::string get_internal_state() const;
  void set_internal_state(std::string const &);

  const std::vector<std::vector<double>> &time_series() const { return m_data; }
  std::vector<std::size_t> shape() const {
    std::vector<std::size_t> shape{m_data.size()};
    auto obs_shape = m_obs->shape();
    shape.insert(shape.end(), obs_shape.begin(), obs_shape.end());
    return shape;
  }
  void clear() { m_data.clear(); }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  std::vector<std::vector<double>> m_data;
};

} // namespace Accumulators

#endif
