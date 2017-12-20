/*
  Copyright (C) 2016,2017 The ESPResSo project

  This file is part of ESPResSo.


  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _ACCUMULATORS_ACCUMULATOR_H
#define _ACCUMULATORS_ACCUMULATOR_H

/* Needed for transform_iterator to work with
   lambdas on older compilers. */
#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "observables/Observable.hpp"

#include <memory>

namespace Accumulators {
namespace ba = boost::accumulators;
typedef ba::accumulator_set<std::vector<double>,
                            ba::stats<ba::tag::mean, ba::tag::variance>>
    acc;

class Accumulator {
public:
  // The accumulator struct has to be initialized with the correct vector size,
  // therefore the order of init is important.
  Accumulator(std::shared_ptr<Observables::Observable> const& obs)
      : m_obs(obs), m_acc(std::vector<double>(obs->n_values())) {}

  int update();
  std::vector<double> get_mean();
  std::vector<double> get_variance();

private:
  std::shared_ptr<Observables::Observable> m_obs;
  acc m_acc;
};

} // namespace Accumulators

#endif
