/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

#include <iostream>
#include <limits>
#include <algorithm>
#include <numeric>

#define BOOST_TEST_MODULE RunningAverage test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../utils/statistics/RunningAverage.hpp"

/** random sequence */
#include "random_sequence.hpp"

using namespace Testing;

BOOST_AUTO_TEST_CASE(simple_tests) {
  Utils::Statistics::RunningAverage<double> avg;

  BOOST_CHECK(avg.n() == 0);
  BOOST_CHECK(avg.avg() == 0.0);
  BOOST_CHECK(avg.var() == 0.0);

  avg.add_sample(5.0);

  BOOST_CHECK(avg.n() == 1);
  BOOST_CHECK(avg.avg() == 5.0);
  BOOST_CHECK(avg.var() == 0.0);

  avg.add_sample(5.0);

  BOOST_CHECK(avg.n() == 2);
  BOOST_CHECK(avg.avg() == 5.0);
  BOOST_CHECK(avg.var() == 0.0);

  avg.clear();

  BOOST_CHECK(avg.n() == 0);
  BOOST_CHECK(avg.avg() == 0.0);
  BOOST_CHECK(avg.var() == 0.0);
}

BOOST_AUTO_TEST_CASE(simple_variance_check) {
  Utils::Statistics::RunningAverage<double> avg;

  avg.add_sample(0.0);
  avg.add_sample(5.0);

  /** Var should be <x**2> - <x>**2 = (0**2 + 5.0**2) / 2. - 2.5**2 = 12.5 - 6.25 = 6.25 */
  BOOST_CHECK(std::abs(avg.avg() - 2.5) <= std::numeric_limits<double>::epsilon());
  BOOST_CHECK(std::abs(avg.var() - 6.25) <= std::numeric_limits<double>::epsilon());

  /** Standard deviation should be sqrt(var()) */
  BOOST_CHECK(std::abs(avg.sig() - std::sqrt(avg.var())) <= std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(mean_and_variance) {
  Utils::Statistics::RunningAverage<double> running_average;
  const size_t sample_size = sizeof(RandomSequence::values) / sizeof(double);

  for(auto const& val: RandomSequence::values) {
    running_average.add_sample(val);
  }

  BOOST_CHECK(running_average.n() == sample_size);  
  
  /** Directly calculate the meann from the data */
  const double m_mean = std::accumulate(std::begin(RandomSequence::values),
                                        std::end(RandomSequence::values),
                                        0.0) / sample_size;

  BOOST_CHECK(std::abs(running_average.avg() - m_mean) <= 1e-12);

  /** Directly calculate the variance from the data */
  double m_var = 0.0;
  for(auto const& val: RandomSequence::values) {
    m_var += (val - m_mean)*(val - m_mean);
  }
  m_var /= sample_size;

  BOOST_CHECK(std::abs(running_average.var() - m_var) <= 1e-12);  
}
