/*
  Copyright (C) 2016-2018 The ESPResSo project

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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "MeanVarianceCalculator.hpp"

namespace Accumulators {
void MeanVarianceCalculator::update() { m_acc(m_obs->operator()()); }

std::vector<double> MeanVarianceCalculator::get_mean() {
  return m_acc.get_mean();
}

std::vector<double> MeanVarianceCalculator::get_variance() {
  return m_acc.get_variance();
}

std::string MeanVarianceCalculator::get_internal_state() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << m_acc;

  return ss.str();
}

void MeanVarianceCalculator::set_internal_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  ia >> m_acc;
}
} // namespace Accumulators
