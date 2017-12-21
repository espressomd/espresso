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
#include "Accumulator.hpp"
#include "partCfg_global.hpp"

namespace Accumulators {
int Accumulator::update() {
  m_acc(m_obs->operator()(partCfg()));
  return 0;
}

std::vector<double> Accumulator::get_mean() {
  return ba::mean(m_acc);
}

std::vector<double> Accumulator::get_variance() {
  return ba::variance(m_acc);
}

}
