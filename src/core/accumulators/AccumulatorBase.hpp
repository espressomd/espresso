/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef CORE_ACCUMULATORS_ACCUMULATORBASE
#define CORE_ACCUMULATORS_ACCUMULATORBASE

namespace Accumulators {

class AccumulatorBase {
public:
  explicit AccumulatorBase(int delta_N = 1) : m_delta_N(delta_N){};
  int &delta_N() { return m_delta_N; };
  virtual ~AccumulatorBase() = default;

  virtual void update() = 0;

private:
  // Number of timesteps between automatic updates.
  int m_delta_N;
};
} // namespace Accumulators

#endif
