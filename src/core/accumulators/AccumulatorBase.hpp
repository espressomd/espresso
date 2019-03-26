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
  void auto_update();
  int &delta_N() { return m_delta_N; };
  virtual ~AccumulatorBase() = default;

private:
  virtual void update() = 0;
  // Number of timesteps between automatic updates.
  int m_delta_N;
  int m_counter = 0;
};

inline void AccumulatorBase::auto_update() {
  if (m_counter % m_delta_N == 0) {
    update();
  }
  ++m_counter;
}
} // namespace Accumulators

#endif
