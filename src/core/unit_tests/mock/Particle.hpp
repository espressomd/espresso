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
#ifndef CORE_UNIT_TESTS_MOCK_PARTICLE_HPP
#define CORE_UNIT_TESTS_MOCK_PARTICLE_HPP

#include <array>

namespace Testing {
class Particle {

public:
  Particle() : m_id(-1) {}
  explicit Particle(int id) : m_id(id) {}

  int identity() const { return m_id; }

  int m_id;
};

class Position {
  std::array<double, 3> m_pos;

public:
  std::array<double, 3> &position() { return m_pos; }
};
} // namespace Testing

#endif
