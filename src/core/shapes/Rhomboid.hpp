/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef __RHOMBOID_HPP
#define __RHOMBOID_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class Rhomboid : public Shape {
public:
  Rhomboid()
      : m_pos({0.0, 0.0, 0.0}), m_a({0.0, 0.0, 0.0}), m_b({0.0, 0.0, 0.0}),
        m_c({0.0, 0.0, 0.0}), m_direction(0.0) {}

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

  Utils::Vector3d &pos() { return m_pos; }
  Utils::Vector3d &a() { return m_a; }
  Utils::Vector3d &b() { return m_b; }
  Utils::Vector3d &c() { return m_c; }
  double &direction() { return m_direction; }

private:
  /** corner of the rhomboid */
  Utils::Vector3d m_pos;
  /** edges adjacent to the corner */
  Utils::Vector3d m_a;
  Utils::Vector3d m_b;
  Utils::Vector3d m_c;
  /** rhomboid direction. (+1 outside -1 inside interaction direction)*/
  double m_direction;
};
} // namespace Shapes

#endif
