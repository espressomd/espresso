/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#ifndef __SPHERE_HPP
#define __SPHERE_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class Sphere : public Shape {
public:
  Sphere() : m_pos({0.0, 0.0, 0.0}), m_rad(0.0), m_direction(1.0) {}

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

  Utils::Vector3d &pos() { return m_pos; }
  double &rad() { return m_rad; }
  double &direction() { return m_direction; }

private:
  Utils::Vector3d m_pos;
  double m_rad;
  double m_direction;
};
} // namespace Shapes

#endif
