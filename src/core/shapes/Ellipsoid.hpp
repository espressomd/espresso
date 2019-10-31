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

#ifndef SHAPES_ELLIPSOID_HPP
#define SHAPES_ELLIPSOID_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class Ellipsoid : public Shape {
public:
  Ellipsoid()
      : m_center({0.0, 0.0, 0.0}), m_semiaxes({1.0, 1.0, 1.0}),
        m_direction(1.0) {}

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

  void set_semiaxis_a(const double &value) { m_semiaxes[0] = value; }
  void set_semiaxis_b(const double &value) {
    m_semiaxes[1] = value;
    m_semiaxes[2] = value;
  }

  // not exposed via interface; used in core unit test
  void set_semiaxis_c(const double &value) { m_semiaxes[2] = value; }

  Utils::Vector3d &center() { return m_center; }
  double &semiaxis_a() { return m_semiaxes[0]; }
  double &semiaxis_b() { return m_semiaxes[1]; }
  double &semiaxis_c() { return m_semiaxes[2]; }
  double &direction() { return m_direction; }

private:
  bool inside_ellipsoid(const Utils::Vector3d &ppos) const;
  double newton_term(const Utils::Vector3d &ppos, const double &l) const;

  Utils::Vector3d m_center;
  Utils::Vector3d m_semiaxes;
  double m_direction;
};
} // namespace Shapes

#endif
