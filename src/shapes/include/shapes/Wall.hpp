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

#ifndef SHAPES_WALL_HPP
#define SHAPES_WALL_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class Wall : public Shape {
public:
  Wall() : m_n({1., 0., 0.}), m_d(0.0) {}

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

  void set_normal(const Utils::Vector3d &normal) {
    m_n = normal;

    m_n.normalize();
  }

  Utils::Vector3d const &n() const { return m_n; }

  double const &d() const { return m_d; }

  double &d() { return m_d; }

private:
  /** normal vector on the plane */
  Utils::Vector3d m_n;
  /** distance of the wall from the origin. */
  double m_d;
};

} /* namespace Shapes */

#endif
