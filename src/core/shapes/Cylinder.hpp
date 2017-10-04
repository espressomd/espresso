/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef __CYLINDER_HPP
#define __CYLINDER_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {
class Cylinder : public Shape {
public:
  Cylinder()
      : m_pos({0.0, 0.0, 0.0}), m_axis({0.0, 0.0, 0.0}), m_length(0.0),
        m_direction(1.0), m_rad(0), m_open(false) {}
  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  Vector3d &pos() { return m_pos; }
  Vector3d &axis() { return m_axis; }

  double &rad() { return m_rad; }
  double &length() { return m_length; }
  double &direction() { return m_direction; }

protected:
  /** center of the cylinder. */
  Vector3d m_pos;
  /** Axis of the cylinder. */
  Vector3d m_axis;
  /** cylinder radius. */
  double m_rad;
  /** cylinder length. */
  double m_length;
  /** cylinder direction. (+1 outside -1 inside interaction direction) */
  double m_direction;
  /** ignore bottom and top cap of cylinder */
  bool m_open;
};
};

#endif
