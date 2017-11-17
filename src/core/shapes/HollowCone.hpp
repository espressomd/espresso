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

#ifndef __HOLLOWCONE_HPP
#define __HOLLOWCONE_HPP

#include "Shape.hpp"

namespace Shapes {

class HollowCone : public Shape {
public:
  HollowCone()
      : m_position_x(0.0), m_position_y(0.0), m_position_z(0.0),
        m_orientation_x(0.0), m_orientation_y(0.0), m_orientation_z(0.0),
        m_outer_radius(0.0), m_inner_radius(0.0), m_width(0.0),
        m_opening_angle(0.0), m_direction(0.0) {}

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  double &position_x() { return m_position_x; }
  double &position_y() { return m_position_y; }
  double &position_z() { return m_position_z; }

  double &orientation_x() { return m_orientation_x; }
  double &orientation_y() { return m_orientation_y; }
  double &orientation_z() { return m_orientation_z; }

  double &outer_radius() { return m_outer_radius; }
  double &inner_radius() { return m_inner_radius; }
  double &width() { return m_width; }
  double &opening_angle() { return m_opening_angle; }
  double &direction() { return m_direction; }

private:
  /** Hollow cone position. */
  double m_position_x;
  double m_position_y;
  double m_position_z;

  /** Hollow cone orientation. */
  double m_orientation_x;
  double m_orientation_y;
  double m_orientation_z;

  /** Hollow cone dimensions. */
  double m_outer_radius;
  double m_inner_radius;
  double m_width;
  double m_opening_angle;

  /** Inside/Outside (+1 outside -1 inside interaction direction)*/
  double m_direction;
};
}

#endif
