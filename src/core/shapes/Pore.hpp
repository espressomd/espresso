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

#ifndef __PORE_HPP
#define __PORE_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {
class Pore : public Shape {
public:
  Pore()
      : m_pos({0.0, 0.0, 0.0}), m_axis({0.0, 0.0, 0.0}), m_rad_left(0.0),
        m_rad_right(0.0), m_smoothing_radius(0.0), m_length(0.0),
        m_outer_rad_left(0.0), m_outer_rad_right(0.0) {}

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  Vector3d &pos() { return m_pos; }
  Vector3d &axis() { return m_axis; }
  double &rad_left() { return m_rad_left; }
  double &rad_right() { return m_rad_right; }
  double &smoothing_radius() { return m_smoothing_radius; }
  double &length() { return m_length; }
  double &outer_rad_left() { return m_outer_rad_left; }
  double &outer_rad_right() { return m_outer_rad_right; }

private:
  /** center of the cylinder. */
  Vector3d m_pos;
  /** Axis of the cylinder .*/
  Vector3d m_axis;
  /** cylinder radius. */
  double m_rad_left;
  double m_rad_right;
  double m_smoothing_radius;
  /** cylinder length. **/
  double m_length;
  double m_outer_rad_left;
  double m_outer_rad_right;
};
}

#endif
