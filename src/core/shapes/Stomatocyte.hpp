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

#ifndef __STOMATOCYTE_HPP
#define __STOMATOCYTE_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {
class Stomatocyte : public Shape {
public:
  Stomatocyte() : m_position({0., 0., 0.}), m_orientation({1., 0., 0.}),
                  m_outer_radius(0.0), m_inner_radius(0.0), m_layer_width(0.0),
                  m_direction(0.0) {}

  int calculate_dist(const double *ppos, double *dist, double *vec) const override;

  Vector3d const &position() const { return m_position; }
  void set_position(Vector3d const &position) {
      m_position = position;
      m_position_x = position[0];
      m_position_y = position[1];
      m_position_z = position[2];
  }

  Vector3d const &orientation() const { return m_orientation; }
  void set_orientation(Vector3d const &orientation) {
      m_orientation = orientation;
      m_orientation_x = orientation[0];
      m_orientation_y = orientation[1];
      m_orientation_z = orientation[2];
  }

  double  &outer_radius()  { return m_outer_radius; }
  double  &inner_radius()  { return m_inner_radius; }
  double  &layer_width()  { return m_layer_width; }
  
  double  &direction()  { return m_direction; }

private:
  Vector3d m_position;
  Vector3d m_orientation;

  /** Stomatocyte position. */
  double m_position_x;
  double m_position_y;
  double m_position_z;

  /** Stomatocyte orientation. */
  double m_orientation_x;
  double m_orientation_y;
  double m_orientation_z;

  /** Stomatocyte dimensions. */
  double m_outer_radius;
  double m_inner_radius;
  double m_layer_width;

  /** Inside/Outside (+1 outside -1 inside interaction direction)*/
  double m_direction;
};
}

#endif
