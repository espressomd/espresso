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

#ifndef __ELLIPSOID_HPP
#define __ELLIPSOID_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {
class Ellipsoid : public Shape {
public:
  Ellipsoid() : m_pos({0.0, 0.0, 0.0}), m_semiaxis_a(1.0), m_semiaxis_b(1.0), m_semiaxis_c(1.0), m_direction(1.0) {}

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  Vector3d &pos() { return m_pos; }
  double &semiaxis_a() { return m_semiaxis_a; }
  double &semiaxis_b() { return m_semiaxis_b; }
  double &semiaxis_c() { return m_semiaxis_c; }
  double &direction() { return m_direction; }

private:
  Vector3d ClosestEllipsoidPoint(Vector3d ppos) const;
  double GetRoot(const double r0, const double r1, const double z0, const double z1, const double z2, double g) const;
  double GetRoot(const double r0, const double z0, const double z1, double g) const;
  double DistancePointEllipse(const double e0, const double e1, const double y0, const double y1, double& x0, double& d1) const;
  double RobustLength(const double v0, const double v1, const double v2) const;
  double RobustLength(const double v0, const double v1) const;

  Vector3d m_pos;
  double m_semiaxis_a;
  double m_semiaxis_b;
  double m_semiaxis_c;
  double m_direction;
};
}

#endif
