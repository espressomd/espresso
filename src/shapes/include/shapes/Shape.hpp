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

#ifndef SHAPES_SHAPE_HPP
#define SHAPES_SHAPE_HPP

#include <utils/Vector.hpp>

namespace Shapes {

class Shape {
public:
  virtual void calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const = 0;
  virtual bool is_inside(Utils::Vector3d const &pos) const {
    Utils::Vector3d vec;
    double dist;
    calculate_dist(pos, dist, vec);
    return dist > 0.0;
  }
  virtual ~Shape() = default;
};

} /* namespace Shapes */

#endif
