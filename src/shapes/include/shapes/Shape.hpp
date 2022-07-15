/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
  /**
   * @brief Calculate the minimum distance and the corresponding
   * distance vector between a given position and the shape.
   * @param[in]  pos  Position for which to calculate the distance.
   * @param[out] dist Minimum distance between @p pos and the shape. Value
   *                  is negative when @c pos is inside the shape, zero if
   *                  @c pos is on the shape surface or positive if @c pos
   *                  is outside the shape.
   * @param[out] vec  Distance vector.
   */
  virtual void calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const = 0;
  /**
   * @brief Check whether the given point is inside the shape or not.
   * For the edge case where the point is on the surface (zero distance),
   * it is considered to be inside the shape.
   * @param pos Position to check.
   */
  virtual bool is_inside(Utils::Vector3d const &pos) const {
    Utils::Vector3d vec;
    double dist;
    calculate_dist(pos, dist, vec);
    return dist <= 0.0;
  }
  /**
   * @brief Rasterize a shape on a regular grid.
   * @param grid_size     Number of grid points in every direction.
   * @param grid_spacing  %Lattice distance.
   * @param grid_offset   %Lattice offset.
   * @return Flattened 3D matrix with 1's inside the shape and 0's outside.
   */
  std::vector<int> rasterize(Utils::Vector3i const &grid_size,
                             double grid_spacing, double grid_offset) const;
  virtual ~Shape() = default;
};

} /* namespace Shapes */

#endif
