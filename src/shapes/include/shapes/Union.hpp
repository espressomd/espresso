/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifndef SHAPES_UNION
#define SHAPES_UNION

#include <boost/algorithm/cxx11/all_of.hpp>

#include <memory>

#include "Shape.hpp"

namespace Shapes {

class Union : public Shape {
public:
  void add(std::shared_ptr<Shapes::Shape> const &s) {
    m_shapes.emplace_back(s);
  }

  void remove(std::shared_ptr<Shapes::Shape> const &s) {
    m_shapes.erase(std::remove(m_shapes.begin(), m_shapes.end(), s),
                   m_shapes.end());
  }

  /**
   * @brief Calculates the minimum of all distances and the corresponding
   * distance vector for a given position and any contained shape.
   * @param[in] pos Position for which to calculate the distance.
   * @param[out] dist Minimum distance between pos and any contained shape.
   * @param[out] vec Distance vector.
   */
  void calculate_dist(Utils::Vector3d const &pos, double &dist,
                      Utils::Vector3d &vec) const override {
    auto dist_compare = [&pos](std::pair<double, Utils::Vector3d> const &res,
                               std::shared_ptr<Shapes::Shape> const &s) {
      double d;
      Utils::Vector3d vec;
      (*s).calculate_dist(pos, d, vec);
      if (d < 0.0)
        throw std::domain_error(
            "Distance to Union not well-defined for given position!");
      if (d < res.first) {
        return std::make_pair(d, vec);
      }
      return res;
    };
    std::tie(dist, vec) =
        std::accumulate(m_shapes.begin(), m_shapes.end(),
                        std::make_pair(std::numeric_limits<double>::infinity(),
                                       Utils::Vector3d{}),
                        dist_compare);
  }

  bool is_inside(Utils::Vector3d const &pos) const override {
    return boost::algorithm::all_of(
        m_shapes, [&pos](auto const &s) { return s->is_inside(pos); });
  }

private:
  std::vector<std::shared_ptr<Shapes::Shape>> m_shapes;
};

} // namespace Shapes

#endif
