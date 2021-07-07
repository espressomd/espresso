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
#ifndef LBBOUNDARIES_LBBOUNDARY_HPP
#define LBBOUNDARIES_LBBOUNDARY_HPP

#include "config.hpp"

#include <shapes/NoWhere.hpp>
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <memory>

namespace LBBoundaries {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
class LBBoundary;
void lb_init_boundaries();
#endif
class LBBoundary {
public:
  LBBoundary()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_velocity(Utils::Vector3d{0, 0, 0}),
        m_force(Utils::Vector3d{0, 0, 0}) {}

  /* Calculate distance from the lbboundary */
  void calc_dist(const Utils::Vector3d &pos, double &dist,
                 Utils::Vector3d &vec) const {
    m_shape->calculate_dist(pos, dist, vec);
  }

  double calc_dist(const Utils::Vector3d &pos) const {
    double dist;
    Utils::Vector3d tmp;
    calc_dist(pos, dist, tmp);
    return dist;
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  void set_velocity(const Utils::Vector3d &velocity) {
    m_velocity = velocity;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    lb_init_boundaries();
#endif
  }
  void reset_force() { m_force = Utils::Vector3d{0, 0, 0}; }

  Shapes::Shape const &shape() const { return *m_shape; }
  Utils::Vector3d &velocity() { return m_velocity; }
  Utils::Vector3d &force() { return m_force; }
  Utils::Vector3d get_force() const;

private:
  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape;
  Utils::Vector3d m_velocity;
  Utils::Vector3d m_force;
};

} /* namespace LBBoundaries */

#endif
