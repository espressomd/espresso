/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef LBBOUNDARIES_LBBOUNDARY_HPP
#define LBBOUNDARIES_LBBOUNDARY_HPP

#include <memory>

#include "config.hpp"
#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"
#include "utils/Vector.hpp"

namespace LBBoundaries {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
int lbboundary_get_force(void *lbb, double *f);
void lb_init_boundaries();
#endif
class LBBoundary {
public:
  LBBoundary()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_velocity(Vector3d{0, 0, 0}), m_force(Vector3d{0, 0, 0}) {
#ifdef EK_BOUNDARIES
    m_charge_density = 0.0;
    m_net_charge = 0.0;
#endif
  }

  /* Calculate distance from the lbboundary */
  void calc_dist(const Vector3d &pos, double *dist, double *vec) const {
    m_shape->calculate_dist(pos, dist, vec);
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  void set_velocity(const Vector3d &velocity) {
    m_velocity = velocity;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    lb_init_boundaries();
#endif
  }
  void reset_force() { m_force = Vector3d{0, 0, 0}; }

  Shapes::Shape const &shape() const { return *m_shape; }
  Vector3d &velocity() { return m_velocity; }
  Vector3d &force() { return m_force; }
  Vector3d get_force() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    double tmp_force[3];
    lbboundary_get_force(this, tmp_force);
    return Vector3d{tmp_force[0], tmp_force[1], tmp_force[2]};
#else
    throw std::runtime_error("Needs LB_BOUNDARIES or LB_BOUNDARIES_GPU.");
#endif
  }

#ifdef EK_BOUNDARIES // TODO: ugly. Better would be a class EKBoundaries,
                     // deriving from LBBoundaries, but that requires completely
                     // different initialization infrastructure.
  void set_charge_density(float charge_density) {
    m_charge_density = charge_density;
  }
  void set_net_charge(float net_charge) { m_net_charge = net_charge; }

  float &charge_density() { return m_charge_density; }
  float &net_charge() { return m_net_charge; }
#endif

private:
  /** Private methods */
  /* The actual boundary */
  std::shared_ptr<::LBBoundaries::LBBoundary>
      m_lbboundary; // TODO probably a brainfart (and named wrong)

  /** Private data members */
  std::shared_ptr<Shapes::Shape>
      m_shape; // TODO: I don't like this being a pointer just to get around the
               // virtual limitations
  Vector3d m_velocity;
  Vector3d m_force;

#ifdef EK_BOUNDARIES // TODO: ugly. Better would be a class EKBoundaries,
                     // deriving from LBBoundaries, but that requires completely
                     // different initialization infrastructure.
  float m_charge_density;
  float m_net_charge;
#endif
};

} /* namespace LBBoundaries */

#endif
