#ifndef LBBOUNDARIES_LBBOUNDARY_HPP
#define LBBOUNDARIES_LBBOUNDARY_HPP

#include <memory>

#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"

namespace LBBoundaries {
class LBBoundary {
public:
  Constraint()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_velocity(Vector3d{0, 0, 0}) {}

  /* Calculate distance from the lbboundary */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  void set_velocity(Vector3d velocity) { m_velocity = velocity; }

  Shapes::Shape const &shape() const { return *m_shape; }
  Vector3d const &velocity() { return *m_velocity; }

private:
  /** Private methods */

  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape;
  Vector3d m_velocity;
};

} /* namespace LBBoundaries */

#endif
