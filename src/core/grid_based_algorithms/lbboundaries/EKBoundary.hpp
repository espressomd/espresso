#ifndef ESPRESSO_EKBOUNDARY_HPP
#define ESPRESSO_EKBOUNDARY_HPP

#include <shapes/NoWhere.hpp>
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <memory>

namespace EKBoundaries {

class EKBoundary {
public:
  EKBoundary() : m_shape(std::make_shared<Shapes::NoWhere>()) {}

  /* Calculate distance from the boundary */
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

  Shapes::Shape const &shape() const { return *m_shape; }

private:
  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape;
};

} // namespace EKBoundaries

#endif // ESPRESSO_EKBOUNDARY_HPP
