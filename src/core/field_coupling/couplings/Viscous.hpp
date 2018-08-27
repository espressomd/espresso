#ifndef CORE_CONSTRAINTS_DETAIL_VISCOUS_HPP
#define CORE_CONSTRAINTS_DETAIL_VISCOUS_HPP

#include "Vector.hpp"

namespace FieldCoupling {
namespace Coupling {
class Viscous {
  double m_gamma;

public:
  static constexpr bool is_linear = true;

  Viscous(double gamma) : m_gamma(gamma) {}
  double &gamma() { return m_gamma; }
  double const &gamma() const { return m_gamma; }

  template <typename Particle>
  Vector3d operator()(Particle const &p, Vector3d const &field) const {
    return m_gamma * (field - p.m.v);
  }
};
} // namespace Coupling
} // namespace FieldCoupling

#endif
