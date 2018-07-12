#ifndef CONSTRAINTS_EXTERNAL_POTENTIAL_HPP
#define CONSTRAINTS_EXTERNAL_POTENTIAL_HPP

#include "detail/Base.hpp"
#include "detail/BindCoupling.hpp"

#include "Vector.hpp"

namespace FieldCoupling {
template <typename Coupling, typename Field>
class PotentialField : public detail::Base<Coupling, Field> {
  using field_gradient_type = typename Field::gradient_type;
  using field_value_type = typename Field::value_type;

public:
  using Base::Base;

  Vector3d force(const Particle &p, const Vector3d &folded_pos) const {
    return m_field.gradient(
        [p, this](const field_gradient_type &field) {
          return m_coupling(p, field);
        },
        folded_pos);
  }

  double energy(const Particle &p, const Vector3d &folded_pos) const {
    return m_field(
        [p, this](const field_value_type &field) {
          return m_coupling(p, field);
        },
        folded_pos);
  }
};
} /* namespace FieldCoupling */

#endif
