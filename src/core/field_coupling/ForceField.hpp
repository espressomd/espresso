#ifndef EXTERNAL_FIELD_FORCE_HPP
#define EXTERNAL_FIELD_FORCE_HPP

#include "detail/Base.hpp"
#include "detail/BindCoupling.hpp"

#include "Vector.hpp"

namespace FieldCoupling {
template <typename Coupling, typename Field>
class ForceField : public detail::Base<Coupling, Field> {
  using Base = detail::Base<Coupling, Field>;

  using Base::m_coupling;
  using Base::m_field;

public:
  using Base::Base;

  template <typename Particle>
  Vector3d force(const Particle &p, const Vector3d &folded_pos) const {
    using detail::make_bind_coupling;
    return m_field(make_bind_coupling(m_coupling, p), folded_pos);
  }
};
} /* namespace FieldCoupling */

#endif
