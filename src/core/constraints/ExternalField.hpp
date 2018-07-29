#ifndef CONSTRAINTS_EXTERNAL_FIELD_HPP
#define CONSTRAINTS_EXTERNAL_FIELD_HPP

#include "Constraint.hpp"
#include "field_coupling/ForceField.hpp"

namespace Constraints {
/**
 * @brief Constraint interface for ExternalField::ForceField.
 */
template <typename Coupling, typename Field>
class ExternalField : public Constraint {
  FieldCoupling::ForceField<Coupling, Field> impl;

public:
  template <typename... Args>
  ExternalField(Args... args) : impl(std::forward<Args>(args)...) {}

  Coupling &coupling() { return impl.coupling(); }
  Field &field() { return impl.field(); }

  ParticleForce force(const Particle &p, Vector3d const &folded_pos) override {
    return impl.force(p, folded_pos);
  }
};
} /* namespace Constraints */

#endif
