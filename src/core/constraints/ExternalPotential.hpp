#ifndef CONSTRAINTS_EXTERNAL_POTENTIAL_HPP
#define CONSTRAINTS_EXTERNAL_POTENTIAL_HPP

#include "Constraint.hpp"
#include "field_coupling/PotentialField.hpp"

namespace Constraints {
/**
 * @brief Constraint interface for ExternalField::PotentialField.
 */
template <typename Coupling, typename Field>
class ExternalPotential : public Constraint {
  FieldCoupling::PotentialField<Coupling, Field> impl;

public:
  template <typename... Args>
  ExternalPotential(Args... args)
      : impl(std::forward<Args>(args)...) {}

  Coupling &coupling() { return impl.coupling(); }
  Field &field() { return impl.field(); }

  ParticleForce force(const Particle &p, Vector3d const &folded_pos) override {
    return impl.force(p, folded_pos);
  }
};
} /* namespace Constraints */
#endif
