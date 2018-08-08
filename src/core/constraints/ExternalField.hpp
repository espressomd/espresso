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

  const Coupling &coupling() const { return impl.coupling(); }
  const Field &field() const { return impl.field(); }

  void add_energy(const Particle &, const Vector3d &,
                  Observable_stat &) const override {}

  ParticleForce force(const Particle &p, Vector3d const &folded_pos) override {
    return impl.force(p, folded_pos);
  }

  bool fits_in_box(Vector3d const &box) const override {
    return impl.field().fits_in_box(box);
  }
};
} /* namespace Constraints */

#endif
