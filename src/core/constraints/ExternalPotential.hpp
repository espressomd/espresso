#ifndef CONSTRAINTS_EXTERNAL_POTENTIAL_HPP
#define CONSTRAINTS_EXTERNAL_POTENTIAL_HPP

#include "Constraint.hpp"
#include "particle_data.hpp"

#include "detail/BindCoupling.hpp"

namespace Constraints {
template <typename Coupling, typename Field>
class ExternalPotential : public Constraint {
  using field_gradient_type = typename Field::gradient_type;
  using field_value_type = typename Field::value_type;

  Coupling m_coupling;
  Field m_field;

public:
  template <typename CouplingRef, typename FieldRef>
  ExternalPotential(CouplingRef &&coupling, FieldRef &&field)
      : m_coupling(coupling), m_field(field) {}

  Coupling &coupling() { return m_coupling; }
  Field &field() { return m_field; }

  void add_force(Particle *p, double *folded_pos) override {
    const Vector3d force = m_field.gradient(
        [&p, this](const field_gradient_type &field) {
          return m_coupling(*p, field);
        },
        Vector3d{folded_pos, folded_pos + 3});

    for (int i = 0; i < 3; i++) {
      p->f.f[i] += force[i];
    }
  }

  double energy(Particle *p, double *folded_pos) const {
    return m_field([&p, this](const field_value_type &field) {
      return m_coupling(*p, field);
    });
  }
};
} /* namespace Constraints */

#endif
