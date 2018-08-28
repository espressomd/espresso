#ifndef CONSTRAINTS_HOMOGENEOUSMAGNETICFIELD_HPP
#define CONSTRAINTS_HOMOGENEOUSMAGNETICFIELD_HPP

#include "Constraint.hpp"
#include "particle_data.hpp"

namespace Constraints {

class HomogeneousMagneticField : public Constraint {
public:
  HomogeneousMagneticField() : m_field({1., 0., 0.}) {}

  void set_H(Vector3d const &H) { m_field = H; }

  Vector3d const &H() const { return m_field; }

  void add_energy(const Particle &p, const Vector3d &folded_pos,
                  Observable_stat &energy) const override;

  ParticleForce force(const Particle &p, const Vector3d &folded_pos) override;

  bool fits_in_box(Vector3d const &box) const override { return true; }

private:
  Vector3d m_field;
};

} /* namespace Constraints */

#endif
