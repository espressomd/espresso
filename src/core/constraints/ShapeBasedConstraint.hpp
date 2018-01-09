#ifndef CONSTRAINTS_SHAPEBASEDCONSTRAINT_HPP
#define CONSTRAINTS_SHAPEBASEDCONSTRAINT_HPP

#include <memory>

#include "Constraint.hpp"
#include "energy.hpp"
#include "particle_data.hpp"
#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"

namespace Constraints {

class ShapeBasedConstraint : public Constraint {
public:
  enum class ReflectionType { NONE, NORMAL, NORMAL_TANGENTIAL };

  ShapeBasedConstraint()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_reflection_type(ReflectionType::NONE), m_penetrable(false),
        m_only_positive(false), m_tuneable_slip(0), m_type(-1) {
    ShapeBasedConstraint::reset_force();
  }

  virtual void add_energy(Particle *p, double *folded_pos,
                  Observable_stat &energy) const override;

  virtual void add_force(Particle *p, double *folded_pos) override;

  /* Calculate distance from the constraint */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  Shapes::Shape const &shape() const { return *m_shape; }

  ReflectionType const &reflection_type() const;

  void reset_force() override { m_local_force = Vector3d{0, 0, 0}; }
  int &only_positive() { return m_only_positive; }
  int &penetrable() { return m_penetrable; }
  int &type() { return m_type; }
  
  void set_type(const int &type) {
    m_type = type;
    make_particle_type_exist(m_type);
  }

  Vector3d total_force() const;

private:
  /** Private methods */
  void reflect_particle(Particle *p, const double *distance_vector,
                        const double *folded_pos) const;

  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape;

  ReflectionType m_reflection_type;
  int m_penetrable;
  int m_only_positive;
  int m_tuneable_slip;
  int m_type;
  Vector3d m_local_force;
};

} /* namespace Constaints */

#endif
