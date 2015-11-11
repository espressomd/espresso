#ifndef __GEOMETRY_CONSTRAINT_HPP
#define __GEOMETRY_CONSTRAINT_HPP

#include "Constraint.hpp"
#include "shapes/Shape.hpp"

namespace Constraints {

enum ReflectionType { REFLECTION_NONE, REFLECTION_NORMAL, REFLECTION_NORMAL_TANGENTIAL };

struct GeometryConstraint : public Constraint {
 public:
  GeometryConstraint(Shapes::Shape &shape, bool _penetrable = false, ReflectionType _reflection_type = REFLECTION_NONE, bool _tuneable_slip = false) : 
      penetrable(_penetrable), tuneable_slip(_tuneable_slip), reflection_type(_reflection_type), m_shape(shape) {}
  virtual ConstraintType type() { return CONSTRAINT_GEOMETRY; }
  virtual std::string name() { return Constraint::name() + std::string("GeometryConstraint::") + m_shape.name(); }
  void reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos);
  bool penetrable;
  bool tuneable_slip;
  ReflectionType reflection_type;
  inline int calculate_dist(const double *ppos, double *dist, double *vec) const { return m_shape.calculate_dist(ppos, dist, vec); }
  void set_shape(Shapes::Shape &shape) { m_shape = shape; }
  const Shapes::Shape &constget_shape() const { return m_shape; }
 protected:
  Shapes::Shape &m_shape;
};
};

#endif
