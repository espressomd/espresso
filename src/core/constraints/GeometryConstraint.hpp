#ifndef __GEOMETRY_CONSTRAINT_HPP
#define __GEOMETRY_CONSTRAINT_HPP

#include "Constraint.hpp"
#include "shapes/ShapeList.hpp"

namespace Constraints {

  enum ReflectionType { REFLECTION_NONE, REFLECTION_NORMAL, REFLECTION_NORMAL_TANGENTIAL };

  struct GeometryConstraint : public Constraint {
  public:
    GeometryConstraint() {};
    GeometryConstraint(int shape, bool _penetrable, ReflectionType _reflection_type, bool _tuneable_slip = false) : 
      penetrable(_penetrable), tuneable_slip(_tuneable_slip), reflection_type(_reflection_type), m_shape_id(shape) {
    	set_shape();
    };
    virtual ConstraintType type() { return CONSTRAINT_GEOMETRY; }
    virtual std::string name() { return Constraint::name() + std::string("GeometryConstraint::") + m_shape->name(); }
    void reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos);
    bool penetrable;
    bool tuneable_slip;
    ReflectionType reflection_type;
    void set_shape();
    int calculate_dist(const double *ppos, double *dist, double *vec) const { return m_shape->calculate_dist(ppos, dist, vec); }
    void set_shape(int id) { m_shape_id = id; set_shape(); }
    const Shapes::Shape &constget_shape() const { return *m_shape; }
  protected:
    int m_shape_id;
    Shapes::pointer_type m_shape;
  };
};

#endif
