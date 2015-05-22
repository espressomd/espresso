#ifndef __GEOMETRY_CONSTRAINT_HPP
#define __GEOMETRY_CONSTRAINT_HPP

#include "Constraint.hpp"
#include "Shape.hpp"

namespace ConstraintClass {

enum ReflectionType { REFLECTION_NONE, REFLECTION_NORMAL, REFLECTION_NORMAL_TANGENTIAL };

struct GeometryConstraint : public Constraint {
public:
GeometryConstraint(ConstraintType _type, Shape &shape, bool _penetrable, ReflectionType _reflection_type) : 
  Constraint(_type), penetrable(_penetrable), reflection_type(_reflection_type), m_shape(shape) {};
bool penetrable;
ReflectionType reflection_type;
protected:
Shape &m_shape;
void reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos);
};

};

#endif
