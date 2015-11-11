#ifndef __INTERACTION_CONSTRAINT_HPP
#define __INTERACTION_CONSTRAINT_HPP

#include "GeometryConstraint.hpp"

#include "interaction_data.hpp"

namespace Constraints {
  struct InteractionConstraint : public GeometryConstraint {
    InteractionConstraint(Shapes::Shape &shape, bool _penetrable, ReflectionType _reflection_type, bool _only_positive, int ia_type) : GeometryConstraint(shape, _penetrable, _reflection_type), only_positive(_only_positive) {
      part_rep.p.type = ia_type;
    }
    virtual ConstraintType type() { return CONSTRAINT_INTERACTION; }
    void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    Particle part_rep;
    bool only_positive;
  };
}

#endif
