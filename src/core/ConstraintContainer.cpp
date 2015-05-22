#include "ConstraintContainer.hpp"

#include "GeometryConstraint.hpp"

int ConstraintContainer::add_constraint(ConstraintClass::Constraint *c) {
  /* Check if c is a GeometryConstraint in which case we have
     to do a deep copy of the Shape object. */
  if(dynamic_cast<ConstraintClass::GeometryConstraint *>(c)) {
    ;
  } else {
    ;
  }
}
