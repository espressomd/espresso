#include "ConstraintContainer.hpp"

#include "GeometryConstraint.hpp"

void ConstraintContainer::init_forces() {
  for(iterator it = begin(); it != end(); ++it) {
    for(int i = 0; i < 3; i++)
      (*it)->total_force[i] = 0.0;
  }
}

int ConstraintContainer::add_constraint(Constraints::Constraint *c) {
  /* Check if c is a GeometryConstraint in which case we have
     to do a deep copy of the Shape object. */
  if(dynamic_cast<Constraints::GeometryConstraint *>(c)) {
    ;
  } else {
    ;
  }
}
