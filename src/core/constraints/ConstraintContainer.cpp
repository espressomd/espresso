#include "ConstraintContainer.hpp"

#include "GeometryConstraint.hpp"

void ConstraintContainer::init_forces() {
  for(iterator it = begin(); it != end(); ++it) {
    for(int i = 0; i < 3; i++)
      (*it)->total_force[i] = 0.0;
  }
}

int ConstraintContainer::add_constraint(Constraints::Constraint *c) {
  /* Add c to local list */
  push_back(c);

  /* c is now the last element */
  c->id = size() - 1;

  /* Check if c is a GeometryConstraint in which case we have
     to do a deep copy of the Shape object. */
  if(dynamic_cast<Constraints::GeometryConstraint *>(c)) {
    ;
  } else {
    ;
  }

  return c->id;
}

void ConstraintContainer::remove_constraint(int i) {
  if(i < size()) {
    Constraints::Constraint *c = operator[](i);
      if(c)
        delete c;

    }
}
