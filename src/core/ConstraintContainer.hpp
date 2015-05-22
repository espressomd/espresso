#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <vector>

class ConstraintContainer : public std::vector<ConstraintClass::Constraint *> {
public:
  int add_constraint(ConstraintClass::Constraint *c);
  bool remove_constraint(int i);
  void add_forces(Particle *p);
  void add_energies(Particle *p);
};
#endif
