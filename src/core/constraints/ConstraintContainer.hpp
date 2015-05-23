#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <vector>

class ConstraintContainer : public std::vector<Constraints::Constraint *> {
public:
  ConstraintContainer() : n_constraints(0) {}
  int add_constraint(Constraints::Constraint *c);
  void remove_constraint(int i);
  void add_forces(Particle *p);
  void add_energies(Particle *p);
  void init_forces();
  int n_constraints;
};
#endif
