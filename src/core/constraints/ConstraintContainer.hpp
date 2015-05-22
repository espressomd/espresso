#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <vector>

class ConstraintContainer : public std::vector<Constraints::Constraint *> {
public:
  int add_constraint(Constraints::Constraint *c);
  bool remove_constraint(int i);
  void add_forces(Particle *p);
  void add_energies(Particle *p);
  void init_forces();
};
#endif
