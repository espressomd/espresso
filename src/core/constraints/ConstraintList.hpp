#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <map>

class ConstraintList : public std::map<int, Constraints::Constraint *> {
public:
  ConstraintList() : n_constraints(0), m_next_id(0) {}
  int add_constraint(Constraints::Constraint *c);
  void remove_constraint(int i);
  void add_forces(Particle *p);
  void add_energies(Particle *p);
  void init_forces();
  int n_constraints;
private:
  int m_next_id;
};
#endif
