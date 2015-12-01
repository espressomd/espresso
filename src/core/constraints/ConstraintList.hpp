#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#include <set>

namespace Constraints {

  class ConstraintList : public std::set<Constraints::Constraint *> {
  public:
    ConstraintList() {}
    void add_forces(Particle *p);
    void add_energies(Particle *p);
    void init_forces();
    double min_dist(double pos[3]);
  };

  extern ConstraintList list;
}
  
#endif
