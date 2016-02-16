#ifndef __CONSTRAINT_CONTAINER_HPP
#define __CONSTRAINT_CONTAINER_HPP

#include "Constraint.hpp"

#ifdef HAVE_CXX11
#include <unordered_set>
#else
#include <set>
#endif

#include <functional>
#include <memory>

namespace Constraints {

#ifdef HAVE_CXX11
typedef std::unordered_set<std::shared_ptr<Constraint>> set_type;
#else
typedef std::set<std::shared_ptr<Constraint>> set_type;
#endif

class ConstraintList : public set_type {
 public:
  enum Action { ADD, DELETE };
  
  void add_forces(Particle *p);
  void add_energies(Particle *p);
  void init_forces();
  double min_dist(double pos[3]);
 private:
  set_type m_set;
};

  extern ConstraintList list;
}
  
#endif
