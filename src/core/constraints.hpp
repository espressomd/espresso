#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "config.hpp"
#include "ObjectRegistry.hpp"


#include <memory>
#include <vector>

#include "energy.hpp"
#include "particle_data.hpp"
#include "constraints/Constraint.hpp"


namespace Constraints {
extern ObjectRegistry<std::vector<std::shared_ptr<Constraint>>> constraints;
}
#ifdef CONSTRAINTS

inline void add_constraints_forces(Particle *p) {
  int image[3];

  fold_position(p->r.p, image);
  for (auto const &c : Constraints::constraints) {
    c->add_force(p, p->r.p);
  }
}

inline void init_constraint_forces() {
  for (auto const &c : Constraints::constraints) {
    c->reset_force();
  }
}

inline void add_constraints_energy(Particle *p) {
  int image[3];

  fold_position(p->r.p, image);
  for (auto const &c : Constraints::constraints) {
    c->add_energy(p, p->r.p, energy);
  }
}

#endif
#endif
