#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "ObjectRegistry.hpp"
#include "config.hpp"

#include <memory>
#include <vector>

#include "constraints/Constraint.hpp"
#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
extern ObjectRegistry<std::vector<std::shared_ptr<Constraint>>> constraints;
}
#ifdef CONSTRAINTS

inline void add_constraints_forces(Particle *p) {

  // Copy position and image count so as not to modify particle when folding
  // position
  Vector3d pos=folded_position(p);
  for (auto const &c : Constraints::constraints) {
    c->add_force(p, pos);
  }
}

inline void init_constraint_forces() {
  for (auto const &c : Constraints::constraints) {
    c->reset_force();
  }
}

inline void add_constraints_energy(Particle *p) {

  //Vector3d pos=folded_position(p);
  for (auto const &c : Constraints::constraints) {
    c->add_energy(p, p->r.p, energy);
  }
}

#endif
#endif
