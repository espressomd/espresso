#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "config.hpp"

#ifdef CONSTRAINTS

#include <memory>
#include <vector>

#include "energy.hpp"
#include "particle_data.hpp"

#include "constraints/Constraints.hpp"

namespace Constraints {
extern Constraints<std::vector<std::shared_ptr<Constraint>>> constraints;
}

inline void add_constraints_forces(Particle *p) {
  double folded_position[3];
  int image[3];

  fold_position(p->r.p, folded_position, image);
  for (auto const &c : Constraints::constraints) {
    c->add_force(p, folded_position);
  }
}

inline void init_constraint_forces() {
  Constraints::constraints.reset_forces();
}

inline void add_constraints_energy(Particle *p) {
  double folded_position[3];
  int image[3];

  fold_position(p->r.p, folded_position, image);
  for (auto const &c : Constraints::constraints) {
    c->add_energy(p, folded_position, energy);
  }
}

#endif
#endif
