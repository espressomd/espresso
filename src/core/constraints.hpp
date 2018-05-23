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
  int image[3];
  double pos[3];
  memcpy(pos, p->r.p, 3 * sizeof(double));
  memcpy(image, p->l.i, 3 * sizeof(int));

  fold_position(pos, image);
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
  auto pos = folded_position(*p);

  for (auto const &c : Constraints::constraints) {
    c->add_energy(p, pos.data(), energy);
  }
}

#endif
#endif
