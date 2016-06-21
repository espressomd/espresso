#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
struct Constraint {
public:
  virtual void add_energy(Particle const *p, double const *folded_pos,
                          Observable_stat &energy) const = 0;
  virtual void add_force(Particle *p, double const *folded_pos) const = 0;
  /* Accumulated force excerted by the constraint */
  mutable std::array<double, 3> total_force;
  /* Human readable name */
  virtual ~Constraint() {}
};
}

#endif
