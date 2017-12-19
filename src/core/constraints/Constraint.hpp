#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:
  virtual void add_energy(Particle *p, double *folded_pos,
                  Observable_stat &energy) const {};

  virtual void add_force(Particle *p, double *folded_pos) = 0;
  virtual void reset_force() {};
};
} /* namespace Constaints */

#endif
