#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:
  virtual void add_energy(const Particle &p, const Vector3d &folded_pos,
                          Observable_stat &energy) const {};
  virtual ParticleForce force(const Particle &p,
                              const Vector3d &folded_pos) = 0;

  virtual void reset_force(){};
};

} /* namespace Constraints */

#endif
