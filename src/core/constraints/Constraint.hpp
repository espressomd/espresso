#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:
  virtual void add_energy(const Particle &p, const Vector3d &folded_pos,
                          Observable_stat &energy) const = 0;
  virtual ParticleForce force(const Particle &p,
                              const Vector3d &folded_pos) = 0;

  virtual void reset_force(){};

  virtual ~Constraint() {}
};

} /* namespace Constraints */

#endif
