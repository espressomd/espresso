#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:

  Constraint() {
    reset_force();
  }

  virtual void add_energy(Particle *p, double *folded_pos,
                  Observable_stat &energy) const = 0;

  virtual void add_force(Particle *p, double *folded_pos) = 0;

  void reset_force() { m_local_force = Vector3d{0, 0, 0}; }
  Vector3d total_force() const;

private:
  Vector3d m_local_force;
};

} /* namespace Constaints */

#endif
