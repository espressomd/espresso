#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:
  /**
   * @brief Add energy contribution of this constraints to energy.
   */
  virtual void add_energy(const Particle &p, const Vector3d &folded_pos,
                          Observable_stat &energy) const = 0;
  /**
   * @brief Return constraint force on particle.
   */
  virtual ParticleForce force(const Particle &p,
                              const Vector3d &folded_pos) = 0;

  /**
   * @brief Check if constraints if compatible with box size.
   */
  virtual bool fits_in_box(Vector3d const &box) const = 0;

  virtual void reset_force(){};

  virtual ~Constraint() = default;
};
} /* namespace Constraints */

#endif
