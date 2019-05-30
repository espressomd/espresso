/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
class Constraint {
public:
  /**
   * @brief Add energy contribution of this constraints to energy.
   *
   * Add constraint energy for particle to observable.
   *
   * @param[in] p The particle to add the energy for.
   * @param[in] folded_pos Folded position of the particle.
   * @param[in] t The time at which the energy should be calculated.
   * @param[out] energy to add the energy to.
   */
  virtual void add_energy(const Particle &p, const Utils::Vector3d &folded_pos,
                          double t, Observable_stat &energy) const = 0;

  /**
   * @brief Calculate the force of the constraint on a particle.
   *
   * Add constraint energy for particle to observable.
   *
   * @param[in] p The particle to calculate the force for.
   * @param[in] folded_pos Folded position of the particle.
   * @param[in] t The time at which the force should be calculated.
   * @return The force on the particle.
   */
  virtual ParticleForce force(const Particle &p,
                              const Utils::Vector3d &folded_pos, double t) = 0;

  /**
   * @brief Check if constraints if compatible with box size.
   */
  virtual bool fits_in_box(Utils::Vector3d const &box) const = 0;

  virtual void reset_force(){};

  virtual ~Constraint() = default;
};
} /* namespace Constraints */

#endif
