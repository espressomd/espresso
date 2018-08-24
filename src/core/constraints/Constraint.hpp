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
