/*
 * Copyright (C) 2010-2020 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXTBASE_HPP
#define SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXTBASE_HPP
/**
 * @file
 * @ref Scafacos::ScafacosContextBase provides the public interface
 * of the ScaFaCoS bridge. It relies on type erasure to hide the
 * ScaFaCoS implementation details from the ESPResSo core. It is
 * implemented by @ref Scafacos::ScafacosContext.
 */

#include "config.hpp"

#if defined(SCAFACOS)

#include <utils/Vector.hpp>

#include <string>

namespace Scafacos {

/**
 * @brief Public interface to the ScaFaCoS context.
 * Implementations of this class will store particle positions and
 * charges/dipoles in flat arrays, as required by ScaFaCoS, as well
 * as the output arrays.
 */
struct ScafacosContextBase {
  virtual ~ScafacosContextBase() = default;
  /** @brief Collect particle data in continuous arrays. */
  virtual void update_particle_data() = 0;
  /** @brief Write forces back to particles. */
  virtual void update_particle_forces() const = 0;
  /** @brief Calculate long-range part of the energy. */
  virtual double long_range_energy() = 0;
  /** @brief Add long-range part of the forces to particles. */
  virtual void add_long_range_force() = 0;
  /** @brief Add near-field pair force. */
  inline void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                             Utils::Vector3d &force) {
    if (dist > get_r_cut())
      return;

    auto const field = get_pair_force(dist);
    auto const fak = q1q2 * field / dist;
    force -= fak * d;
  }
  /** @brief Calculate near-field pair energy. */
  inline double pair_energy(double q1q2, double dist) {
    if (dist > get_r_cut())
      return 0.;

    return q1q2 * get_pair_energy(dist);
  }
  /** @brief Reinitialize number of particles, box shape and periodicity. */
  virtual void update_system_params() = 0;
  virtual double get_r_cut() const = 0;
  virtual double get_pair_force(double dist) const = 0;
  virtual double get_pair_energy(double dist) const = 0;
  virtual std::string get_method_and_parameters() = 0;
};

} // namespace Scafacos
#endif // SCAFACOS
#endif // SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXTBASE_HPP
