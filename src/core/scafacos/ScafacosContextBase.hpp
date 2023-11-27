/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

/**
 * @file
 * @ref ScafacosContextBase provides the public interface
 * of the ScaFaCoS bridge. It relies on type erasure to hide the
 * ScaFaCoS implementation details from the ESPResSo core. It is
 * implemented by @ref ScafacosContext.
 */

#pragma once

#include "config/config.hpp"

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLAR)

#include <utils/Vector.hpp>

#include <string>
#include <vector>

/**
 * @brief Public interface to the ScaFaCoS context.
 * Implementations of this class will store particle positions and
 * charges/dipoles in flat arrays, as required by ScaFaCoS, as well
 * as the output arrays.
 */
struct ScafacosContextBase {
  ScafacosContextBase() = default;
  virtual ~ScafacosContextBase() = default;
  /** @brief Collect particle data in continuous arrays. */
  virtual void update_particle_data() = 0;
  /** @brief Write forces back to particles. */
  virtual void update_particle_forces() const = 0;
  /** @brief Calculate long-range part of the energy. */
  virtual double long_range_energy() = 0;
  /** @brief Add long-range part of the forces to particles. */
  virtual void add_long_range_forces() = 0;
  /** @brief Reinitialize number of particles, box shape and periodicity. */
  virtual void update_system_params() = 0;
  virtual std::string get_method() const = 0;
  virtual std::string get_parameters() const = 0;
  virtual void sanity_checks() const = 0;

  static std::vector<std::string> available_methods();
  static void sanity_check_method(std::string const &method_name);
};

#endif // SCAFACOS or SCAFACOS_DIPOLAR
