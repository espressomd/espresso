/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#include "actor/optional.hpp"
#include "actor/traits.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <type_traits>

namespace Dipoles {

struct Solver {
#ifdef DIPOLES
  struct Implementation;
  /// @brief Pointer-to-implementation.
  std::unique_ptr<Implementation> impl;
  /// @brief Whether to reinitialize the solver on observable calculation.
  bool reinit_on_observable_calc;

  void sanity_checks() const;
  double cutoff() const;

  void on_observable_calc();
  void on_dipoles_change();
  void on_boxl_change();
  void on_node_grid_change();
  void on_periodicity_change();
  void on_cell_structure_change();
  void on_particle_change() { reinit_on_observable_calc = true; }

  void calc_pressure_long_range() const;
  void calc_long_range_force(ParticleRange const &particles) const;
  double calc_energy_long_range(ParticleRange const &particles) const;
#ifdef DIPOLE_FIELD_TRACKING
  void calc_long_range_field(ParticleRange const &particles) const;
#endif
  Solver();
#else  // DIPOLES
  Solver() = default;
#endif // DIPOLES

  using ShortRangeForceKernel =
      std::function<ParticleForce(Particle const &, Particle const &,
                                  Utils::Vector3d const &, double, double)>;
  using ShortRangeEnergyKernel =
      std::function<double(Particle const &, Particle const &,
                           Utils::Vector3d const &, double, double)>;

  inline std::optional<ShortRangeForceKernel> pair_force_kernel() const;
  inline std::optional<ShortRangeEnergyKernel> pair_energy_kernel() const;
};

#ifdef DIPOLES
Solver const &get_dipoles();
#endif

} // namespace Dipoles
