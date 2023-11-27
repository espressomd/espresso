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

namespace Coulomb {

struct Solver {
#ifdef ELECTROSTATICS
  struct Implementation;
  /// @brief Pointer-to-implementation.
  std::unique_ptr<Implementation> impl;
  /// @brief Whether to reinitialize the solver on observable calculation.
  bool reinit_on_observable_calc;

  Utils::Vector9d
  calc_pressure_long_range(ParticleRange const &particles) const;

  void sanity_checks() const;
  double cutoff() const;

  void on_observable_calc();
  void on_coulomb_change();
  void on_boxl_change();
  void on_node_grid_change();
  void on_periodicity_change();
  void on_cell_structure_change();
  void on_particle_change() { reinit_on_observable_calc = true; }

  void calc_long_range_force(ParticleRange const &particles) const;
  double calc_energy_long_range(ParticleRange const &particles) const;
  Solver();
#else  // ELECTROSTATICS
  Solver() = default;
#endif // ELECTROSTATICS

  using ShortRangeForceKernel =
      std::function<Utils::Vector3d(double, Utils::Vector3d const &, double)>;
  using ShortRangeForceCorrectionsKernel =
      std::function<void(Particle &, Particle &, double)>;
  using ShortRangePressureKernel = std::function<Utils::Matrix<double, 3, 3>(
      double, Utils::Vector3d const &, double)>;
  using ShortRangeEnergyKernel =
      std::function<double(Particle const &, Particle const &, double,
                           Utils::Vector3d const &, double)>;

  inline std::optional<ShortRangeForceKernel> pair_force_kernel() const;
  inline std::optional<ShortRangePressureKernel> pair_pressure_kernel() const;
  inline std::optional<ShortRangeEnergyKernel> pair_energy_kernel() const;
  inline std::optional<ShortRangeForceCorrectionsKernel>
  pair_force_elc_kernel() const;
};

#ifdef ELECTROSTATICS
Solver const &get_coulomb();
#endif

} // namespace Coulomb
