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

#include "actor/traits.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <functional>
#include <memory>
#include <type_traits>

#ifdef ELECTROSTATICS
// forward declarations
struct DebyeHueckel;
struct ReactionField;
struct ICCStar;
#ifdef P3M
struct CoulombP3M;
#ifdef CUDA
struct CoulombP3MGPU;
#endif
struct ElectrostaticLayerCorrection;
#endif
struct CoulombMMM1D;
#ifdef MMM1D_GPU
class CoulombMMM1DGpu;
#endif
#ifdef SCAFACOS
struct CoulombScafacos;
#endif

using ElectrostaticsActor =
    boost::variant<std::shared_ptr<DebyeHueckel>,
#ifdef P3M
                   std::shared_ptr<CoulombP3M>,
#ifdef CUDA
                   std::shared_ptr<CoulombP3MGPU>,
#endif // CUDA
                   std::shared_ptr<ElectrostaticLayerCorrection>,
#endif // P3M
                   std::shared_ptr<CoulombMMM1D>,
#ifdef MMM1D_GPU
                   std::shared_ptr<CoulombMMM1DGpu>,
#endif // MMM1D_GPU
#ifdef SCAFACOS
                   std::shared_ptr<CoulombScafacos>,
#endif // SCAFACOS
                   std::shared_ptr<ReactionField>>;

using ElectrostaticsExtension = boost::variant<std::shared_ptr<ICCStar>>;
#endif // ELECTROSTATICS

namespace Coulomb {

namespace detail {
bool flag_all_reduce(bool flag);
} // namespace detail

struct Solver {
#ifdef ELECTROSTATICS
  /// @brief Main electrostatics solver.
  boost::optional<ElectrostaticsActor> solver;
  /// @brief Extension that modifies the solver behavior.
  boost::optional<ElectrostaticsExtension> extension;
  /// @brief Whether to reinitialize the solver on observable calculation.
  bool reinit_on_observable_calc = false;

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

  inline boost::optional<ShortRangeForceKernel> pair_force_kernel() const;
  inline boost::optional<ShortRangePressureKernel> pair_pressure_kernel() const;
  inline boost::optional<ShortRangeEnergyKernel> pair_energy_kernel() const;
  inline boost::optional<ShortRangeForceCorrectionsKernel>
  pair_force_elc_kernel() const;
};

#ifdef ELECTROSTATICS
Solver const &get_coulomb();
#endif

/** @brief Check if the system is charge-neutral. */
void check_charge_neutrality(double relative_tolerance);

} // namespace Coulomb
