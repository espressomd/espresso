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
#ifndef SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXT_HPP
#define SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXT_HPP
/** @file
 *  Provide a C-like interface for ScaFaCoS.
 */

#include "config.hpp"

#if defined(SCAFACOS)

#include "Scafacos.hpp"

#include <utils/Vector.hpp>

#include <string>
#include <vector>

#if defined(SCAFACOS_DIPOLES) && !defined(FCS_ENABLE_DIPOLES)
#error                                                                         \
    "SCAFACOS_DIPOLES requires dipoles support in scafacos library (FCS_ENABLE_DIPOLES)."
#endif

namespace Scafacos {

/** Encapsulation for the particle data needed by ScaFaCoS */
struct ScafacosContext : Scafacos {
  using Scafacos::Scafacos;
  virtual ~ScafacosContext() {}
  /** @brief Collect particle data in continuous arrays as required
   *  by ScaFaCoS.
   */
  virtual void update_particle_data() = 0;
  /** @brief Write forces back to particles. */
  virtual void update_particle_forces() const = 0;
  virtual double long_range_energy() = 0;
  virtual void add_long_range_force() = 0;
  /** @brief Reinitialize number of particles, box shape and periodicity. */
  void update_system_params();
  void add_pair_force(double q1q2, Utils::Vector3d const &d, double dist,
                      Utils::Vector3d &force) {
    if (dist > r_cut())
      return;

    auto const field = Scafacos::pair_force(dist);
    auto const fak = q1q2 * field / dist;
    force -= fak * d;
  }
  double pair_energy(double q1q2, double dist) {
    if (dist > r_cut())
      return 0.;

    return q1q2 * Scafacos::pair_energy(dist);
  }
  std::string get_method_and_parameters();

protected:
  /** Outputs */
  std::vector<double> fields, potentials;
};

struct ScafacosContextCoulomb : ScafacosContext {
  using ScafacosContext::ScafacosContext;
  void update_particle_data() override;
  void update_particle_forces() const override;
  double long_range_energy() override;
  void add_long_range_force() override {
    update_particle_data();
    run(charges, positions, fields, potentials);
    update_particle_forces();
  }
  void tune() { Scafacos::tune(charges, positions); }

private:
  /** Inputs */
  std::vector<double> positions, charges;
};

#ifdef SCAFACOS_DIPOLES
struct ScafacosContextDipoles : ScafacosContext {
  using ScafacosContext::ScafacosContext;
  void update_particle_data() override;
  void update_particle_forces() const override;
  double long_range_energy() override;
  void add_long_range_force() override {
    update_particle_data();
    run_dipolar(dipoles, positions, fields, potentials);
    update_particle_forces();
  }

private:
  /** Inputs */
  std::vector<double> positions, dipoles;
};
#endif

} // namespace Scafacos
#endif // SCAFACOS
#endif // SRC_CORE_ELECTROSTATICS_MAGNETOSTATICS_SCAFACOSCONTEXT_HPP
