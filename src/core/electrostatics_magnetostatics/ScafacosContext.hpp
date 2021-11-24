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
/**
 * @file
 * @ref Scafacos::ScafacosContextBase implements the interface of the
 * ScaFaCoS bridge. It is further derived for the coulombic and dipolar
 * versions of ScaFaCoS.
 */

#include "config.hpp"

#if defined(SCAFACOS)

#include "Scafacos.hpp"
#include "electrostatics_magnetostatics/ScafacosContextBase.hpp"

#include <utils/Vector.hpp>

#include <string>
#include <vector>

#if defined(SCAFACOS_DIPOLES) && !defined(FCS_ENABLE_DIPOLES)
#error                                                                         \
    "SCAFACOS_DIPOLES requires dipoles support in scafacos library (FCS_ENABLE_DIPOLES)."
#endif

namespace Scafacos {

/** Encapsulation for the particle data needed by ScaFaCoS */
struct ScafacosContext : ScafacosContextBase, Scafacos {
  using Scafacos::Scafacos;
  using ScafacosContextBase::ScafacosContextBase;
  ~ScafacosContext() override = default;
  void update_system_params() override;
  double get_pair_force(double dist) const override {
    return Scafacos::pair_force(dist);
  }
  double get_pair_energy(double dist) const override {
    return Scafacos::pair_energy(dist);
  }
  std::string get_method_and_parameters() override;
  double get_r_cut() const override { return Scafacos::r_cut(); }

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
  void tune();
  void set_r_cut_and_tune(double r_cut);

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
