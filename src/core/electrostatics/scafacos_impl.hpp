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

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_SCAFACOS_IMPL_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_SCAFACOS_IMPL_HPP

#include "config.hpp"

#ifdef SCAFACOS

#include "electrostatics/scafacos.hpp"

#include "scafacos/ScafacosContext.hpp"

#include <scafacos/Coulomb.hpp>

#include <boost/range/numeric.hpp>

#include <vector>

struct CoulombScafacosImpl : public CoulombScafacos,
                             public ScafacosContext<::Scafacos::Coulomb> {
  using ScafacosContext = ScafacosContext<::Scafacos::Coulomb>;
  using CoulombScafacos::CoulombScafacos;
  using ScafacosContext::ScafacosContext;

  void update_particle_data() override;
  void update_particle_forces() const override;

  double long_range_energy() override {
    update_particle_data();
    run(charges, positions, fields, potentials);
    return 0.5 * prefactor * boost::inner_product(charges, potentials, 0.);
  }

  void add_long_range_forces() override {
    update_particle_data();
    run(charges, positions, fields, potentials);
    update_particle_forces();
  }

  double get_pair_force(double dist) const override {
    return ScafacosContext::pair_force(dist);
  }

  double get_pair_energy(double dist) const override {
    return ScafacosContext::pair_energy(dist);
  }

  void set_near_field_delegation(bool delegate) override {
    return ScafacosContext::set_near_field_delegation(delegate);
  }

  bool get_near_field_delegation() const override {
    return ScafacosContext::get_near_field_delegation();
  }

  double get_r_cut() const override { return ScafacosContext::r_cut(); }

private:
  /** Inputs */
  std::vector<double> positions, charges;
  /** Outputs */
  std::vector<double> fields, potentials;

  /** Determine runtime for a specific cutoff */
  double time_r_cut(double r_cut);

  /** Determine the optimal cutoff by bisection */
  void tune_r_cut();
  void tune_impl() override;
  void set_r_cut_and_tune(double r_cut) {
    update_particle_data();
    set_r_cut(r_cut);
    ScafacosContext::tune(charges, positions);
  }
};

#endif // SCAFACOS
#endif
