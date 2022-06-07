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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_SCAFACOS_IMPL_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_SCAFACOS_IMPL_HPP

#include "config.hpp"

#ifdef SCAFACOS_DIPOLES

#include "magnetostatics/scafacos.hpp"

#include "scafacos/ScafacosContext.hpp"

#include <scafacos/Dipoles.hpp>

#include <boost/range/numeric.hpp>

#include <vector>

#if !defined(FCS_ENABLE_DIPOLES)
#error                                                                         \
    "SCAFACOS_DIPOLES requires dipoles support in scafacos library (FCS_ENABLE_DIPOLES)."
#endif

struct DipolarScafacosImpl : public DipolarScafacos,
                             public ScafacosContext<::Scafacos::Dipoles> {
  using ScafacosContext = ScafacosContext<::Scafacos::Dipoles>;
  using DipolarScafacos::DipolarScafacos;
  using ScafacosContext::ScafacosContext;

  void update_particle_data() override;
  void update_particle_forces() const override;

  double long_range_energy() override {
    update_particle_data();
    run(dipoles, positions, fields, potentials);
    return -0.5 * prefactor * boost::inner_product(dipoles, potentials, 0.0);
  }

  void add_long_range_forces() override {
    update_particle_data();
    run(dipoles, positions, fields, potentials);
    update_particle_forces();
  }

private:
  /** Inputs */
  std::vector<double> positions, dipoles;
  /** Outputs */
  std::vector<double> fields, potentials;
};

#endif // SCAFACOS_DIPOLES
#endif
