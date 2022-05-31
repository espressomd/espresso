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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_SCAFACOS_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_SCAFACOS_HPP

#include "config.hpp"

#ifdef SCAFACOS_DIPOLES

#include "scafacos/ScafacosContextBase.hpp"

#include <memory>
#include <string>

struct DipolarScafacos : virtual public ScafacosContextBase {
  ~DipolarScafacos() override = default;
  double prefactor = 0.;

  void on_activation() { update_system_params(); }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { update_system_params(); }
  void on_node_grid_change() const {}
  void on_periodicity_change() { update_system_params(); }
  void on_cell_structure_change() const {}
  void init() const {}

  void sanity_checks() const override {}
};

std::shared_ptr<DipolarScafacos>
make_dipolar_scafacos(std::string const &method, std::string const &parameters);

#endif // SCAFACOS_DIPOLES
#endif
