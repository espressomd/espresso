/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/** \file
 *
 *  Implementation of \ref thermalized_bond.hpp
 */

#include "thermalized_bond.hpp"
#include "event.hpp"

int n_thermalized_bonds = 0;

ThermalizedBond::ThermalizedBond(double temp_com, double gamma_com,
                                 double temp_distance, double gamma_distance,
                                 double r_cut) {
  this->temp_com = temp_com;
  this->gamma_com = gamma_com;
  this->temp_distance = temp_distance;
  this->gamma_distance = gamma_distance;
  this->r_cut = r_cut;

  pref1_com = -1.;
  pref2_com = -1.;
  pref1_dist = -1.;
  pref2_dist = -1.;

  n_thermalized_bonds++;
  on_thermostat_param_change();
}
