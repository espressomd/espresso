/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#include "thermalized_bond_utils.hpp"
#include "integrate.hpp"
#include "thermalized_bond.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"

#include <boost/variant.hpp>

void thermalized_bond_init(double time_step) {
  for (auto &kv : bonded_ia_params) {
    if (auto *t = boost::get<ThermalizedBond>(&(*kv.second))) {
      t->pref1_com = t->gamma_com;
      t->pref2_com = sqrt(24.0 * t->gamma_com / time_step * t->temp_com);
      t->pref1_dist = t->gamma_distance;
      t->pref2_dist =
          sqrt(24.0 * t->gamma_distance / time_step * t->temp_distance);
    }
  }
}
