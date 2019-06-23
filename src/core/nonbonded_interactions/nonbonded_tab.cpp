/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *
 *  Implementation of \ref nonbonded_tab.hpp
 */
#include "nonbonded_interactions/nonbonded_tab.hpp"

#ifdef TABULATED
#include "communication.hpp"

#include <utils/constants.hpp>

int tabulated_set_params(int part_type_a, int part_type_b, double min,
                         double max, std::vector<double> const &energy,
                         std::vector<double> const &force) {
  auto data = get_ia_param_safe(part_type_a, part_type_b);
  assert(max >= min);
  assert((max == min) || force.size() > 1);
  assert(force.size() == energy.size());

  data->TAB.maxval = max;
  data->TAB.minval = min;
  if (max == min)
    data->TAB.invstepsize = 0;
  else
    data->TAB.invstepsize = static_cast<double>(force.size() - 1) / (max - min);

  data->TAB.force_tab = force;
  data->TAB.energy_tab = energy;

  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif
