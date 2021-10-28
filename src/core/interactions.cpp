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
#include "communication.hpp"

#include "TabulatedPotential.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/bonded_tab.hpp"

#include "serialization/IA_parameters.hpp"

#include <boost/mpi.hpp>

void mpi_bcast_ia_params_local(int i, int j) {
  boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);
  on_short_range_ia_change();
}

REGISTER_CALLBACK(mpi_bcast_ia_params_local)

void mpi_bcast_ia_params(int i, int j) {
  mpi_call_all(mpi_bcast_ia_params_local, i, j);
}
