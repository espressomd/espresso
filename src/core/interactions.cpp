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

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "event.hpp"

#include "serialization/IA_parameters.hpp"

#include <boost/mpi.hpp>
#include <utils/mpi/cart_comm.hpp>

#include <mpi.h>

inline bool is_tabulated_bond(BondedInteraction const type) {
  return (type == BONDED_IA_TABULATED_DISTANCE or
          type == BONDED_IA_TABULATED_ANGLE or
          type == BONDED_IA_TABULATED_DIHEDRAL);
}

void mpi_bcast_all_ia_params_slave() {
  boost::mpi::broadcast(comm_cart, ia_params, 0);
}

REGISTER_CALLBACK(mpi_bcast_all_ia_params_slave)

void mpi_bcast_all_ia_params() { mpi_call_all(mpi_bcast_all_ia_params_slave); }

void mpi_bcast_ia_params_slave(int i, int j) {
  if (j >= 0) {
    // non-bonded interaction parameters
    boost::mpi::broadcast(comm_cart, *get_ia_param(i, j), 0);
  } else {
    // bonded interaction parameters
    if (this_node) {
      make_bond_type_exist(i); // realloc bonded_ia_params on slave nodes!
      if (is_tabulated_bond(bonded_ia_params[i].type)) {
        delete bonded_ia_params[i].p.tab.pot;
      }
    }
    MPI_Bcast(&(bonded_ia_params[i]), sizeof(Bonded_ia_parameters), MPI_BYTE, 0,
              comm_cart);
    // for tabulated potentials we have to send the tables extra
    if (is_tabulated_bond(bonded_ia_params[i].type)) {
      if (this_node) {
        bonded_ia_params[i].p.tab.pot = new TabulatedPotential();
      }
      boost::mpi::broadcast(comm_cart, *bonded_ia_params[i].p.tab.pot, 0);
    }
  }

  on_short_range_ia_change();
}

REGISTER_CALLBACK(mpi_bcast_ia_params_slave)

void mpi_bcast_ia_params(int i, int j) {
  mpi_call_all(mpi_bcast_ia_params_slave, i, j);
}

REGISTER_CALLBACK(realloc_ia_params)

void mpi_realloc_ia_params(int ns) { mpi_call_all(realloc_ia_params, ns); }
