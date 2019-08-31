/*
  Copyright (C) 2012-2018 The ESPResSo project

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
 * Routines to calculate the outward direction of the
 * membrane using a particle quadruple (one particle and its 3 strategically
 * placed neighbors)
 */

#include "out_direction.hpp"

#ifdef MEMBRANE_COLLISION
#include "communication.hpp"

#include <utils/constants.hpp>

// set out_direction parameters
int oif_out_direction_set_params(int bond_type) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_OIF_OUT_DIRECTION;
  bonded_ia_params[bond_type].num = 3;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

#endif
