/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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
/** \file stretching_force.cpp
 *
 *  Implementation of \ref stretching_force.hpp
 */

#include "stretching_force.hpp"
#include "communication.hpp"

/// set the parameters for the stretching_force potential
int stretching_force_set_params(int bond_type, double r0, double ks)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.stretching_force.ks = ks;
  bonded_ia_params[bond_type].p.stretching_force.r0 = r0;
  
  bonded_ia_params[bond_type].type = BONDED_IA_STRETCHING_FORCE;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 
  
  return ES_OK;
}
