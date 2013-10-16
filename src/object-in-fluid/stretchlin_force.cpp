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
/** \file stretchlin_force.cpp
 *
 *  Implementation of \ref stretchlin_force.hpp
 */

#include "stretchlin_force.hpp"
#include "communication.hpp"

/// set the parameters for the stretchlin_force potential
int stretchlin_force_set_params(int bond_type, double r0, double kslin)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.stretchlin_force.kslin = kslin;
  bonded_ia_params[bond_type].p.stretchlin_force.r0 = r0;
  
  bonded_ia_params[bond_type].type = BONDED_IA_STRETCHLIN_FORCE;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 
  
  return ES_OK;
}
