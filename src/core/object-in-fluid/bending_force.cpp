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

/** \file bending_force.hpp Routines to calculate the bending_force energy or/and
 *  and force for a particle quadruple (two triangles that have 2 particles in common)
*/

#include "communication.hpp"
#include "bending_force.hpp"

/// set bending_force parameters
int bending_force_set_params(int bond_type, double phi0, double kb)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_BENDING_FORCE;
  bonded_ia_params[bond_type].num  = 3;
  bonded_ia_params[bond_type].p.bending_force.phi0 = phi0;
  bonded_ia_params[bond_type].p.bending_force.kb = kb;

  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

