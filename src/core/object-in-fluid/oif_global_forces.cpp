/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
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

/** \file oif_global_forces.hpp
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "oif_global_forces.hpp"
#include "communication.hpp"


/** set parameters for the OIF_GLOBAL_FORCES potential. 
*/
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g, double V0, double kv)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.oif_global_forces.ka_g = ka_g;
  bonded_ia_params[bond_type].p.oif_global_forces.A0_g = A0_g;
  bonded_ia_params[bond_type].p.oif_global_forces.V0 = V0;
  bonded_ia_params[bond_type].p.oif_global_forces.kv = kv;

  bonded_ia_params[bond_type].type = BONDED_IA_OIF_GLOBAL_FORCES;
  bonded_ia_params[bond_type].num = 2;				
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

