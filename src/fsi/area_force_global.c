/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

/** \file area_force_global.h
 *  Routines to calculate the AREA_FORCE_GLOBAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.c
*/

#include "fsi/area_force_global.h"
#include "communication.h"


/** set parameters for the AREA_FORCE_GLOBAL potential. 
*/
int area_force_global_set_params(int bond_type, double A0_g, double ka_g)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.area_force_global.ka_g = ka_g;
  bonded_ia_params[bond_type].p.area_force_global.A0_g = A0_g;

  bonded_ia_params[bond_type].type = BONDED_IA_AREA_FORCE_GLOBAL;
  bonded_ia_params[bond_type].num = 2;				
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

