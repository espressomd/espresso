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

/** \file oif_local_forces.hpp
 *  Routines to calculate the OIF_LOCAL_FORCES
 *  for a particle quadruple (two neighboring triangles with common edge). (Dupin2007)
 *  \ref forces.cpp
 */

#include "communication.hpp"
#include "oif_local_forces.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

/** set parameters for the OIF_LOCAL_FORCES potential.*/

int oif_local_forces_set_params(int bond_type, double r0, double ks, double kslin, double phi0, double kb, double A01, double A02, double kal)
{
  if(bond_type < 0)
    return ES_ERROR;

  //create bond class
  bond_container.set_bond_by_type(bond_type,
				  Utils::make_unique<Bond::OifLocalForces>(phi0, kb, r0, ks, kslin,
									   A01, A02, kal));

  return ES_OK;
}

