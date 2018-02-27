/*
  Copyright (C) 2012,2013,2016 The ESPResSo project
  
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

/** \file out_direction.hpp Routines to calculate the outward direction of the membrane
 *  using a particle quadruple (one particle and its 3 strategically placed neighbors)
*/

#include "../communication.hpp"
#include "out_direction.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

#ifdef MEMBRANE_COLLISION

// set out_direction parameters
int out_direction_set_params(int bond_type)
{
  if(bond_type < 0)
    return ES_ERROR;

  //create bond classs
  bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::MembraneCollision>());

  return ES_OK;
}

#endif
