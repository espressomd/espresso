/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file thermalized_bond.cpp
 *
 *  Implementation of \ref thermalized_bond.hpp
 */

#include "thermalized_bond.hpp"
#include "communication.hpp"
#include "global.hpp"
#include "interaction_data.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

int n_thermalized_bonds = 0;

int thermalized_bond_set_params(int bond_type, double temp_com, double gamma_com, double temp_distance, double gamma_distance, double r_cut)
{
  if (bond_type < 0)
    return ES_ERROR;

  n_thermalized_bonds += 1;

  //create bond class
  bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::ThermalizedBond>
				  (temp_com, gamma_com, temp_distance, gamma_distance, r_cut,
				   gamma_com/time_step, sqrt(24.0*gamma_com/time_step*temp_com),
				   gamma_distance/time_step,
				   sqrt(24.0*gamma_distance/time_step*temp_distance)));
  
  mpi_bcast_parameter(FIELD_THERMALIZEDBONDS);
  return ES_OK;
}

void thermalized_bond_heat_up() {
  double pref_scale = sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void thermalized_bond_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void thermalized_bond_init()
{

  bond_container.thermalized_bond_init();
  
}


void thermalized_bond_update_params(double pref_scale) {
  
  bond_container.thermalized_bond_update_params(pref_scale);
  
}
