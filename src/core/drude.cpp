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
/** \file drude.cpp
 *
 *  Implementation of \ref drude.hpp
 */

#include "drude.hpp"
#include "communication.hpp"

#ifdef DRUDE

int drude_set_params(int bond_type, double temp_core, double gamma_core, double temp_drude, double gamma_drude, double k, double mass_drude, double r_cut)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.drude.temp_core = temp_core;
  bonded_ia_params[bond_type].p.drude.gamma_core = gamma_core;
  bonded_ia_params[bond_type].p.drude.temp_drude = temp_drude;
  bonded_ia_params[bond_type].p.drude.gamma_drude = gamma_drude;
  bonded_ia_params[bond_type].p.drude.k = k;
  bonded_ia_params[bond_type].p.drude.mass_drude = mass_drude;
  bonded_ia_params[bond_type].p.drude.r_cut = r_cut;

  bonded_ia_params[bond_type].type = BONDED_IA_DRUDE;

  bonded_ia_params[bond_type].num = 1;

  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

/*
void drude_recalc_params()
{



}
*/

#endif
