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
/** \file subt_lj.cpp
 *
 *  Implementation of \ref subt_lj.hpp
 */
#include "subt_lj.hpp"

#ifdef LENNARD_JONES
#include "communication.hpp"

int subt_lj_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return ES_ERROR;
  
  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.subt_lj.k = k;
  bonded_ia_params[bond_type].p.subt_lj.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_SUBT_LJ;  
  bonded_ia_params[bond_type].p.subt_lj.r2 = SQR(bonded_ia_params[bond_type].p.subt_lj.r);
  bonded_ia_params[bond_type].num = 1;

  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

#endif

