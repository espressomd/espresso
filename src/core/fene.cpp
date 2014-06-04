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
/** \file fene.cpp
 *
 *  Implementation of \ref fene.hpp
 */

#include "fene.hpp"
#include "communication.hpp"

/// set the parameters for the fene potential
int fene_set_params(int bond_type, double k, double drmax, double r0)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.fene.k = k;
  bonded_ia_params[bond_type].p.fene.drmax = drmax;
  bonded_ia_params[bond_type].p.fene.r0 = r0;

  bonded_ia_params[bond_type].p.fene.drmax2 = SQR(bonded_ia_params[bond_type].p.fene.drmax);
  bonded_ia_params[bond_type].p.fene.drmax2i = 1.0/bonded_ia_params[bond_type].p.fene.drmax2;

  bonded_ia_params[bond_type].type = BONDED_IA_FENE;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 
  
  return ES_OK;
}
