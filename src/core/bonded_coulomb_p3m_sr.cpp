/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file bonded_coulomb.cpp
 *
 *  Implementation of \ref bonded_coulomb.hpp
 */

#include "bonded_coulomb_p3m_sr.hpp"
#include "communication.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object

#ifdef ELECTROSTATICS

int bonded_coulomb_p3m_sr_set_params(int bond_type, double q1q2)
{
  if(bond_type < 0)
    return ES_ERROR;

  //create bond class
  bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::BondedCoulombP3MSR>(q1q2));

  return ES_OK;
}

#endif
