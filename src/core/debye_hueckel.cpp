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
/** \file debye_hueckel.cpp
 *
 *  Implementation of \ref debye_hueckel.hpp
 */
#include "utils.hpp"
#include "debye_hueckel.hpp"
#include "communication.hpp"

#ifdef ELECTROSTATICS

int dh_set_params(double kappa, double r_cut)
{
  if(dh_params.kappa < 0.0)
    return -1;

  if(dh_params.r_cut < 0.0)
    return -2;

  dh_params.kappa = kappa;
  dh_params.r_cut = r_cut;

  mpi_bcast_coulomb_params();

  return 1;
}

#endif
