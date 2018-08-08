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

/** \file ljcos2.hpp
 *
 *  Routines to calculate the lennard-jones with cosine tail energy and/or  force 
 *  for a particle pair.  Cosine tail is different from that in ljcos.hpp
 *  Used for attractive tail/tail interactions in lipid bilayer calculations
 *  \ref forces.cpp
*/
#include "ljcos2.hpp"

#ifdef LJCOS2
#include <cmath>

#include "communication.hpp"

int ljcos2_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double offset,
				      double w)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->LJCOS2_eps    = eps;
  data->LJCOS2_sig    = sig;
  data->LJCOS2_offset = offset;
  data->LJCOS2_w      = w;

  /* calculate dependent parameters */
  data->LJCOS2_rchange = pow(2,1/6.)*sig;
  data->LJCOS2_cut     = w + data->LJCOS2_rchange;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif /* ifdef LJCOS2 */
