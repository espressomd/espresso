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
/** \file bmhtf-nacl.cpp
 *
 *  Implementation of \ref bmhtf-nacl.hpp
 */
#include "bmhtf-nacl.hpp"

#ifdef BMHTF_NACL
#include "communication.hpp"

int BMHTF_set_params(int part_type_a, int part_type_b,
		     double A, double B, double C,
		     double D, double sig, double cut)
{
  double shift, dist2, pw6;
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  dist2 = cut*cut;
  pw6 = dist2*dist2*dist2;
  shift = -(A*exp(B*(sig - cut)) - C/pw6 - D/pw6/dist2);

  data->BMHTF_A   = A;
  data->BMHTF_B   = B;
  data->BMHTF_C   = C;
  data->BMHTF_D   = D;
  data->BMHTF_sig = sig;
  data->BMHTF_cut = cut;
  data->BMHTF_computed_shift = shift;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif
