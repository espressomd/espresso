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
/** \file morse.cpp
 *
 *  Implementation of \ref morse.hpp
 */
#include "morse.hpp"
#include "communication.hpp"

#ifdef MORSE

int morse_set_params(int part_type_a, int part_type_b,
		     double eps, double alpha,
		     double rmin, double cut)
{
  double add1, add2;
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->MORSE_eps   = eps;
  data->MORSE_alpha = alpha;
  data->MORSE_rmin  = rmin;
  data->MORSE_cut   = cut;

  /* calculate dependent parameter */
  add1 = exp(-2.0*data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  add2 = 2.0*exp(-data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  data->MORSE_rest = data->MORSE_eps * (add1 - add2);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif
