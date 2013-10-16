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
/** \file gb.cpp
 *
 *  Implementation of \ref gb.hpp
 */
#include "gb.hpp"
#include "communication.hpp"

#ifdef GAY_BERNE

int gay_berne_set_params(int part_type_a, int part_type_b,
			 double eps, double sig, double cut,
			 double k1, double k2,
			 double mu, double nu)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->GB_eps    = eps;
  data->GB_sig    = sig;
  data->GB_cut    = cut;
  data->GB_k1     = k1;
  data->GB_k2     = k2;
  data->GB_mu     = mu;
  data->GB_nu     = nu;
 
  /* Calculate dependent parameters */

  data->GB_chi1 = ((data->GB_k1*data->GB_k1) - 1) / ((data->GB_k1*data->GB_k1) + 1);
  data->GB_chi2 = (pow(data->GB_k2,(1/data->GB_mu))-1)/(pow(data->GB_k2,(1/data->GB_mu))+1);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif
