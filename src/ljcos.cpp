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

#include "utils.hpp"
#include "ljcos.hpp"

#ifdef LJCOS
#include "communication.hpp"

int lj_cos_set_params(int part_type_a, int part_type_b,
		      double eps, double sig, double cut,
		      double offset)
{
  double facsq;
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);
  
  if (!data) return ES_ERROR;

  data->LJCOS_eps    = eps;
  data->LJCOS_sig    = sig;
  data->LJCOS_cut    = cut;
  data->LJCOS_offset = offset;

  /* Calculate dependent parameters */
  facsq = driwu2*SQR(sig);
  data->LJCOS_rmin = sqrt(driwu2)*sig;
  data->LJCOS_alfa = PI/(SQR(data->LJCOS_cut)-facsq);
  data->LJCOS_beta = PI*(1.-(1./(SQR(data->LJCOS_cut)/facsq-1.)));

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  
  return ES_OK;
}
#endif
