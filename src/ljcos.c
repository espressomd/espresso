/*
  Copyright (C) 2010,2011 The ESPResSo project
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

#include "utils.h"
#include "ljcos.h"

#ifdef LJCOS
#include "parser.h"
#include "communication.h"

int lj_cos_set_params(int part_type_a, int part_type_b,
		      double eps, double sig, double cut,
		      double offset)
{
  IA_parameters *data, *data_sym;

  double facsq;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJCOS should be symmetrically */
  data_sym->LJCOS_eps    = data->LJCOS_eps    = eps;
  data_sym->LJCOS_sig    = data->LJCOS_sig    = sig;
  data_sym->LJCOS_cut    = data->LJCOS_cut    = cut;
  data_sym->LJCOS_offset = data->LJCOS_offset = offset;

  /* Calculate dependent parameters */
  facsq = driwu2*SQR(sig);
  data_sym->LJCOS_rmin = data->LJCOS_rmin = sqrt(driwu2)*sig;
  data_sym->LJCOS_alfa = data->LJCOS_alfa = PI/(SQR(data->LJCOS_cut)-facsq);
  data_sym->LJCOS_beta = data->LJCOS_beta = PI*(1.-(1./(SQR(data->LJCOS_cut)/facsq-1.)));

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);
  
  return TCL_OK;
}
#endif
