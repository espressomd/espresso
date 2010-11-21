/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef COMFIXED_H
#define COMFIXED_H

#include "utils.h"

/** \file comfixed.h
 *  Routines to enable comfixed
 */
#ifdef COMFIXED

MDINLINE int comfixed_set_params(int part_type_a, int part_type_b, int flag)
{
  Particle *p;
  int i, j, np, c;
  Cell *cell;
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym)
    return 1;

  if (n_nodes > 1)
    return 2;

  if (PERIODIC(0) || PERIODIC(1) || PERIODIC(2)) {
    return 3;
  }

  /* COMFIXED should be symmetrically */
  data_sym->COMFIXED_flag    = data->COMFIXED_flag    = flag;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if(p[i].p.type==part_type_a) {
	for(j = 0; j < 3; j++) {
	  p[i].m.v[j] = 0.;
	  p[i].f.f[j] = 0.;
	}
      }
    }
  }
  
  return 0;
}

MDINLINE int tclprint_to_result_comfixedIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%d",data->COMFIXED_flag);
  Tcl_AppendResult(interp, "comfixed ", buffer, " ", (char *) NULL);
    
  return TCL_OK;
}

MDINLINE int tclcommand_inter_parse_comfixed(Tcl_Interp * interp,
			     int part_type_a, int part_type_b,
			     int argc, char ** argv)
{
  int flagc;
  	
  if (argc != 2) {
    Tcl_AppendResult(interp, "comfixed needs 1 parameters: "
		     "<comfixed_flag> ", (char *) NULL);
    return 0;
  }
	 
  if (part_type_a != part_type_b) {
    Tcl_AppendResult(interp, "comfixed must be among same type interactions", (char *) NULL);
    return 0;
  }

  /* copy comfixed parameters */
  if ((! ARG_IS_I(1, flagc)) ) {
    Tcl_AppendResult(interp, "comfixed needs 1 INTEGER parameter: "
		     "<comfixed_flag>", (char *) NULL);
    return 0;
  }

  switch (comfixed_set_params(part_type_a, part_type_b, flagc)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "works only with a single CPU", (char *) NULL);
    return 0;
  case 3:
    Tcl_AppendResult(interp, "works only with non periodic BC", (char *) NULL);
    return 0;
  }

   
   return 2;
}

MDINLINE void calc_comfixed()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  IA_parameters *ia_params;
  int t0;
  int j;
  double fsum0[3], type_mass;

  for (t0=0; t0<n_particle_types; t0++) {
    ia_params = get_ia_param(t0,t0);
    if(ia_params->COMFIXED_flag == 1) {
      type_mass=0.0;
      for(j = 0; j < 3; j++) {
	fsum0[j]= 0.;
      }
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
	  if(p[i].p.type==t0) {
	    type_mass += PMASS(p[i]);
      	    for(j = 0; j < 3; j++) {
	      fsum0[j] += p[i].f.f[j];
	    }
	  }
        }
      }
      
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
	  if(p[i].p.type==t0) {
      	    for(j = 0; j < 3; j++) {
	      p[i].f.f[j] -= fsum0[j]/type_mass*PMASS(p[i]);
	    }
	  }
        }
      }
    }
  }  
}
#endif

#endif
