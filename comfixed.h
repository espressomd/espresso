// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.

#include "utils.h"

/** \file comfixed.h
 *  Routines to enable comfixed
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
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
  
  if (!data || !data_sym) {
    return TCL_ERROR;
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
  
  return TCL_OK;
}

MDINLINE int printcomfixedIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

    sprintf(buffer,"%d",data->COMFIXED_flag);
    Tcl_AppendResult(interp, "comfixed ", buffer, " ", (char *) NULL);
    
    return TCL_OK;
}

MDINLINE int comfixed_parser(Tcl_Interp * interp,
			   int part_type_a, int part_type_b,
			   int argc, char ** argv, int *change)
{
  int flagc;
  	
  if (argc != 2) {
    Tcl_AppendResult(interp, "comfixed needs 1 parameters: "
		     "<comfixed_flag> ", (char *) NULL);
    return TCL_ERROR;
  }
	 
	if (part_type_a != part_type_b) {
	  Tcl_AppendResult(interp, "comfixed must be among same type interactions", (char *) NULL);
	  return TCL_ERROR;
	}

  /* copy comfixed parameters */
  if ((! ARG_IS_I(1, flagc)) )
  {
	  Tcl_AppendResult(interp, "comfixed needs 1 INTEGER parameter: "
			"<comfixed_flag>", (char *) NULL);
	  return TCL_ERROR;
  }

  *change = 2;

  if (comfixed_set_params(part_type_a, part_type_b, flagc) == TCL_ERROR) {
	  Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
	  return TCL_ERROR;
  }
}

MDINLINE void calc_comfixed()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  IA_parameters *ia_params;
  int t0,type_count0;
  int j;
	double fsum0[3];

  for (t0=0; t0<n_particle_types; t0++) {
    ia_params = get_ia_param(t0,t0);
    if(ia_params->COMFIXED_flag == 1) {
  	  type_count0=0;
      for(j = 0; j < 3; j++) {
			  fsum0[j]= 0.;
			}
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
	 	      if(p[i].p.type==t0) {
			      type_count0 ++;
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
				      p[i].f.f[j] -= fsum0[j]/type_count0;
			      }
		      }
        }
      }
    }
  }  
}
#endif
