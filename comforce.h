// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file comforce.h
 *  Routines to enable comforce
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
*/
#include "utils.h"
#include "particle_data.h"
#include "statistics.h"
#include "parser.h"

#ifdef COMFORCE

MDINLINE int comforce_set_params(int part_type_a, int part_type_b,
				 int flag, int dir, double force, double fratio)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym)
    return 1;

  if (n_nodes > 1)
    return 2;


  /* COMFORCE should be symmetrically */
  data_sym->COMFORCE_flag    = data->COMFORCE_flag    = flag;
  data_sym->COMFORCE_dir    = data->COMFORCE_dir    = dir;
  data_sym->COMFORCE_force    = data->COMFORCE_force    = force;
  data_sym->COMFORCE_fratio    = data->COMFORCE_fratio    = fratio;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return 0;
}

MDINLINE int printcomforceIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%d",data->COMFORCE_flag);
  Tcl_AppendResult(interp, "comforce ", buffer, " ", (char *) NULL);
  sprintf(buffer,"%d",data->COMFORCE_dir);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COMFORCE_force, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->COMFORCE_fratio, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    
  return TCL_OK;
}

MDINLINE int comforce_parser(Tcl_Interp * interp,
			     int part_type_a, int part_type_b,
			     int argc, char ** argv)
{
  int flag, dir, change; 
  double force, fratio;

  if (argc != 5) {
    Tcl_AppendResult(interp, "comforce needs 4 parameters: "
		     "<comforce_flag> <comforce_dir> <comforce_force> <comforce_fratio>",
		     (char *) NULL);
    return 0;
  }
  
  if (part_type_a == part_type_b) {
    Tcl_AppendResult(interp, "comforce needs 2 different types ", (char *) NULL);
    return 0;
  }

  /* copy comforce parameters */
  if ((! ARG_IS_I(1, flag)) || (! ARG_IS_I(2, dir)) || (! ARG_IS_D(3, force)) || (! ARG_IS_D(4, fratio)) ) {
    Tcl_AppendResult(interp, "comforce needs 2 INTEGER 1 DOUBLE parameter: "
		     "<comforce_flag> <comforce_dir> <comforce_force> <comforce_fratio>", (char *) NULL);
    return 0;
  }
    
  change = 5;
    
  switch (comforce_set_params(part_type_a, part_type_b, flag, dir, force, fratio)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "works only with a single CPU", (char *) NULL);
    return 0;
  }
  return change;
}

MDINLINE void calc_comforce()
{
  int t0,t1,k, j;
  IA_parameters *ia_params;
  double com0[3], com1[3], MofImatrix[9], diff[3];
  double vect0[3], vect1[3], eva[3], eve[3], fvect[3];
  Particle *p;
  int i, np, c;
  Cell *cell;
  
  for (t0=0; t0<n_particle_types-1; t0++) {
    for (t1=t0+1; t1<n_particle_types; t1++) {
      ia_params = get_ia_param(t0,t1);
      if(ia_params->COMFORCE_flag == 1) {
	      centermass(t0,com0);
	      centermass(t1,com1);
	      for (i = 0; i < 3; i++) {
		      diff[i]=com1[i]-com0[i];
	      }
        momentofinertiamatrix(t0, MofImatrix);
        k=calc_eigenvalues_3x3(MofImatrix, eva);
        /* perpendicular force */
        if(ia_params->COMFORCE_dir == 1) {
	        k=calc_eigenvector_3x3(MofImatrix,eva[0],eve);
          /*By doing two vector products find radial axis along the target system */
	        vector_product(eve,diff,vect0);
	        vector_product(vect0,eve,vect1);
  
          /* normalize vect1, return is fvect */
	        unit_vector(vect1,fvect);
        } else {
        /* parallel force */
	        k=calc_eigenvector_3x3(MofImatrix,eva[0],fvect);
        }
        
        /* orient it along the com vector */
        if (scalar(fvect,diff) < 0.) {
	        for (i = 0; i < 3; i++) {
		        fvect[i] = -fvect[i];
	        }
        }
        
        /* Now apply the force */
        for (c = 0; c < local_cells.n; c++) {
          cell = local_cells.cell[c];
          p  = cell->part;
          np = cell->n;
          for(i = 0; i < np; i++) {
	 	        if(p[i].p.type==t0) {
      	      for(j = 0; j < 3; j++) {
				        p[i].f.f[j] -= ia_params->COMFORCE_fratio * ia_params->COMFORCE_force * fvect[j];
			        }
		        }
	 	        if(p[i].p.type==t1) {
      	      for(j = 0; j < 3; j++) {
				        p[i].f.f[j] +=  ia_params->COMFORCE_force * fvect[j];
			        }
		        }
          }
        }
        /*end of force application */
      }
    }
  }

}
#endif
