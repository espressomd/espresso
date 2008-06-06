// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#ifndef MOL_CUT_H
#define MOL_CUT_H

#include "virtual_sites.h"
#include "particle_data.h"
#include "interaction_data.h"
#include <tcl.h>

#ifdef MOL_CUT

#define CUTOFF_CHECK(cond) ((ia_params->mol_cut_type!=0)||(cond))

MDINLINE int checkIfParticlesInteractViaMolCut(Particle *p1, Particle *p2,IA_parameters *data){
   if (data->mol_cut_type==0){
      return 1;
   }
   else{
      double com_dist=get_mol_dist(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

MDINLINE int checkIfParticlesInteractViaMolCut_cfg(Particle *p1, Particle *p2,IA_parameters *data){
   if (data->mol_cut_type==0){
      return 1;
   }
   else{
      double com_dist=get_mol_dist_cfg(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

MDINLINE int printmolcutIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%i",data->mol_cut_type);
  Tcl_AppendResult(interp, "molcut ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, data->mol_cut_cutoff, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}

MDINLINE int molcut_set_params(int part_type_a, int part_type_b,int mol_cut_type,double mol_cut_cutoff)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  data->mol_cut_type         = data_sym->mol_cut_type         = mol_cut_type;
  data->mol_cut_cutoff         = data_sym->mol_cut_cutoff         = mol_cut_cutoff;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int molcut_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for molcut */
  int mol_cut_type;
  double mol_cut_cutoff;
  int change;

  /* get interaction type */
  if (argc < 3) {
    Tcl_AppendResult(interp, "molcut needs 2 parameter: "
		     "<type> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy parameters */
  if (! ARG_IS_I(1, mol_cut_type)) {
    Tcl_AppendResult(interp, "type must be int",
		     (char *) NULL);
    return TCL_ERROR;
  }
  if (! ARG_IS_D(2, mol_cut_cutoff)) {
    Tcl_AppendResult(interp, "cutoff must be double",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 3;
	
  if (! ((mol_cut_type==0) || (mol_cut_type==1)) ) {
    Tcl_AppendResult(interp, "type must be 0 or 1", (char *) NULL);
    return 0;
  }
  if (molcut_set_params(part_type_a, part_type_b,mol_cut_type,mol_cut_cutoff) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}
#else
#define CUTOFF_CHECK(cond) (cond)
#endif


#endif
