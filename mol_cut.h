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
   else if (p1->p.mol_id == p2->p.mol_id) {
      return 1;
   }
   else{
      double com_dist=get_mol_dist(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

MDINLINE int checkIfParticlesInteractViaMolCut_partcfg(Particle *p1, Particle *p2,IA_parameters *data){
   if (data->mol_cut_type==0){
      return 1;
   }
   else if (p1->p.mol_id == p2->p.mol_id) {
      return 1;
   }
   else{
      double com_dist=get_mol_dist_partcfg(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

MDINLINE int tclprint_to_result_molcutIA(Tcl_Interp *interp, int i, int j)
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

MDINLINE int tclcommand_inter_parse_molcut(Tcl_Interp * interp,
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
