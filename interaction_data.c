/** \file interaction_data.c
    Implementation of \ref interaction_data.h "interaction_data.h"
 */
#include <string.h>
#include <stdlib.h>
#include "interaction_data.h"
#include "debug.h"
#include "communication.h"
#include "p3m.h"

/****************************************
 * variables
 *****************************************/

IA_parameters *ia_params = NULL;
int n_particle_types = 0;

int n_bonded_ia = 0;
Bonded_ia_parameters *bonded_ia_params = NULL;

double max_cut;

double lj_force_cap = 0.0;

/*****************************************
 * functions
 *****************************************/

/** Initialize interaction parameters. */
void initialize_ia_params(IA_parameters *params) {
  params->LJ_eps =
    params->LJ_sig =
    params->LJ_cut =
    params->LJ_shift =
    params->LJ_offset = 
    params->LJ_capradius = 0;
  
  params->ramp_cut =
    params->ramp_force = 0;
}

/** Copy interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
  dst->LJ_eps = src->LJ_eps;
  dst->LJ_sig = src->LJ_sig;
  dst->LJ_cut = src->LJ_cut;
  dst->LJ_shift = src->LJ_shift;
  dst->LJ_offset = src->LJ_offset;
  dst->LJ_capradius = src->LJ_capradius;

  dst->ramp_cut = src->ramp_cut;
  dst->ramp_force = src->ramp_force;  
}

/** returns non-zero if particles of type i and j have a nonbonded interaction */
int checkIfParticlesInteract(int i, int j) {
  IA_parameters *data = get_ia_param(i, j);

  if (data->LJ_eps != 0)
    return 1;
  
  if (data->ramp_cut > 0)
    return 1;

  return 0;
}

/** This function increases the LOCAL ia_params field
    to the given size. This function is not exported
    since it does not do this on all nodes. Use
    make_particle_type_exist for that.
*/
void realloc_ia_params(int nsize)
{
  int i, j;
  IA_parameters *new_params;

  if (nsize <= n_particle_types)
    return;

  new_params = (IA_parameters *) malloc(nsize*nsize*sizeof(IA_parameters));
  if (ia_params) {
    /* if there is an old field, copy entries and delete */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++) {
	if ((i < n_particle_types) && (j < n_particle_types))
	  copy_ia_params(&new_params[i*nsize + j],
			 &ia_params[i*n_particle_types + j]);
	else
	  initialize_ia_params(&new_params[i*nsize + j]);
      }
    free(ia_params);
  }
  else {
    /* new field, just init */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++)
	initialize_ia_params(&new_params[i*nsize + j]);
  }

  n_particle_types = nsize;
  ia_params = new_params;
}

void make_particle_type_exist(int type)
{
  int ns = type + 1;
  if (ns <= n_particle_types)
    return;

  mpi_bcast_n_particle_types(ns);
}

void make_bond_type_exist(int type)
{
  int i, ns = type + 1;
  if(ns <= n_bonded_ia)
    return;
  bonded_ia_params = (Bonded_ia_parameters *)realloc(bonded_ia_params,
						     ns*sizeof(Bonded_ia_parameters));
  for (i = n_bonded_ia; i < ns; i++)
    bonded_ia_params[i].type = BONDED_IA_NONE;

  n_bonded_ia = ns;
}

int printBondedIAToResult(Tcl_Interp *interp, int i)
{
  Bonded_ia_parameters *params = &bonded_ia_params[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];

  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (params->type) {
  case BONDED_IA_FENE:
    Tcl_PrintDouble(interp, params->p.fene.k_fene, buffer);
    Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.fene.r_fene, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_ANGLE:
    Tcl_PrintDouble(interp, params->p.angle.bend, buffer);
    Tcl_AppendResult(interp, "angle ", buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_DIHEDRAL:
    Tcl_AppendResult(interp, "dihedral",(char *) NULL);
    return (TCL_OK);
  case BONDED_IA_NONE:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
		     (char *) NULL);
    return (TCL_ERROR);
  default:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction type",(char *) NULL);
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

int printNonbondedIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  if (!data) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "interaction does not exist",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  sprintf(buffer, "%d %d ", i, j);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  if (data->LJ_eps != 0) {
    Tcl_PrintDouble(interp, data->LJ_eps, buffer);
    Tcl_AppendResult(interp, "lennard-jones ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_offset, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_capradius, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  }
  if (data->ramp_cut > 0) {
    Tcl_PrintDouble(interp, data->ramp_cut, buffer);
    Tcl_AppendResult(interp, "ramp ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->ramp_force, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
  }
  return (TCL_OK);
}

void calc_maximal_cutoff()
{
  int i, j;
  max_cut = -1.0;
  /* bonded */
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      if(max_cut < bonded_ia_params[i].p.fene.r_fene)
	max_cut = bonded_ia_params[i].p.fene.r_fene;
      break;
    default:
      break;
    }
  }
  /* non bonded */
  for (i = 0; i < n_particle_types; i++)
     for (j = i; j < n_particle_types; j++) {
       	if (checkIfParticlesInteract(i, j)) {
	  IA_parameters *data = get_ia_param(i, j);
	  if (data->LJ_eps != 0) {
	    if(max_cut < (data->LJ_cut+data->LJ_offset) ) 
	      max_cut = (data->LJ_cut+data->LJ_offset);
	  }
	  if(max_cut < data->ramp_cut)
	    max_cut = data->ramp_cut;
	}
     }
  /* real cpace electrostatic */
  if(p3m.bjerrum > 0.0) {
    if(max_cut < p3m.r_cut)
      max_cut = p3m.r_cut;
  }

}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j;

  /* no argument -> print all interaction informations. */
  if (argc == 1) {
    int start = 1;
    for (i = 0; i < n_bonded_ia; i++) {
      if (bonded_ia_params[i].type != BONDED_IA_NONE) {
	if (start) {
	  Tcl_AppendResult(interp, "{", (char *)NULL);
	  start = 0;
	}
	else
	  Tcl_AppendResult(interp, " {", (char *)NULL);
	printBondedIAToResult(interp, i);
	Tcl_AppendResult(interp, "}", (char *)NULL);
      }
    }
    for (i = 0; i < n_particle_types; i++)
      for (j = i; j < n_particle_types; j++) {
	if (checkIfParticlesInteract(i, j)) {
	  if (start) {
	    Tcl_AppendResult(interp, "{", (char *)NULL);
	    start = 0;
	  }
	  else
	    Tcl_AppendResult(interp, " {", (char *)NULL);
	  printNonbondedIAToResult(interp, i, j);
	  Tcl_AppendResult(interp, "}", (char *)NULL);
	}
      }
    return (TCL_OK);
  }

  /* first argument should be an integer */
  if (Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) 
    return (TCL_ERROR);

  /* if second argument is not an integer -> bonded interaction.  */
  if (argc==2 || Tcl_GetInt(interp, argv[2], &j) == TCL_ERROR) {
    Tcl_ResetResult(interp);

    if (argc < 2) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " <type>  ?interaction? ?values?\"",
		       (char *) NULL);
      return (TCL_ERROR);
    }
    if(i < 0) {
      Tcl_AppendResult(interp, "interaction type must be nonnegative",
		       (char *) NULL);
      return (TCL_ERROR);
    }
    make_bond_type_exist(i);

    if (argc == 2) {
      char buffer[TCL_INTEGER_SPACE];
      /* print specific interaction information */
      if(i<n_bonded_ia) {
	printBondedIAToResult(interp, i);
	return (TCL_OK);
      }
      else {
	sprintf(buffer, "%d", i);
	Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
			 (char *) NULL);
	return (TCL_ERROR);
      }
    }

    if (argc == 3 && !strncmp(argv[2], "num", strlen(argv[2]))) {
      char buffer[TCL_INTEGER_SPACE];
      /* print interaction partner number */
      if(i<n_bonded_ia) {
	Bonded_ia_parameters *params = &bonded_ia_params[i];
	sprintf(buffer, "%d", params->num);
	Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
			 (char *) NULL);
	return (TCL_OK);
      }
      else {
	sprintf(buffer, "%d", i);
	Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
			 (char *) NULL);
	return (TCL_ERROR);
      }
    }

    /* set interaction parameters */

    argc -= 2;
    argv += 2;

    if (!strncmp(argv[0], "FENE", strlen(argv[0])) || 
	!strncmp(argv[0], "fene", strlen(argv[0]))) {
      /* set new FENE interaction type */
      if (argc != 3) {
	Tcl_AppendResult(interp, "fene needs 2 parameters: "
			 "<k_fene> <r_fene>", (char *) NULL);
	return (TCL_ERROR);
      }
      /* copy FENE parameters */
      if ((Tcl_GetDouble(interp, argv[1], &(bonded_ia_params[i].p.fene.k_fene)) == 
	   TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[2], &(bonded_ia_params[i].p.fene.r_fene))  == 
	   TCL_ERROR) ) 
	return (TCL_ERROR);
      bonded_ia_params[i].type = BONDED_IA_FENE;
      bonded_ia_params[i].num  = 1;
      /* broadcast interaction parameters */
      mpi_bcast_ia_params(i,-1); 
    }
    else if (!strncmp(argv[0], "angle", strlen(argv[0]))) {
      /* set new ANGLE interaction type */
      if (argc != 2) {
	Tcl_AppendResult(interp, "angle needs 1 parameter: "
			 "<bend> ", (char *) NULL);
	return (TCL_ERROR);
      }
      /* copy ANGLE parameters */
      if ((Tcl_GetDouble(interp, argv[1], &(bonded_ia_params[i].p.angle.bend)) == 
	   TCL_ERROR) )
	return (TCL_ERROR);
      bonded_ia_params[i].type = BONDED_IA_ANGLE;
      bonded_ia_params[i].num = 2;
      /* broadcast interaction parameters */
      mpi_bcast_ia_params(i,-1); 
    }
    else if (!strncmp(argv[0], "dihedral", strlen(argv[0]))) {
 
      /* remove this lines after the implementation of dihedral */
      Tcl_AppendResult(interp, "Do not use interaction type \"", argv[0],
		       "\"", (char *) NULL);
      return (TCL_ERROR);

      bonded_ia_params[i].type = BONDED_IA_DIHEDRAL;
      bonded_ia_params[i].num = 0;
      mpi_bcast_ia_params(i,-1); 
      return (TCL_OK);
    }
    else {
      Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		       "\"", (char *) NULL);
      return (TCL_ERROR);
    }
  }

  /* second argument is an integer -> non bonded interaction. */
  else {
    IA_parameters *data, *data_sym;

    if (argc < 3) {
      Tcl_AppendResult(interp, "wrong # args:  should be \"",
		       argv[0], " <type 1> <type 2> ?interaction? ?values?\"",
		       (char *) NULL);
      return (TCL_ERROR);
    }

    make_particle_type_exist(i);
    make_particle_type_exist(j);
    
    data     = get_ia_param(i, j);
    data_sym = get_ia_param(j, i);

    if (!data || !data_sym) {
      Tcl_AppendResult(interp, "particle types must be nonnegative",
		       (char *) NULL);
      return (TCL_ERROR);
    }
    /* print specific interaction information */
    if (argc == 3) {
      return printNonbondedIAToResult(interp, i, j);
    }    
    /* set interaction parameters */
    argc -= 3;
    argv += 3;
    
    while (argc > 0) {
      if (!strncmp(argv[0], "lennard-jones", strlen(argv[0]))) {
	/* set new lennard-jones interaction type */
	if (argc < 6) {
	  Tcl_AppendResult(interp, "lennard-jones needs 5 parameters: "
			   "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
			   (char *) NULL);
	return (TCL_ERROR);
	}
	/* copy lennard-jones parameters */
	if ((Tcl_GetDouble(interp, argv[1], &data->LJ_eps) == TCL_ERROR) ||
	    (Tcl_GetDouble(interp, argv[2], &data->LJ_sig)  == TCL_ERROR) ||
	    (Tcl_GetDouble(interp, argv[3], &data->LJ_cut)  == TCL_ERROR) ||
	    (Tcl_GetDouble(interp, argv[4], &data->LJ_shift)   == TCL_ERROR) ||
	    (Tcl_GetDouble(interp, argv[5], &data->LJ_offset)  == TCL_ERROR))
	  return (TCL_ERROR);
	
	/* LJ should be symmetrically */
	data_sym->LJ_eps = data->LJ_eps;
	data_sym->LJ_sig = data->LJ_sig;
	data_sym->LJ_cut = data->LJ_cut;
	data_sym->LJ_shift   = data->LJ_shift;
	data_sym->LJ_offset  = data->LJ_offset;
	argc -= 6;
	argv += 6;
	/* broadcast interaction parameters */
	mpi_bcast_ia_params(i, j);
	mpi_bcast_ia_params(j, i);
      }
      else if (!strncmp(argv[0], "ramp", strlen(argv[0]))) {
	/* set new ramp interaction type */
	if (argc < 3) {
	  Tcl_AppendResult(interp, "ramp potential needs 2 parameters: "
			   "<ramp_cut> <ramp_force>",
			   (char *) NULL);
	  return (TCL_ERROR);
	}
	
	if ((Tcl_GetDouble(interp, argv[1], &data->ramp_cut) == TCL_ERROR) ||
	    (Tcl_GetDouble(interp, argv[2], &data->ramp_force)  == TCL_ERROR))
	  return (TCL_ERROR);
	
	/* ramp should be symmetrically */
	data_sym->ramp_cut = data->ramp_cut;
	data_sym->ramp_force = data->ramp_force;
	argc -= 3;
	argv += 3;
	/* broadcast interaction parameters */	
	mpi_bcast_ia_params(i, j);
	mpi_bcast_ia_params(j, i);
      }
      else {
	Tcl_AppendResult(interp, "unknown interaction type \"", argv[3],
			 "\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }
  return (TCL_OK);
}

int lj_force_cap_callback(Tcl_Interp *interp, void *data)
{
  lj_force_cap = *(double *)data;

  mpi_bcast_parameter(FIELD_LJFORCECAP);
  mpi_bcast_event(INTERACTION_CHANGED);
  return (TCL_OK);
}
