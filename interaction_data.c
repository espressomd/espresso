/** \file interaction_data.c
    Implementation of \ref interaction_data.h "interaction_data.h"
 */
#include <string.h>
#include <stdlib.h>
#include "config.h"
#include "debug.h"
#include "interaction_data.h"
#include "communication.h"
#include "grid.h"
#include "p3m.h"
#include "debye_hueckel.h"
#include "lj.h"

/****************************************
 * variables
 *****************************************/
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

Coulomb_parameters coulomb = { 0.0, COULOMB_NONE };

Debye_hueckel_params dh_params = { 0.0, 0.0, 0.0, 0.0 };

int n_bonded_ia = 0;
Bonded_ia_parameters *bonded_ia_params = NULL;

double max_cut;

double lj_force_cap = 0.0;

#ifdef CONSTRAINTS
int n_constraints       = 0;
Constraint *constraints = NULL;
#endif 

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
    Tcl_PrintDouble(interp, params->p.fene.k, buffer);
    Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.fene.r, buffer);
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

int printCoulombIAToResult(Tcl_Interp *interp) 
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  if (coulomb.bjerrum == 0.0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "coulomb 0.0",
		     (char *) NULL);
    return (TCL_OK);
  }
  Tcl_PrintDouble(interp, coulomb.bjerrum, buffer);
  Tcl_AppendResult(interp, "coulomb ", buffer, " ", (char *) NULL);
  if (coulomb.method == COULOMB_P3M) {
    Tcl_PrintDouble(interp, p3m.r_cut, buffer);
    Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.mesh[0]);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.cao);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.alpha, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.accuracy, buffer);
    Tcl_AppendResult(interp, buffer, "} ", (char *) NULL);

    Tcl_PrintDouble(interp, p3m.epsilon, buffer);
    Tcl_AppendResult(interp, "{coulomb epsilon ", buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.inter);
    Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[0], buffer);
    Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
  }
  else if (coulomb.method == COULOMB_DH) {
    Tcl_PrintDouble(interp, dh_params.kappa, buffer);
    Tcl_AppendResult(interp, "dh ", buffer, " ",(char *) NULL);
    Tcl_PrintDouble(interp, dh_params.r_cut, buffer);
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
      if(max_cut < bonded_ia_params[i].p.fene.r)
	max_cut = bonded_ia_params[i].p.fene.r;
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
  /* real space electrostatic */
  if(coulomb.method == COULOMB_P3M     && max_cut < p3m.r_cut) 
    max_cut = p3m.r_cut;
  else if(coulomb.method == COULOMB_DH && max_cut < dh_params.r_cut) 
    max_cut = dh_params.r_cut;
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
    if(coulomb.bjerrum > 0.0) {
      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);
      printCoulombIAToResult(interp);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
    if(lj_force_cap != 0.0) {
      char buffer[TCL_DOUBLE_SPACE];

      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);
      if (lj_force_cap == -1.0)
	Tcl_AppendResult(interp, "ljforcecap individual");
      else {
	Tcl_PrintDouble(interp, lj_force_cap, buffer);
	Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
      }
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
    return (TCL_OK);
  }

  if(!strncmp(argv[1], "ljforcecap", strlen(argv[1]))) {
    if (argc == 2) {
      char buffer[TCL_DOUBLE_SPACE];
      if (lj_force_cap == -1.0)
	Tcl_AppendResult(interp, "ljforcecap individual", (char *) NULL);
      else {
	Tcl_PrintDouble(interp, lj_force_cap, buffer);
	Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
      }
      return TCL_OK;
    }

    if (argc > 3) {
      Tcl_AppendResult(interp, "inter ljforcecap takes at most 1 parameter", (char *) NULL);      
      return TCL_ERROR;
    }

    if (!strncmp(argv[2], "individual", strlen(argv[2])))
      lj_force_cap = -1.0;
    else if (Tcl_GetDouble(interp, argv[2], &lj_force_cap) == TCL_ERROR || lj_force_cap < 0) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		       (char *) NULL);
      return (TCL_ERROR);  
    }

    if (lj_force_cap != -1.0)
      mpi_lj_cap_forces(lj_force_cap);
    return (TCL_OK);
  }

  if(!strncmp(argv[1], "coulomb", strlen(argv[1]))) {
    /* print coulomb interaction parameters */
    if(argc < 3) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      printCoulombIAToResult(interp);
      Tcl_AppendResult(interp, "}", (char *)NULL);
      return (TCL_OK);
    }

    /* check bjerrum length */
    if ((Tcl_GetDouble(interp, argv[2], &(coulomb.bjerrum)) == TCL_ERROR)) {
      Tcl_ResetResult(interp);
      /* if p3m is set already you can set additional optional p3m parameters */
      if (coulomb.method == COULOMB_P3M) {
	argc -= 2;
	argv += 2;
	while (argc > 0) {
	  /* p3m parameter: inter */
	  if (!strncmp(argv[0], "n_interpol", strlen(argv[0]))) {
	    if(argc<2) return (TCL_ERROR);
	    if ((Tcl_GetInt(interp, argv[1], &(p3m.inter)) == TCL_ERROR)) 
	      return (TCL_ERROR);
	    if (p3m.inter < 0) return (TCL_ERROR);
	    argc -= 2;
	    argv += 2;
	  }
	  /* p3m parameter: mesh_off */
	  else if (!strncmp(argv[0], "mesh_off", strlen(argv[0]))) {
	    if(argc<4) return (TCL_ERROR);
	    if ((Tcl_GetDouble(interp, argv[1], &(p3m.mesh_off[0])) == TCL_ERROR) ||
		(Tcl_GetDouble(interp, argv[2], &(p3m.mesh_off[1])) == TCL_ERROR) ||
		(Tcl_GetDouble(interp, argv[3], &(p3m.mesh_off[2])) == TCL_ERROR) ) 
	      return (TCL_ERROR);
	    if(p3m.mesh_off[0] < 0.0 || p3m.mesh_off[0] > 1.0 ||
	       p3m.mesh_off[1] < 0.0 || p3m.mesh_off[1] > 1.0 ||
	       p3m.mesh_off[2] < 0.0 || p3m.mesh_off[2] > 1.0 )
	      return (TCL_ERROR);
	    argc -= 4;
	    argv += 4;
	  }
	  /* p3m parameter: epsilon */
	  else if(!strncmp(argv[0], "epsilon", strlen(argv[0]))) {
	    if(argc<2) return (TCL_ERROR);
	    if ((Tcl_GetDouble(interp, argv[1], &(p3m.epsilon)) == TCL_ERROR))
	      return (TCL_ERROR);
	    argc -= 2;
	    argv += 2;	    
	  }
	  else {
	    Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
	    return (TCL_ERROR);
	  }
	}
      }
      else {
	Tcl_AppendResult(interp, "Expect: inter coulomb <bjerrum> || <p3m-parameter>",(char *) NULL);
	return (TCL_ERROR);
      }
      mpi_bcast_coulomb_params();
      return (TCL_OK);
    }
    
    /* check bjerrum length */
    if(coulomb.bjerrum < 0.0) {
      Tcl_AppendResult(interp, "bjerrum length must be positiv",(char *) NULL);
      return (TCL_ERROR);
    }
    argc -= 3;
    argv += 3;

    if(coulomb.bjerrum == 0.0) {
      coulomb.method = COULOMB_NONE;
      mpi_bcast_coulomb_params();
      return (TCL_OK);
    }

    /* check number of parameters */
    if(argc < 2) {
      Tcl_AppendResult(interp, "wrong # args for inter coulomb.", (char *) NULL);
      return (TCL_ERROR);
    }

    /* check method */
    if(!strncmp(argv[0], "p3m", strlen(argv[0])) || 
       !strncmp(argv[0], "P3M", strlen(argv[0])) ) {
      argc -= 1;
      argv += 1;
      
      coulomb.method = COULOMB_P3M;
      p3m.bjerrum    = coulomb.bjerrum;
      
#ifdef PARTIAL_PERIODIC
      if(periodic[0]==0 || periodic[1]==0 || periodic[2]==0) {
	Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb P3M", (char *) NULL);
	return (TCL_ERROR);  
      }
#endif
      
      if(Tcl_GetDouble(interp, argv[0], &(p3m.r_cut)) == TCL_ERROR) {
	Tcl_ResetResult(interp);
	/* must be tune, tune parameters */
	if(strncmp(argv[0], "tune", strlen(argv[0]))) {
	  Tcl_AppendResult(interp, "Unknown p3m parameter \"",argv[0],"\"",(char *) NULL);
	  return (TCL_ERROR);  
	}
	argc -= 1;
	argv += 1;
	/* parse tune parameters */
	while(argc > 0) {
	  if(!strncmp(argv[0], "r_cut", strlen(argv[0]))) {
	    if(argc > 1) {
	      if(Tcl_GetDouble(interp, argv[1], &(p3m.r_cut)) == TCL_ERROR) 
		return (TCL_ERROR); 
	    }
	    else return (TCL_ERROR); 
	  }
	  else if(!strncmp(argv[0], "mesh", strlen(argv[0]))) {
	    if(argc > 1) { 
	      if(Tcl_GetInt(interp, argv[1], &(p3m.mesh[0])) == TCL_ERROR) 
		return (TCL_ERROR);
	    } 
	    else return (TCL_ERROR); 
	  }
	  else if(!strncmp(argv[0], "cao", strlen(argv[0]))) {
	    if(argc > 1) {
	      if(Tcl_GetInt(interp, argv[1], &(p3m.cao)) == TCL_ERROR) 
		return (TCL_ERROR);
	    } 
	    else return (TCL_ERROR); 
	  }
	  else if(!strncmp(argv[0], "accuracy", strlen(argv[0]))) {
	    if(argc > 1) { 
	      if(Tcl_GetDouble(interp, argv[1], &(p3m.accuracy)) == TCL_ERROR) 
		return (TCL_ERROR);
	    } 
	    else return (TCL_ERROR); 
	  }
	  else {
	    Tcl_AppendResult(interp, "Unkwon p3m tune parameter \"",argv[0],"\"",(char *) NULL);
	    return (TCL_ERROR);  
	  }
	  argc-=2;
	  argv+=2;
	}
	/* check tune parameters */
	if(p3m.accuracy <= 0.0) {
	  Tcl_AppendResult(interp, "You have to give p3m tune a positive value for accuracy!",(char *) NULL);
	  return (TCL_ERROR);
	}
	if(p3m.r_cut < 0.0) {
	  Tcl_AppendResult(interp, "p3m tune r_cut must be positiv.",(char *) NULL);
	  return (TCL_ERROR);
	}
	if(p3m.mesh[0] < 0) {
	  Tcl_AppendResult(interp, "p3m tune mesh must be positiv.",(char *) NULL);
	  return (TCL_ERROR);
	}
	p3m.mesh[2] = p3m.mesh[1] = p3m.mesh[0];
	if(p3m.cao < 0 || p3m.cao > 7) {
	  Tcl_AppendResult(interp, "p3m tune cao must be in between 1 and 7",(char *) NULL);
	  return (TCL_ERROR);
	}
	if(P3M_tune_parameters(interp) == TCL_ERROR) 
	  return (TCL_ERROR);
      }
      else {
	/* parse p3m parameters */
	if(argc<3) {
	  Tcl_AppendResult(interp, "p3m needs at least 3 parameters: <r_cut> <mesh> <cao> [alpha] [accuracy]",(char *) NULL);
	  return (TCL_ERROR);  
	}
	if(Tcl_GetInt(interp, argv[1], &(p3m.mesh[0])) == TCL_ERROR ||
	   Tcl_GetInt(interp, argv[2], &(p3m.cao    )) == TCL_ERROR ) 
	  return (TCL_ERROR);
	
	if(p3m.mesh[0] < 0) {
	  Tcl_AppendResult(interp, "p3m mesh must be positiv.",(char *) NULL);
	  return (TCL_ERROR);
	}
	p3m.mesh[2] = p3m.mesh[1] = p3m.mesh[0];
	if(p3m.cao < 1 || p3m.cao > 7) {
	  Tcl_AppendResult(interp, "p3m cao must be between 1 and 7.",(char *) NULL);
	  return (TCL_ERROR);
	}
	if(p3m.cao > p3m.mesh[0] )  {
	  Tcl_AppendResult(interp, "p3m cao can not be larger than p3m mesh",(char *) NULL);
	  return (TCL_ERROR);
	}
	if(argc>3) {
	  if(Tcl_GetDouble(interp, argv[3], &(p3m.alpha)) == TCL_ERROR )
	    return (TCL_ERROR);
	}
	else {
	  Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",(char *) NULL);
	  return (TCL_ERROR);  
	}
	if(argc>4) {
	  if(Tcl_GetDouble(interp, argv[4], &(p3m.accuracy)) == TCL_ERROR )
	    return (TCL_ERROR);
	}
      }
    }
    else if (!strncmp(argv[0], "dh ", strlen(argv[0])) ||
	     !strncmp(argv[0], "DH ", strlen(argv[0])) ) {
      argc -= 1;
      argv += 1;
      
      coulomb.method = COULOMB_DH;
      dh_params.bjerrum = coulomb.bjerrum;

      if(argc < 2) {
	Tcl_AppendResult(interp, "Not enough parameters: inter coulomb dh <kappa> <r_cut>", (char *) NULL);
	return (TCL_ERROR);
      }
      if(Tcl_GetDouble(interp, argv[0], &(dh_params.kappa)) == TCL_ERROR ||
	 Tcl_GetDouble(interp, argv[1], &(dh_params.r_cut)) == TCL_ERROR) 
	return (TCL_ERROR);
	
      if(dh_params.kappa < 0.0) {
	Tcl_AppendResult(interp, "dh kappa must be positiv.",(char *) NULL);
	return (TCL_ERROR);
      }
      if(dh_params.r_cut < 0.0) {
	Tcl_AppendResult(interp, "dh r_cut must be positiv.",(char *) NULL);
	return (TCL_ERROR);
      }
    }
    else {
      coulomb.bjerrum = 0.0;
      coulomb.method  = COULOMB_NONE;
      /* communicate bjerrum length */
      Tcl_AppendResult(interp, "Do not know coulomb method \"",argv[1],"\": coulomb switched off",(char *) NULL);
      return (TCL_ERROR);
    }
    mpi_bcast_coulomb_params();
    return (TCL_OK);
  }

  if (Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "don't know about \"",argv[1],"\"",(char *) NULL);
    return (TCL_ERROR);
  }
  
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
      if ((Tcl_GetDouble(interp, argv[1], &(bonded_ia_params[i].p.fene.k)) == 
	   TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[2], &(bonded_ia_params[i].p.fene.r))  == 
	   TCL_ERROR) ) 
	return (TCL_ERROR);
      bonded_ia_params[i].type = BONDED_IA_FENE;
      bonded_ia_params[i].p.fene.r2 = SQR(bonded_ia_params[i].p.fene.r);
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
	if (argc >= 1 && (Tcl_GetDouble(interp, argv[0], &data->LJ_capradius) == TCL_OK)) {
	  data_sym->LJ_capradius = data->LJ_capradius;
	  argc--;
	  argv++;
	}
	/* broadcast interaction parameters */
	mpi_bcast_ia_params(i, j);
	mpi_bcast_ia_params(j, i);
	if (lj_force_cap != -1.0)
	  mpi_lj_cap_forces(lj_force_cap);
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
	Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
			 "\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }
  return (TCL_OK);
}

#ifdef CONSTRAINTS
int printConstraintToResult(Tcl_Interp *interp, int i)
{
  Constraint con = constraints[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];

  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (con.type) {
  case CONSTRAINT_WAL:
    Tcl_PrintDouble(interp, con.c.wal.n[0], buffer);
    Tcl_AppendResult(interp, "wall ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.wal.n[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.wal.n[2], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.wal.d, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    break;
  case CONSTRAINT_SPH:
    Tcl_PrintDouble(interp, con.c.sph.pos[0], buffer);
    Tcl_AppendResult(interp, "sphere ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.sph.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.sph.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con.c.sph.rad, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    break;
  case CONSTRAINT_CYL:
    Tcl_AppendResult(interp, "cylinder not implemented", (char *) NULL);
    return (TCL_OK);
  default:
    Tcl_AppendResult(interp, "unknown constraint type", (char *) NULL);
    return (TCL_OK);
  }

  Tcl_PrintDouble(interp, con.LJ_eps, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, con.LJ_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, con.LJ_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, con.LJ_shift, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return (TCL_OK);
 
}
#endif

int constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef CONSTRAINTS
  int i, del_num = -1;
  double norm;
  Constraint *con=NULL;

  /* no argument -> print constraint information */
  if (argc == 1) {
    if(n_constraints>0) Tcl_AppendResult(interp, "{", (char *)NULL);
    for (i = 0; i < n_constraints; i++) {
      if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
      printConstraintToResult(interp, i);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
    return (TCL_OK);
  }
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    if(argc < 10) {
      Tcl_AppendResult(interp, "wrong # args! Usage: constraint wall <nx> <ny> <nz> <d> <lj_eps> <lj_sig> <lj_cut> <lj_shift>", 
		       (char *) NULL);
      return (TCL_ERROR);
    }
    n_constraints++;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
    con = &constraints[n_constraints-1];
    con->type = CONSTRAINT_WAL;
    if(Tcl_GetDouble(interp, argv[2], &(con->c.wal.n[0])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[3], &(con->c.wal.n[1])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[4], &(con->c.wal.n[2])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[5], &(con->c.wal.d)) == TCL_ERROR)
      return (TCL_ERROR);
    norm = SQR(con->c.wal.n[0])+SQR(con->c.wal.n[1])+SQR(con->c.wal.n[2]);
    for(i=0;i<3;i++) con->c.wal.n[i] /= sqrt(norm);
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    if(argc < 10) {
      Tcl_AppendResult(interp, "wrong # args! Usage: constraint sphere <px> <py> <pz> <rad> <lj_eps> <lj_sig> <lj_cut> <lj_shift>", (char *) NULL);
      return (TCL_ERROR);
    }
    n_constraints++;
    constraints = realloc(constraints,n_constraints*sizeof(Constraint));
    con = &constraints[n_constraints-1];
    con->type = CONSTRAINT_SPH;
     if(Tcl_GetDouble(interp, argv[2], &(con->c.sph.pos[0])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[3], &(con->c.sph.pos[1])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[4], &(con->c.sph.pos[2])) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[5], &(con->c.sph.rad)) == TCL_ERROR)
      return (TCL_ERROR);
    if(con->c.sph.rad < 0.0) return (TCL_ERROR);
   }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all constraints */
      n_constraints = 0;
      constraints = realloc(constraints,n_constraints*sizeof(Constraint));
    }
    if(Tcl_GetInt(interp, argv[2], &(del_num)) == TCL_ERROR) return (TCL_ERROR);
    if(del_num >= n_constraints) {
      Tcl_AppendResult(interp, "Can not delete non existing constraint",(char *) NULL);
      return (TCL_OK);
    }
  }
  else {
    Tcl_AppendResult(interp, "unknown feature of constraints! Usage:\n",(char *) NULL);
    Tcl_AppendResult(interp, "constraint wall <nx> <ny> <nz> <d> <lj_eps> <lj_sig> <lj_cut> <lj_shift>\n",(char *) NULL);
    Tcl_AppendResult(interp, "constraint sphere <px> <py> <pz> <rad> <lj_eps> <lj_sig> <lj_cut> <lj_shift>\n",(char *) NULL);
    Tcl_AppendResult(interp, "constraint delete <num>",(char *) NULL);
    return (TCL_ERROR);
  }

  /* set the LJ parameters */
  if(del_num<0) {
    if(Tcl_GetDouble(interp, argv[6], &(con->LJ_eps)) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[7], &(con->LJ_sig)) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[8], &(con->LJ_cut)) == TCL_ERROR ||
       Tcl_GetDouble(interp, argv[9], &(con->LJ_shift)) == TCL_ERROR)
    return (TCL_ERROR);
    if(con->LJ_eps<0.0 || con->LJ_sig < 0.0 || con->LJ_cut < 0.0) {
      Tcl_AppendResult(interp, "Wrong constraint LJ parameters",(char *) NULL);
      return (TCL_ERROR);
    }
  }


  mpi_bcast_constraint(del_num);
  return (TCL_OK);
#else
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}

