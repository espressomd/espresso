#include "interaction_data.h"
#include "global.h"
#include "communication.h"
#include <string.h>
#include <stdlib.h>

IA_parameters *safe_get_ia_param(int i, int j)
{
  if ((i < 0) || (j < 0))
    return NULL;

  /* expand array if necessary */
  realloc_ia_params(((i > j) ? i : j) + 1);

  return &ia_params[i*n_particle_types + j];
}

void realloc_ia_params(int nsize)
{
  int i, j;
  IA_parameters *new_params;

  if (nsize <= n_particle_types)
    return;

  if (ia_params)
    free(ia_params);
  new_params = (IA_parameters *) malloc(nsize*nsize*sizeof(IA_parameters));
  for (i = 0; i < nsize; i++)
    for (j = 0; j < nsize; j++) {
      if ((i < n_particle_types) && (i < n_particle_types))
	copy_ia_params(&new_params[i*n_particle_types + j],
			   &ia_params[i*n_particle_types + j]);
      else
	initialize_ia_params(&new_params[i*n_particle_types + j]);
    }

  n_particle_types = nsize;
  ia_params = new_params;
}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j;
  IA_parameters *data, *data_sym;

  if (argc < 3) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) ||
      (Tcl_GetInt(interp, argv[2], &j) == TCL_ERROR))
    return (TCL_ERROR);

  data     = safe_get_ia_param(i, j);
  data_sym = safe_get_ia_param(j, i);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (argc == 3) {
    /* print interaction information */
    char buffer[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, data->LJ_eps, buffer);
    Tcl_AppendResult(interp, "{lennard-jones ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_offset, buffer);
    Tcl_AppendResult(interp, buffer, "} ", (char *) NULL);

    Tcl_PrintDouble(interp, data->ramp_cut, buffer);
    Tcl_AppendResult(interp, "{ramp ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->ramp_force, buffer);
    Tcl_AppendResult(interp, buffer, "}", (char *) NULL);
    return (TCL_OK);
  }

  /* set interaction parameters */
  argc -= 3;
  argv += 3;

  while (argc > 0) {
    if (!strncmp(argv[0], "lennard-jones", strlen(argv[0]))) {
      if (argc < 6) {
	Tcl_AppendResult(interp, "lennard-jones needs 5 parameters: "
			 "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      
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

      mpi_bcast_ia_params(i, j);
      mpi_bcast_ia_params(j, i);
    }
    else if (!strncmp(argv[0], "ramp", strlen(argv[0]))) {
      if (argc < 3) {
	Tcl_AppendResult(interp, "lennard-jones needs 2 parameters: "
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

      mpi_bcast_ia_params(i, j);
      mpi_bcast_ia_params(j, i);
    }
    else {
      Tcl_AppendResult(interp, "unknown interaction type \"", argv[3],
		       "\"", (char *)NULL);
      return (TCL_ERROR);
    }
  }

  return (TCL_OK);
}

int niatypes_callback(Tcl_Interp *interp, void *data)
{
  n_interaction_types = *(int *)data;
  return (TCL_OK);
}
