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
#ifndef HERTZIAN_H
#define HERTZIAN_H

/** \file hertzian.h
 *  Routines to calculate the Hertzian energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef HERTZIAN
#include "mol_cut.h"

MDINLINE int tclprint_to_result_HertzianIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->Hertzian_eps, buffer);
  Tcl_AppendResult(interp, "hertzian ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->Hertzian_sig, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

MDINLINE int hertzian_set_params(int part_type_a, int part_type_b,
				 double eps, double sig)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* Hertzian should be symmetrically */
  data->Hertzian_eps = data_sym->Hertzian_eps = eps;
  data->Hertzian_sig = data_sym->Hertzian_sig = sig;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int tclcommand_inter_parse_hertzian(Tcl_Interp * interp,
			     int part_type_a, int part_type_b,
			     int argc, char ** argv)
{
  /* parameters needed for Hertzian */
  double eps, sig;
  int change;

  /* get interaction type */
  if (argc < 3) {
    Tcl_AppendResult(interp, "Hertzian potential needs 2 parameters: "
		     "<epsilon> <sigma>", (char *) NULL);
    return 0;
  }

  /* copy parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))) {
    Tcl_AppendResult(interp, "Hertzian potential needs 2 parameters: "
		     "<epsilon> <sigma>", (char *) NULL);
    return TCL_ERROR;
  }
  change = 3;
  
  Tcl_ResetResult(interp);
  if (hertzian_set_params(part_type_a, part_type_b,
			  eps, sig) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}


/** Calculate Hertzian force between particle p1 and p2 */
MDINLINE void add_hertzian_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				      double d[3], double dist, double dist2, double force[3])
{
  double fac;
  int j;
  if (CUTOFF_CHECK(dist < ia_params->Hertzian_sig)) {
    fac = 5./2.*ia_params->Hertzian_eps/ia_params->Hertzian_sig *
      pow(1 - dist/ia_params->Hertzian_sig, 3./2.)/dist;

    for(j=0;j<3;j++)
      force[j] += fac * d[j];
  }
}

/** calculate Lennard jones energy between particle p1 and p2. */
MDINLINE double hertzian_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				     double d[3], double dist, double dist2)
{
  if (CUTOFF_CHECK(dist < ia_params->Hertzian_sig)) {
    return ia_params->Hertzian_eps *
      pow(1 - dist/ia_params->Hertzian_sig, 5./2.);
  }
  return 0.0;
}

#endif /* ifdef HERTZIAN */
#endif
