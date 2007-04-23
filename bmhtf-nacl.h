// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef BMHTF_NACL_H
#define BMHTF_NACL_H

/** \file bmhtf-nacl.h
 *  Routines to calculate the Born-Meyer-Huggins-Tosi-Fumi energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/

#ifdef BMHTF_NACL

MDINLINE int printBMHTFIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_AppendResult(interp, "bmhtf-nacl ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_A, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_B, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_C, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_D, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->BMHTF_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  
  return TCL_OK;
}

MDINLINE int BMHTF_set_params(int part_type_a, int part_type_b,
			      double A, double B, double C,
			      double D, double sig, double cut)
{
  IA_parameters *data, *data_sym;
  double shift, dist2, pw6;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* should be symmetrically */
  data->BMHTF_A   = data_sym->BMHTF_A   = A;
  data->BMHTF_B   = data_sym->BMHTF_B   = B;
  data->BMHTF_C   = data_sym->BMHTF_C   = C;
  data->BMHTF_D   = data_sym->BMHTF_D   = D;
  data->BMHTF_sig = data_sym->BMHTF_sig = sig;
  data->BMHTF_cut = data_sym->BMHTF_cut = cut;
  dist2 = cut*cut;
  pw6 = dist2*dist2*dist2;
  shift = -(A*exp(B*(sig - cut)) - C/pw6 - D/pw6/dist2);

  data->BMHTF_computed_shift =
    data_sym->BMHTF_computed_shift = shift;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int BMHTF_parser(Tcl_Interp * interp,
			  int part_type_a, int part_type_b,
			  int argc, char ** argv)
{
  double A, B, C, D, sig, cut;

  if (argc < 7) {
    Tcl_AppendResult(interp, "BMHTF NaCl potential needs 6 parameters: "
		     "<A> <B> <C> <D> <sigma> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy smooth step parameters */
  if ((! ARG_IS_D(1, A))    ||
      (! ARG_IS_D(2, B))    ||
      (! ARG_IS_D(3, C))    ||
      (! ARG_IS_D(4, D))    ||
      (! ARG_IS_D(5, sig))  ||
      (! ARG_IS_D(6, cut)   )) {
    Tcl_AppendResult(interp, "BMHTF NaCl potential needs 6 parameters: "
		     "<A> <B> <C> <D> <sigma> <cutoff>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (BMHTF_set_params(part_type_a, part_type_b,
		       A, B, C, D, sig, cut) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return 7;
}

/** Calculate smooth step force between particle p1 and p2 */
MDINLINE void add_BMHTF_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				   double d[3], double dist, double dist2, double force[3])
{
  int j;
  double pw8, fac = 0.0;
  if(dist < ia_params->BMHTF_cut) {
    pw8 = dist2*dist2*dist2*dist2;
    fac = ia_params->BMHTF_A*ia_params->BMHTF_B*
      exp(ia_params->BMHTF_B*(ia_params->BMHTF_sig - dist))/dist -
      6*ia_params->BMHTF_C/pw8 - 8*ia_params->BMHTF_D/pw8/dist2;

    for(j=0;j<3;j++) force[j] += fac * d[j];
  }
}

/** calculate smooth step potential energy between particle p1 and p2. */
MDINLINE double BMHTF_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				  double d[3], double dist, double dist2)
{
  double pw6;
 
  if(dist < ia_params->BMHTF_cut) {
    pw6 = dist2*dist2*dist2;
    return ia_params->BMHTF_A*
      exp(ia_params->BMHTF_B*(ia_params->BMHTF_sig - dist)) -
      ia_params->BMHTF_C/pw6 - ia_params->BMHTF_D/pw6/dist2 + ia_params->BMHTF_computed_shift;
  }
  return 0.0;
}

#endif
#endif
